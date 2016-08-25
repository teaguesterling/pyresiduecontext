from __future__ import division, print_function

import argparse
import functools
import gzip
import logging
import math
import os
import sys
from collections import namedtuple

import numpy as np
from Bio.PDB.PDBParser import PDBParser as _RealPDBParser
from Bio.PDB.Polypeptide import three_to_one as residue_three_code_to_one_code


def get_pdb_id_from_filename(filename):
    pdbid = os.path.basename(filename).split('.')[0]
    return pdbid


class PDBParser:
    def __init__(self):
        self._parser = _RealPDBParser()

    def open(self, source):
        if isinstance(source, str):
            if source.endswith('.gz'):
                ctx = gzip.open(source, 'r')
            else:
                ctx = open(source, 'r')
        else:
            ctx = source

        return ctx

    def guess_name(self, source, default="Unknown"):
        if isinstance(source, str):
            name = get_pdb_id_from_filename(source)
        elif hasattr(source, 'name'):
            name = get_pdb_id_from_filename(source.name)
        else:
            name = default

        return name

    def load_structure(self, source, name=None):
        if name is None:
            name = self.guess_name(source)
        with self.open(source) as f:
            structure = self._parser.get_structure(name, f)
        return structure


class ArgumentParser(argparse.ArgumentParser):
    source = None
    chain = None
    output = sys.stdout
    inverse_output = None

    bin_divisor = 6
    radial_bin_size = 1.0

    merge_low_bins = False

    def __init__(self, *args, **kwargs):
        stdout = kwargs.pop('stdout', sys.stdout)
        super(ArgumentParser, self).__init__(*args, **kwargs)

        self.add_argument('source', type=argparse.FileType('r'))
        self.add_argument('chain', type=str)
        self.add_argument('-o', '--output', type=argparse.FileType('w'), default=stdout)
        self.add_argument('-O', '--inverse-output', dest='inverse_output', type=argparse.FileType('w'), default=None)

        self.add_argument('-b', '--bindivisor', dest='bin_divisor', type=int, default=self.bin_divisor)
        self.add_argument('-a', '--rbinsize', dest='radial_bin_size', type=float, default=self.radial_bin_size)
        self.add_argument('-M', '--no-merge-low-bins', dest='merge_low_bins',
                          action='store_const', const=not self.merge_low_bins, default=self.merge_low_bins)


    def parse_args(self, args=None, namespace=None):
        super(ArgumentParser, self).parse_args(args[1:] if args is sys.argv else args, namespace=self)

        return self


class InvalidContextBuilderState(Exception):
    pass


class ContextBuilder:
    state = -1

    Context = namedtuple('ResidueContext', (
        'index',
        'name',
        'serial',
        'aacode',
        'size',
        'x',
        'y',
        'z',
        'histogram',
        'inverted',
    ))

    def __init__(self, params, logger=logging):
        self.params = params

        self.radial_bin_size = self.params.radial_bin_size
        self.bin_divisor = self.params.bin_divisor
        self.bin_angle = np.pi / self.bin_divisor

        self.one_radial_step_bins = 2 * self.bin_divisor * self.bin_divisor

        # Step between 1 and -1 by 2 * bin_divisor increments
        # Then compute the arc cosine at each point as a bin boundary
        self.inv_zenith_step = -2 / self.bin_divisor
        self.zenith_boundaries = np.arccos(np.arange(1, -1, self.inv_zenith_step))[:self.bin_divisor]

        self.logger = logger
        self._parser = PDBParser()
        self.state = 0

    def _quantize_radial_bin(self, r):
        if r < math.e:
            ln = 0
        else:
            ln = math.ceil(math.log(r))
        frac = ln / self.radial_bin_size
        quantized = int(frac)
        return quantized

    def _qunatize_radial_bins(self, rs):
        take = rs > np.e
        ln = np.zeros(rs.shape, dtype=np.float)
        ln[take] = np.ceil(np.log(rs[take]))
        if self.radial_bin_size == 1:
            frac = ln
        else:
            frac = ln / self.radial_bin_size
        quantized = frac.astype(np.int)
        return quantized

    def _generate_context_name(self, atom):
        residue = atom.parent
        res_id = residue.id[1]
        name = 'd{}{}.{:d}'.format(self.params.source.name, self.params.chain, res_id)
        return name

    def load_structure(self):
        self.logger.debug("1) Loading structure")
        self.structure = self._parser.load_structure(self.params.source)
        model = self.structure.get_list()[0]  # Only consider the first model
        self.chain = [chain for chain in model if chain.get_id() == self.params.chain][0]
        self.logger.info("Loaded {} from {} (chain: )".format(
            self.structure.id,
            self.params.source.name,
            self.chain.get_id())
        )
        self.state = 1
        return self.structure

    def extract_alpha_carbons(self):
        if self.state < 1:
            raise InvalidContextBuilderState("Cannot extract alpha carbons before structure is loaded")
        self.logger.debug("2) Extracting Alpha Carbons")

        self.alpha_carbons = []
        for res in self.chain.get_residues():
            for atom in res:
                if atom.get_name() == 'CA':
                    self.alpha_carbons.append(atom)

        self.num_residues = len(self.alpha_carbons)

        self.logger.info("Loaded {} alpha carbons".format(self.num_residues))
        self.state = 2
        return self.alpha_carbons

    def build_coordinate_matrix(self):
        if self.state < 2:
            raise InvalidContextBuilderState("Cannot create coordinate matrix without alpha carbons")
        self.logger.debug("3) Building coordinate matrix")

        self.coords = np.empty((self.num_residues, 3), dtype=np.float)
        self.coords[..., :] = [atom.get_coord() for atom in self.alpha_carbons]

        self.extents = np.array([
            self.coords.min(axis=0),
            self.coords.max(axis=0)
        ])

        self.state = 3
        return self.coords

    def calculate_displacement_matrices(self):
        if self.state < 3:
            raise InvalidContextBuilderState("Cannot calculate displacement matrices without coordinates")
        self.logger.debug("4) Calculating NxN displacements")

        # Create N copies of coordinate matrix and subtract off a different row from each
        # Done all at once for speed
        self.displacements = self.coords[..., np.newaxis, :] - np.tile(self.coords, (self.num_residues, 1, 1))

        self.state = 4
        return self.displacements

    def calculate_spherical_coordinates(self):
        if self.state < 4:
            raise InvalidContextBuilderState("Cannot calculate spherical coordinates without displacements")
        self.logger.debug("5) Converting to spherical coordinates")

        self.spherical_coords = np.empty((self.num_residues,) + self.coords.shape, dtype=np.float)
        x, y, z = self.displacements[..., 0], self.displacements[..., 1], self.displacements[..., 2]
        xy2 = np.sum(self.displacements[..., :2] ** 2, axis=-1)
        # Radii (High order bin dimension)
        self.spherical_coords[..., :, 0] = np.sqrt(xy2 +  z ** 2)

        # Zenith (theta for legacy reasons; middle bin dimension)
        self.spherical_coords[..., :, 1] = np.arctan2(np.sqrt(xy2), z)

        # Azimuth (phi for legacy reasons; low order bin dimension)
        self.spherical_coords[..., :, 2] = np.arctan2(y, x) + np.pi

        # Force-zero identity
        self.spherical_coords[np.diag_indices(self.num_residues)] = 0

        self.max_radii = self.spherical_coords[..., 0].max(axis=0)
        self.max_radius = self.max_radii.max()

        self.state = 5

        return self.spherical_coords

    def quantize_to_bins(self):
        if self.state < 5:
            raise InvalidContextBuilderState("Cannot bin without spherical coordinates")
        self.logger.debug("6) Quantizing to bins")
        # Bin indices are small numbers (typically between 0 and 6)
        # It's possible to have larger
        self.bins = np.empty(self.spherical_coords.shape, dtype=np.uint8)

        self.max_radial_bin = int(math.ceil(math.log(self.max_radius)) / self.radial_bin_size)
        points = np.arange(0, self.max_radial_bin * self.radial_bin_size, self.radial_bin_size)
        self.radial_bin_boundaries = np.exp(points + 1)

        # Radial binning
        #self.bins[..., 0] = self._qunatize_radial_bins(self.spherical_coords[..., 0])
        radii = self.spherical_coords[..., 0].copy()   # searchsorted is super fast if haystack array is contiguous
        self.bins[..., 0] = np.searchsorted(self.radial_bin_boundaries, radii)
        self.max_radial_bins = self.bins[..., 0].max(axis=0)
        self.num_bins = self.max_radial_bin * self.one_radial_step_bins

        # Zenith binning
        # Again searchsorted is super fast if array is contiguous
        zeniths = self.spherical_coords[..., 1].copy().flatten()
        zenith_bins = np.searchsorted(self.zenith_boundaries, zeniths)
        zenith_bins[zenith_bins > 0] -= 1
        self.bins[..., 1] = zenith_bins.reshape(self.spherical_coords[..., 1].shape)


        # Azmuth binning
        self.bins[..., 2] = self.spherical_coords[..., 2] / self.bin_angle

        self.logger.info("Highest radial bin index: {} (Histogram size: {})".format(
            self.max_radial_bin,
            self.num_bins
        ))
        self.state = 6

        return self.bins

    def assign_bin_indices(self):
        if self.state < 6:
            raise InvalidContextBuilderState("Cannot generate bin indices without binned values")

        self.logger.debug("7) Indexing bins")
        # Number of bins is 2 x bin_divisor x bin_divisor x highest_bin_index
        # By default this is 2 x 6 x 6 x [dynamic] and is usually no larger than 360
        self.indices = np.zeros((self.num_residues, self.num_residues), dtype=np.uint16)

        zenith_multiplier = 2 * self.bin_divisor
        radius_multiplier = self.one_radial_step_bins
        multipliers = np.array([
            radius_multiplier,
            zenith_multiplier,
            1
        ])
        self.indices[...] = self.bins[...].dot(multipliers)

        self.state = 7

        return self.indices

    def generate_inverted_histograms(self):
        if self.state < 7:
            raise InvalidContextBuilderState("Cannot generate inverted histograms without indices")

        self.logger.debug("8) Generating inverted histograms")

        self.inverted_histograms = self.indices.transpose()

        return self.inverted_histograms

    def generate_histograms(self):
        if self.state < 7:
            raise InvalidContextBuilderState("Cannot generate histograms without indices")

        self.logger.debug("8) Generating histograms")

        builder = functools.partial(np.bincount, minlength=self.num_bins)
        self.histograms = np.apply_along_axis(builder, -1, self.indices)ZZZZZ
        self.histograms = self.histograms.astype(np.uint16)

        self.state = 8

        return self.histograms

    def adjust_bins(self):
        if self.state < 8:
            raise InvalidContextBuilderState("Cannot adjust histograms before they have been generated")

        if self.params.merge_low_bins:
            self.logger.debug("9) Merging first and second histogram bins (rarely anything within 2.7A)")
            self.full_histograms = self.histograms
            skip = self.one_radial_step_bins
            new_histograms = self.histograms[..., skip:]
            new_histograms[..., skip:] += self.histograms[..., :skip]
            self.histograms = new_histograms

            self.full_inverted_histograms = self.inverted_histograms
            self.inverted_histograms = self.inverted_histograms - 1
            self.inverted_histograms[self.inverted_histograms < 0] = 0
        else:
            self.logger.debug("9) Skipping low/small bin merge")
            self.full_histograms = self.histograms
            self.full_inverted_histograms = self.inverted_histograms

        self.state = 10

    def generate_contexts(self):
        if self.state < 8:
            raise InvalidContextBuilderState("Cannot create contexts before histograms have been generated")

        self.logger.debug("10) Generating Residue Contexts")

        ContextClass = self.Context
        self.contexts = []
        for index in range(self.num_residues):
            atom = self.alpha_carbons[index]
            residue = atom.parent
            name = self._generate_context_name(atom)

            x, y, z = atom.get_coord()
            histogram = self.histograms[index]
            inverted = self.inverted_histograms[index]
            context = ContextClass(
                index=index,
                name=name,
                serial=atom.get_serial_number(),
                aacode=residue_three_code_to_one_code(residue.get_resname()),
                size=self.max_radii[index],
                x=x,
                y=y,
                z=z,
                histogram=histogram,
                inverted=inverted
            )
            self.contexts.append(context)

        self.state = 11

        return  self.contexts

    def run(self):
        self.load_structure()
        self.extract_alpha_carbons()
        self.build_coordinate_matrix()
        self.calculate_displacement_matrices()
        self.calculate_spherical_coordinates()
        self.quantize_to_bins()
        self.assign_bin_indices()
        self.generate_inverted_histograms()
        self.generate_histograms()
        self.adjust_bins()
        self.generate_contexts()




def main(args, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    parser = ArgumentParser()
    params = parser.parse_args(args)
    builder = ContextBuilder(params)
    builder.run()

    for context in builder.contexts:
        print(context.name, context.serial, context.aacode, context.size, end=' ')
        print('|', end=' ')
        print(context.x, context.y, context.z, end='')
        print('|', end=' ')
        for v in context.histogram:
            print(v, end=' ')
        print('|', end=' ')
        for index, v in enumerate(context.inverted):
            if index == context.index:
                print(-1, end=' ')
            else:
                print(v, end=' ')
        print()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
