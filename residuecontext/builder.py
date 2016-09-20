#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import argparse

import logging
import math
import sys

from Bio.PDB.Entity import Entity
from collections import namedtuple

import numpy as np
from numpy import ma as ma
from Bio.PDB.Polypeptide import three_to_one as residue_three_code_to_one_code

import time


from residuecontext.pdbfiles import PDBParser


def ball_mask(shape):
    shape = np.array(shape)
    even = np.any(shape % 2 == 0)
    centers = shape // 2
    radius = centers.min()
    indices = np.array([np.arange(0, d) for d in shape])
    grid = np.array(np.meshgrid(*indices, indexing='ij'))
    if even:
        #radius -= np.sqrt(2) / 2
        grid += 1
    repeater = tuple([slice(None, None)] + [np.newaxis] * len(shape))
    centered = grid - centers.T[repeater]
    dists = np.linalg.norm(centered.T, axis=-1)
    masked = dists > radius
    return masked


def array_from_grid(grid):
    if hasattr(grid, 'phiArray'):
        return array_from_phi(grid)
    elif hasattr(grid, 'vdwArray'):
        return array_from_vdw(grid)
    elif hasattr(grid, 'saArray'):
        return array_from_sa(grid)


def array_from_phi(phi):
    array = np.array(phi.phiArray).reshape(phi.gridDimension, phi.gridDimension, phi.gridDimension)
    return array


def array_from_vdw(vdw):
    array = np.array(vdw.vdwArray).reshape(vdw.gridDimension, vdw.gridDimension, vdw.gridDimension)
    return array


def array_from_sa(sa):
    array = np.array(sa.saArray).reshape(sa.gridDimension, sa.gridDimension, sa.gridDimension)
    return array



def get_grid_sphere_at_point(phi, center, boxStep=None, boxSize=np.e ** 2):
    step = 1 / phi.scale
    dim = phi.gridDimension
    grid_center = np.array(phi.oldmid)
    min_corner = grid_center - dim // 2 * step
    max_corner = grid_center + dim // 2 * step
    extents = np.array([
        min_corner,
        max_corner
    ])
    boxStep = boxStep or step
    #GP1 = [np.linspace(X0, X1, gridDim) for X0, X1 in gridExtents.transpose()]

    array = array_from_grid(phi)
    grid_coords = [np.arange(C - boxSize, C + boxSize, boxStep) for C in center]

    XYZ = np.transpose(grid_coords)
    X, Y, Z = XYZ.transpose()
    cXYZ = np.zeros((3, len(X), len(Y), len(Z)))
    cBox = np.zeros((len(X), len(Y), len(Z)))

    # Slow for debugging
    t = time.time()
    for i, x in enumerate(X):
        for j, y in enumerate(Y):
            for k, z in enumerate(Z):
                idx = phi.getIndices((x, y, z))
                #rXYZ[i, j, k] = np.floor((np.array([x, y, z]) - G2m) / Gstep2).astype(np.int)
                #d2XYZ[:, i, j, k] = (i * boxStep, j * boxStep, k * boxStep)
                cXYZ[:, i, j, k] = x, y, z
                cBox[i,j,k] = array[idx]

    t = time.time()
    # XXX: Removing sphere masking
    #pSphere = np.sign(cBox) * np.clip(np.log(np.abs(cBox)), 0, None)
    #pSphere = ma.array(pSphere, mask=ball_mask(pSphere.shape))
    pSphere = np.clip(cBox, -7, 7)

    return cXYZ, pSphere


def create_grid_histogram(phi, center, bins=50, extents=None, params=None, **kwargs):
    coords, values = get_grid_sphere_at_point(phi, center,
                                              boxStep=kwargs.pop('boxStep', None),
                                              boxSize=kwargs.pop('boxSize', np.e ** 2))

    params = params or SphericalHistogramParams(**kwargs)
    grid_points = coords.transpose()
    points = grid_points.reshape(-1, 3)

    builder = SphericalHistogramBuilder(params)
    builder.set_coordinates([center], points)

    histogram, inverse = builder.run()

    if extents is None:
        extents = values.min(), values.max()

    num_spherical_bins = histogram.shape[-1]
    inner_bin_boundaries = np.linspace(extents[0], extents[1], bins+1)
    bin_histograms = np.empty((num_spherical_bins, bins))

    clipped = np.clip(values, extents[0], extents[1])
    flattened = clipped.ravel()
    mappings = inverse.transpose()[0]

    for bin_index in range(num_spherical_bins):
        in_bin = mappings == bin_index
        selection = flattened[in_bin]
        #XXX: Removing sphere masking
        #in_sphere = selection[~selection.mask]
        in_sphere = selection
        num_in_sphere = len(in_sphere)
        if num_in_sphere == 0:
            bin_histograms[bin_index] = 0
        else:
            bin_histogram, _boundaries = np.histogram(in_sphere, bins=inner_bin_boundaries)
            bin_histograms[bin_index] = bin_histogram / num_in_sphere

    return bin_histograms


def create_paired_grid_histogram(vdwA, vdwB, center, bins=50, extents=None, params=None, **kwargs):
    coordsA, valuesA = get_grid_sphere_at_point(vdwA, center,
                                              boxStep=kwargs.pop('boxStep', None),
                                              boxSize=kwargs.pop('boxSize', np.e ** 2))

    coordsB, valuesB = get_grid_sphere_at_point(vdwB, center,
                                                boxStep=kwargs.pop('boxStep', None),
                                                boxSize=kwargs.pop('boxSize', np.e ** 2))

    params = params or SphericalHistogramParams(**kwargs)

    grid_pointsA = coordsA.transpose()
    pointsA = grid_pointsA.reshape(-1, 3)

    grid_pointsB = coordsB.transpose()
    pointsB = grid_pointsB.reshape(-1, 3)

    builder = SphericalHistogramBuilder(params)

    builder.set_coordinates([center], pointsA)
    histogramA, inverseA = builder.run()

    builder.set_coordinates([center], pointsB)
    histogramB, inverseB = builder.run()

    if extents is None:
        extentsA = valuesA.min(), valuesA.max()
        extentsB = valuesB.min(), valuesB.max()
    elif isinstance(extents, list):
        extentsA, extentsB = extents
    else:
        extentsA = extents
        extentsB = extents

    num_spherical_bins = histogramA.shape[-1]
    half_bins = bins//2  # B is negative A is positive

    inner_bin_boundaries = np.linspace(extentsA[0], extentsA[1], half_bins+1)
    bin_histograms = np.empty((num_spherical_bins, bins))

    clippedA = np.clip(valuesA, extentsA[0], extentsA[1])
    flattenedA = clippedA.ravel()
    mappingsA = inverseA.transpose()[0]

    clippedB = np.clip(valuesB, extentsB[0], extentsB[1])
    flattenedB = clippedB.ravel()
    mappingsB = inverseB.transpose()[0]

    for bin_index in range(num_spherical_bins):
        in_binA = mappingsA == bin_index
        in_binB = mappingsB == bin_index
        selectionA = flattenedA[in_binA]
        selectionB = flattenedB[in_binB]

        in_sphereA = selectionA
        num_in_sphereA = len(in_sphereA)

        in_sphereB = selectionB
        num_in_sphereB = len(in_sphereB)

        if num_in_sphereA == 0 and num_in_sphereB == 0:
            bin_histograms[bin_index] = 0
        else:
            bin_histogramA, _boundaries = np.histogram(in_sphereA, bins=inner_bin_boundaries)
            bin_histogramB, _boundaries = np.histogram(in_sphereB, bins=inner_bin_boundaries)
            bin_histogramB_rev = bin_histogramB[::-1]  # All but "0" bin but backwards

            bin_histogram = np.zeros(bins)
            bin_histogram[half_bins:] = bin_histogramA / num_in_sphereA
            bin_histogram[:half_bins] = bin_histogramB_rev / num_in_sphereB

            bin_histograms[bin_index] = bin_histogram

    return bin_histograms


class InvalidContextBuilderState(Exception):
    pass


class SphericalHistogramBuilder:
    def __init__(self, params, *points, **kwargs):
        self.params = params

        self.radial_bin_size = self.params.radial_bin_size
        self.bin_divisor = self.params.bin_divisor

        self.logger = kwargs.pop('logger', logging)
        self.state = 0

        self._sphere_templates = None

        if points:
            self.set_coordinates(*points)

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

    def set_coordinates(self, points, environmens=None):

        points = np.array(points)
        if environmens is None:
            environmens = points
        else:
            environmens = np.array(environmens)

        self.point_coords = points
        self.environment_coords = environmens

        self.num_points = len(self.point_coords)
        self.num_environments = len(self.environment_coords)

        self.extents = np.array([
            self.point_coords.min(axis=0),
            self.point_coords.max(axis=0)
        ])

        self.state = 3

        return self.point_coords, self.environment_coords

    def calculate_displacement_matrices(self):
        if self.state < 3:
            raise InvalidContextBuilderState("Cannot calculate displacement matrices without coordinates")
        self.logger.debug("4) Calculating NxN displacements")

        # Tile N instances of coordinate matrix and subtract off a different row from each
        # Done all at once for speed
        per_point_environments = np.tile(self.environment_coords, (self.num_points, 1, 1))
        self.displacements = per_point_environments - self.point_coords[..., np.newaxis, :]

        self.state = 4
        return self.displacements

    def calculate_spherical_coordinates(self):
        if self.state < 4:
            raise InvalidContextBuilderState("Cannot calculate spherical coordinates without displacements")
        self.logger.debug("5) Converting to spherical coordinates")

        coords_shape = (self.num_points,) + self.environment_coords.shape
        self.spherical_coords = np.empty(coords_shape, dtype=np.float)
        x, y, z = self.displacements[..., 0], self.displacements[..., 1], self.displacements[..., 2]
        xy2 = np.sum(self.displacements[..., :2] ** 2, axis=-1)
        # Radii (High order bin dimension)
        self.spherical_coords[..., :, 0] = np.sqrt(xy2 +  z ** 2)

        # Zenith (theta for legacy reasons; middle bin dimension)
        self.spherical_coords[..., :, 1] = np.arctan2(np.sqrt(xy2), z)

        # Azimuth (phi for legacy reasons; low order bin dimension)
        self.spherical_coords[..., :, 2] = np.arctan2(y, x) + np.pi

        # Force-zero identity
        if self.point_coords is self.environment_coords:
            self.spherical_coords[np.diag_indices(self.num_points)] = 0

        self.max_radii = self.spherical_coords[..., 0].max(axis=0)
        self.max_radius = self.max_radii.max()

        self.state = 5

        return self.spherical_coords

    def prepare_bins(self):
        if self.state < 5:
            raise InvalidContextBuilderState("Cannot create bins without spherical coordinates")
        self.logger.debug("6) Calculating bin parameters")

        self.bin_angle = np.pi / self.bin_divisor

        # Number of bins per radial step (e.g. bins per sphere)
        self.one_radial_step_bins = 2 * self.bin_divisor * self.bin_divisor

        # Step between 1 and -1 by 2 * bin_divisor increments
        # Then compute the arc cosine at each point as a bin boundary
        self.inv_zenith_step = -2 / self.bin_divisor
        self.zenith_boundaries = np.arccos(np.arange(1, -1, self.inv_zenith_step))[1:]

        self.azimuth_boundaries = np.roll(np.arange(0, 2 * np.pi, self.bin_angle), self.bin_divisor)

        self.max_radial_bin = int(math.ceil(math.log(self.max_radius)) / self.radial_bin_size)
        self.radial_bin_boundaries = np.exp(np.arange(
            self.radial_bin_size,
            self.radial_bin_size * (self.max_radial_bin + 1),
            self.radial_bin_size
        ))

        zenith_points = list(self.zenith_boundaries)  # Leave off north pole and both in manually
        azimuth_points = list(self.azimuth_boundaries)
        num_zeniths = len(zenith_points)
        num_azimuths = len(azimuth_points)
        vertices_per_sphere = num_zeniths * num_azimuths + 2  # Poles
        vertices = np.zeros((0, 3))
        faces = np.zeros((0, 4), dtype=np.int)
        connectivity_offset = 0

        vertex_indices = np.arange(vertices_per_sphere - 2)   # Skip Poles
        sphere_faces = np.transpose([
            vertex_indices,
            num_azimuths * (vertex_indices // num_azimuths) + ((vertex_indices + 1) % num_azimuths),
            np.clip(num_azimuths + vertex_indices, 1, vertices_per_sphere - 2),
            np.clip(num_azimuths + (
                num_azimuths * (vertex_indices // num_azimuths) + ((vertex_indices + 1) % num_azimuths)
            ), 1, vertices_per_sphere - 2)
        ])
        sphere_faces += 1  # Add north pole
        sphere_faces = np.vstack((
            np.concatenate((sphere_faces[:num_azimuths, :2], np.zeros((12, 1)), np.zeros((12, 1))), axis=1),
            sphere_faces,
        ))
        spherical_vertices = np.transpose(
            np.transpose(np.meshgrid(zenith_points, azimuth_points)).reshape(-1, 2)
        )
        for radius in self.radial_bin_boundaries:
            no_pole_vertices = np.transpose([
                radius * np.cos(spherical_vertices[1]) * np.sin(spherical_vertices[0]),
                radius * np.sin(spherical_vertices[1]) * np.sin(spherical_vertices[0]),
                radius * np.cos(spherical_vertices[0])
            ])
            cartesian_vertices = np.vstack((
                (0, 0, radius),
                no_pole_vertices,
                (0, 0, -radius),
            ))
            vertices = np.vstack((vertices, cartesian_vertices))
            faces = np.vstack((faces, sphere_faces + connectivity_offset))
            connectivity_offset += vertices_per_sphere

        self.context_template = {
            'vertices': vertices,
            'faces': faces,
            'vertices_per_sphere': vertices_per_sphere,
            'faces_per_sphere': len(sphere_faces),
        }
        self.state = 6
        return self.context_template

    def quantize_to_bins(self):
        if self.state < 6:
            raise InvalidContextBuilderState("Cannot quantize without bins")
        self.logger.debug("7) Quantizing to bins")

        # Bin indices are small numbers (typically between 0 and 6)
        # It's possible to have larger
        self.bins = np.empty(self.spherical_coords.shape, dtype=np.uint8)

        # Radial binning
        #self.bins[..., 0] = self._qunatize_radial_bins(self.spherical_coords[..., 0])
        radii = self.spherical_coords[..., 0].copy()   # searchsorted is super fast if haystack array is contiguous
        # Zenith binning
        # Again searchsorted is super fast if array is contiguous
        zeniths = self.spherical_coords[..., 1].copy()
        #zenith_bins = np.searchsorted(self.zenith_boundaries, zeniths)
        #zenith_bins[zenith_bins > 0] -= 1
        #self.bins[..., ] = (
        #    np.searchsorted(self.radial_bin_boundaries, radii, side='right'),
        #    np.clip(np.searchsorted(self.zenith_boundaries, zeniths, side='right') - 1, 0, len(self.zenith_boundaries)),
        #    self.spherical_coords[..., 2] / self.bin_angle
        #)
        # self.bins[..., 1] = zenith_bins.reshape(self.spherical_coords[..., 1].shape)
        self.bins[..., 0] = np.searchsorted(self.radial_bin_boundaries, radii, side='right')
        self.bins[..., 1] = np.searchsorted(self.zenith_boundaries, zeniths)
        self.bins[..., 2] = self.spherical_coords[..., 2] / self.bin_angle

        self.max_radial_bins = self.bins[..., 0].max(axis=0)
        self.num_bins = self.max_radial_bin * self.one_radial_step_bins
        self.logger.info("Highest radial bin index: {} (Histogram size: {})".format(
            self.max_radial_bin,
            self.num_bins
        ))
        self.state = 7

        return self.bins

    def assign_bin_indices(self):
        if self.state < 7:
            raise InvalidContextBuilderState("Cannot generate bin indices without binned valus")

        self.logger.debug("8) Indexing bins")
        # Number of bins is 2 x bin_divisor x bin_divisor x highest_bin_index
        # By default this is 2 x 6 x 6 x [dynamic] and is usually no larger than 360
        self.indices = np.zeros((self.num_points, self.num_environments), dtype=np.uint16)

        zenith_multiplier = 2 * self.bin_divisor
        radius_multiplier = self.one_radial_step_bins
        multipliers = np.array([
            radius_multiplier,
            zenith_multiplier,
            1
        ])
        self.indices[...] = self.bins[...].dot(multipliers)
        if self.point_coords is self.environment_coords:
            self.indices = ma.masked_array(self.indices, mask=np.eye(len(self.indices), dtype=np.bool))

        self.state = 8

        return self.indices

    def generate_inverted_histograms(self):
        if self.state < 8:
            raise InvalidContextBuilderState("Cannot generate inverted histograms without indices")

        self.logger.debug("9) Generating inverted histograms")

        self.inverted_histograms = self.indices.transpose()

        self.state = 9

        return self.inverted_histograms

    def generate_histograms(self):
        if self.state < 8:
            raise InvalidContextBuilderState("Cannot generate histograms without indices")

        self.logger.debug("10) Generating histograms")

        if self.point_coords is self.environment_coords:
            builder = lambda arr: np.bincount(arr[~arr.mask], minlength=self.num_bins)
            self.histograms = ma.apply_along_axis(builder, -1, self.indices)
            self.histograms = self.histograms.data.astype(np.uint16)
        else:
            builder = lambda arr: np.bincount(arr, minlength=self.num_bins)
            self.histograms = np.apply_along_axis(builder, -1, self.indices)
            self.histograms = self.histograms.astype(np.uint16)

        self.state = 10

        return self.histograms


    def adjust_bins(self):
        if self.state < 10:
            raise InvalidContextBuilderState("Cannot adjust histograms before they have been generated")

        self.full_histograms = self.histograms
        self.full_inverted_histograms = self.inverted_histograms

        if self.params.merge_low_bins:
            self.logger.debug("11) Merging first and second histogram bins (rarely anything within 2.7A)")
            self.full_histograms = self.histograms
            skip = self.one_radial_step_bins
            new_histograms = self.histograms[..., skip:]
            new_histograms[..., :skip] += self.histograms[..., :skip]
            self.histograms = new_histograms

            self.full_inverted_histograms = self.inverted_histograms
            self.inverted_histograms = self.inverted_histograms - 1
            self.inverted_histograms[self.inverted_histograms < 0] = 0
        else:
            self.logger.debug("11) Skipping low/small bin merge")

        self.state = 11

    def run(self):
        self.calculate_displacement_matrices()
        self.calculate_spherical_coordinates()
        self.prepare_bins()
        self.quantize_to_bins()
        self.assign_bin_indices()
        self.generate_inverted_histograms()
        self.generate_histograms()
        self.adjust_bins()

        return self.histograms, self.inverted_histograms


class GridContextBuilder:
    state = -1

    Context = namedtuple('Context', (
        'index',
        'name',
        'serial',
        'aacode',
        'x',
        'y',
        'z',
        'histograms',
    ))


class ContextBuilder:
    state = -1

    Context = namedtuple('Context', (
        'index',
        'name',
        'serial',
        'resnum',
        'aacode',
        'size',
        'x',
        'y',
        'z',
        'histogram',
        'inverted',
    ))

    class ContextCollection:
        def __init__(self, name, bin_data, contexts, surface=None):
            self.name = name
            self.bin_data = bin_data
            self.contexts = contexts
            self.surface = surface

        @property
        def size(self):
            return len(self.contexts)

    def __init__(self, params, logger=logging):
        self.params = params

        self.logger = logger
        self.builder = SphericalHistogramBuilder(params, logger=logger)
        self._parser = PDBParser()
        self.state = 0

    def _generate_context_name(self, atom):
        residue = atom.parent
        res_id = residue.id[1]
        res_name = residue.get_resname()
        name = 'd{}{}.{}{:d}'.format(self.structure_name, self.params.chain, res_name, res_id)
        return name

    def load_structure(self):
        self.logger.debug("1) Loading structure")
        if isinstance(self.params.source, Entity):
            self.structure = self.params.source
        else:
            self.structure = self._parser.load_structure(self.params.source, name=self.params.name)
        if self.structure.get_level() == 'C':
            self.chain = self.structure
        else:
            model = self.structure.get_list()[0]  # Only consider the first model
            if self.params.chain is None:
                self.chain = list(model)[0]
            else:
                self.chain = [chain for chain in model if chain.get_id() == self.params.chain][0]

        if self.structure.id != PDBParser.NAME_UNKONWN:
            self.structure_name = self.structure.id
        else:
            self.structure_name = "____"
        self.logger.info("Loaded {} from {} (chain: )".format(
            self.structure.id,
            self.structure_name,
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

    def execute_builder(self):
        if self.state < 3:
            raise InvalidContextBuilderState("Cannot calculate displacement matrices without coordinates")

        self.builder.set_coordinates(self.coords)
        histogram, inverse_histogram = self.builder.run()

        self.histogram = histogram
        self.inverse_histogram = inverse_histogram

        self.state = 4

        return self.histogram

    def generate_residue_contexts(self):
        if self.state < 4:
            raise InvalidContextBuilderState("Cannot create contexts before histograms have been generated")

        self.logger.debug("5) Generating Residue Contexts")

        ContextClass = self.Context
        self.contexts = []
        for index in range(self.num_residues):
            atom = self.alpha_carbons[index]
            residue = atom.parent
            name = self._generate_context_name(atom)

            x, y, z = atom.get_coord()
            histogram = self.builder.histograms[index]
            inverted = self.builder.inverted_histograms[index]
            context = ContextClass(
                index=index,
                name=name,
                serial=atom.get_serial_number(),
                resnum=residue.get_id()[1],
                aacode=residue_three_code_to_one_code(residue.get_resname()),
                size=self.builder.max_radii[index],
                x=x,
                y=y,
                z=z,
                histogram=histogram,
                inverted=inverted,
            )
            self.contexts.append(context)

        self.state = 5

        return  self.contexts

    def generate_chain_context(self):
        if self.state < 5:
            raise InvalidContextBuilderState("Cannot create chain context without residue contexts")

        self.logger.debug("5) Generating Residue Contexts")

        self.chain_context = self.ContextCollection(
            name='d{}{}'.format(self.structure_name, self.params.chain),
            contexts=self.contexts,
            bin_data=dict(
                num_bins=self.builder.num_bins,
                num_bins_per_radial_unit=self.builder.one_radial_step_bins,
                bin_divisor=self.builder.bin_divisor,
                bin_angle=self.builder.bin_angle,
                radial_bin_size=self.builder.radial_bin_size,
                radial_bin_boundaries=self.builder.radial_bin_boundaries.tolist(),
                zenith_boundaries=self.builder.zenith_boundaries.tolist(),
                azimuth_boundaries=(np.arange(0, 2 * np.pi, np.pi / self.builder.bin_divisor)).tolist(),
                max_radial_bin=self.builder.max_radial_bin,
                max_radius=self.builder.max_radius,
                extents=self.extents.tolist(),
            ),
            surface=self.builder.context_template
        )

        self.state = 6

        return self.chain_context

    def run(self):
        self.load_structure()
        self.extract_alpha_carbons()
        self.build_coordinate_matrix()
        self.execute_builder()
        self.generate_residue_contexts()
        self.generate_chain_context()

        return self.chain_context


class SphericalHistogramParams(argparse.ArgumentParser):
    bin_divisor = 6
    radial_bin_size = 1.0

    merge_low_bins = False
    exclude_outer_bins = True

    def __init__(self, *args, **kwargs):
        for key in type(self).__dict__:
            if key in kwargs:
                self.__dict__[key] = kwargs.pop(key)
        super(SphericalHistogramParams, self).__init__(*args, **kwargs)

        self.add_argument('-b', '--bindivisor', dest='bin_divisor', type=int, default=self.bin_divisor)
        self.add_argument('-a', '--rbinsize', dest='radial_bin_size', type=float, default=self.radial_bin_size)
        self.add_argument('-E', '--no-exclude-outer-bins', dest='exclude_outer_bins',
                          action='store_const', const=not self.exclude_outer_bins, default=self.exclude_outer_bins)
        self.add_argument('-M', '--no-merge-low-bins', dest='merge_low_bins',
                          action='store_const', const=not self.merge_low_bins, default=self.merge_low_bins)

    def parse_args(self, args=None, namespace=None):
        super(SphericalHistogramParams, self).parse_args(args[1:] if args is sys.argv else args, namespace=self)
        return self


class ResidueContextParams(SphericalHistogramParams):
    source = None
    chain = None
    output = sys.stdout
    inverse_output = None

    name = None

    def __init__(self, *args, **kwargs):
        stdout = kwargs.pop('stdout', sys.stdout)
        super(ResidueContextParams, self).__init__(*args, **kwargs)
        self.add_argument('source', type=argparse.FileType('r'))
        self.add_argument('chain', type=str)
        self.add_argument('-o', '--output', type=argparse.FileType('w'), default=stdout)
        self.add_argument('-O', '--inverse-output', dest='inverse_output', type=argparse.FileType('w'), default=None)
        self.add_argument('-N', '--name', dest='name', default=None)


def main(args, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    parser = ResidueContextParams()
    params = parser.parse_args(args)
    builder = ContextBuilder(params)
    builder.run()

    for context in builder.contexts:
        print(context.name, context.serial, context.aacode, context.size, end=' ')
        print('|', end=' ')
        print(context.x, context.y, context.z, end=' ')
        print('|', end=' ')
        for v in context.histogram:
            print(v, end=' ')
        print('|', end=' ')
        for index, v in enumerate(context.inverted):
            print(v, end=' ')
        print()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
