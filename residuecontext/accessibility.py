#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six import StringIO

import numpy as np
import sh

from residuecontext import phigrid
from residuecontext.pdbfiles import split_ext_gz
from residuecontext.prepare_structure import get_charged_pdb
from residuecontext.accessibilitygrid import load_solv_grid

SOLVMAP_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'solvmap')

# QNIFFT has an issue with (longer) absolute paths, it seems
SOLVMAP_PARAM_FILE = "INSEV"
SOLVMAP_INPUT_FILE = 'input.pdb'

SOLVMAP_BOX_FILE = 'box'
SOLVMAP_SPACING = 1
SOLVMAP_ATOM_RADII = "1.60,1.65,1.90,1.90,1.90,1.60"
SOLVMAP_PROBE_RADIUS = 1.4
SOLVMAP_ASSUMED_LIGAND_RADIUS = 2.6

SOLVMAP_OUTSEV = 'OUTSEV'
SOLVMAP_OUTPUT = 'output.solv-heavy'

SOLVMAP_PARAMS = """
{pdb}
{output}
{radii}
{probe}
{spacing}
{box}
{ligand}
""".format(
    pdb=os.path.join('.', SOLVMAP_INPUT_FILE),
    output=os.path.join('.', SOLVMAP_OUTPUT),
    radii=SOLVMAP_ATOM_RADII,
    probe=SOLVMAP_PROBE_RADIUS,
    spacing=SOLVMAP_SPACING,
    ligand=SOLVMAP_ASSUMED_LIGAND_RADIUS,
    box=SOLVMAP_BOX_FILE
).strip()


BOX_TEMPLATE = """
HEADER    CORNERS OF BOX {x0:8.3f}{y0:8.3f}{z0:8.3f}{x1:8.3f}{y1:8.3f}{z1:8.3f}
REMARK    CENTER (X Y Z) {cx:8.3f}{cy:8.3f}{cz:8.3f}
REMARK    DIMENSIONS (X Y Z) {dx:8.3f}{dy:8.3f}{dz:8.3f}
ATOM      1  DUA BOX     1    {x0:8.3f}{y0:8.3f}{z0:8.3f}
ATOM      2  DUB BOX     1    {x1:8.3f}{y0:8.3f}{z0:8.3f}
ATOM      3  DUC BOX     1    {x1:8.3f}{y0:8.3f}{z1:8.3f}
ATOM      4  DUD BOX     1    {x0:8.3f}{y0:8.3f}{z1:8.3f}
ATOM      5  DUE BOX     1    {x0:8.3f}{y1:8.3f}{z0:8.3f}
ATOM      6  DUF BOX     1    {x1:8.3f}{y1:8.3f}{z0:8.3f}
ATOM      7  DUG BOX     1    {x1:8.3f}{y1:8.3f}{z1:8.3f}
ATOM      8  DUH BOX     1    {x0:8.3f}{y1:8.3f}{z1:8.3f}
CONECT    1    2    4    5
CONECT    2    1    3    6
CONECT    3    2    4    7
CONECT    4    1    3    8
CONECT    5    1    6    8
CONECT    6    2    5    7
CONECT    7    3    6    8
CONECT    8    4    5    7
""".strip() + '\n'


solvmap_cmd = sh.Command(os.path.join(SOLVMAP_DIR, 'solvmap'))


def run_solvmap(pdb, box, outsev=None, solv=None, log=None):
    try:
        exec_dir = tempfile.mkdtemp()
        param_file = os.path.join(exec_dir, SOLVMAP_PARAM_FILE)
        cmd = solvmap_cmd.bake(param_file)

        with open(param_file, 'w') as f:
            f.write(SOLVMAP_PARAMS)

        os.symlink(
            pdb,
            os.path.join(exec_dir, SOLVMAP_INPUT_FILE)
        )
        os.symlink(
            box,
            os.path.join(exec_dir, SOLVMAP_BOX_FILE)
        )
        logging.info("Running sovlmap {0!s} in {1}".format(cmd, exec_dir))
        cmd(_err_to_out=True, _cwd=exec_dir, _out=log)

        if outsev is not None:
            shutil.move(
                os.path.join(exec_dir, SOLVMAP_OUTSEV),
                outsev
            )
        if solv is not None:
            shutil.move(
                os.path.join(exec_dir, SOLVMAP_OUTPUT),
                solv
            )
    except Exception:
        raise
    else:
        shutil.rmtree(exec_dir)


def get_solvation_grid(code, chain=None, model=0, alignment_id=None,
                       force=False, log=None,
                       box=None):
    protonated = get_charged_pdb(code,
                                 chain=chain,
                                 model=model,
                                 alignment_id=alignment_id,
                                 force=force)
    base, ext = split_ext_gz(protonated)
    box_file = base + '.box'
    solv = base + '.solv-heavy'
    outsev = os.path.join(os.path.dirname(base), SOLVMAP_OUTSEV)

    if not os.path.exists(box_file) or box is not None:
        with open(box_file, 'w') as f:
            write_solventaccessibility_box(f, box)

    if force or not os.path.exists(solv):
        run_solvmap(protonated, box_file,
                    outsev=outsev,
                    solv=solv,
                    log=log)

    solv_grid = load_solv_grid(solv)
    return solv_grid


def write_solventaccessibility_box(io, box):
    if isinstance(box, phigrid.phi):
        extents = box.getMinsMaxs()
        center = box.oldmid
        box_def = {
            'x0': extents[0][0],
            'x1': extents[1][0],
            'y0': extents[0][1],
            'y1': extents[1][1],
            'z0': extents[0][2],
            'z1': extents[1][2],
            'dx': extents[1][0] - extents[0][0],
            'dy': extents[1][1] - extents[0][1],
            'dz': extents[1][2] - extents[0][2],
            'cx': center[0],
            'cy': center[1],
            'cz': center[2],
        }

    elif len(box) == 3:
        box_def = {
            'x0': box[0][0],
            'x1': box[0][1],
            'y0': box[1][0],
            'y1': box[1][1],
            'z0': box[2][0],
            'z1': box[2][1],
            'dx': box[0][1] - box[0][0],
            'dy': box[1][1] - box[1][0],
            'dz': box[2][1] - box[2][0],
            'cx': (box[0][1] - box[0][0]) / 2,
            'cy': (box[1][1] - box[1][0]) / 2,
            'cz': (box[2][1] - box[2][0]) / 2
        }
    elif len(box) == 2:
        box_def = {
            'x0': box[0][0],
            'x1': box[1][0],
            'y0': box[0][1],
            'y1': box[1][1],
            'z0': box[0][2],
            'z1': box[1][2],
            'dx': box[1][0] - box[0][0],
            'dy': box[1][1] - box[0][1],
            'dz': box[1][2] - box[0][2],
            'cx': (box[1][0] - box[0][0]) / 2,
            'cy': (box[1][1] - box[0][1]) / 2,
            'cz': (box[1][2] - box[0][2]) / 2
        }
    elif len(box) == 6:
        box_def = {
            'x0': box[0],
            'x1': box[3],
            'y0': box[1],
            'y1': box[4],
            'z0': box[2],
            'z1': box[5],
            'dx': box[3] - box[0],
            'dy': box[4] - box[1],
            'dz': box[5] - box[2],
            'cx': (box[3] - box[0]) / 2,
            'cy': (box[4] - box[1]) / 2,
            'cz': (box[5] - box[2]) / 2
        }

    box_data = BOX_TEMPLATE.format(**box_def)
    io.write(box_data)
