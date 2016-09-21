#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six import StringIO

import numpy as np
import sh

try:
    from residuecontext import phigrid
    from residuecontext.pdbfiles import split_ext_gz
    from residuecontext.prepare_structure import get_charged_pdb
    from residuecontext.vdwgrid import load_vdw_grids
except ImportError:

    import phigrid
    from pdbfiles import split_ext_gz
    from prepare_structure import get_charged_pdb
    from vdwgrid import load_vdw_grids

CHEMGRID_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'chemgrid')

# QNIFFT has an issue with (longer) absolute paths, it seems
CHEMGRID_PARAM_FILE = "INCHEM"
CHEMGRID_INPUT_FILE = 'input.pdb'
CHEMGRID_CHARGE_FILE = 'prot.table.ambcrg.ambH'
CHEMGRID_VDW_FILE = 'vdw.parms.amb.mindock'

CHEMGRID_CHARGE_PATH = os.path.join(CHEMGRID_DIR, CHEMGRID_CHARGE_FILE)
CHEMGRID_VDW_PATH = os.path.join(CHEMGRID_DIR, CHEMGRID_VDW_FILE)
CHEMGRID_BOX_FILE = 'box'
CHEMGRID_SPACING = 0.2
CHEMGRID_ES_TYPE = 1  # Not important
CHEMGRID_ES_SCALE = 4  # Not important
CHEMGRID_BUMP_POLAR = 2.3
CHEMGRID_BUMP_NONPOLAR = 2.6
CHEMGRID_OUTPUT_PREFIX = 'output'
CHEMGRID_OUTPUT_BMP = CHEMGRID_OUTPUT_PREFIX + '.bmp'
CHEMGRID_OUTPUT_VDW = CHEMGRID_OUTPUT_PREFIX + '.vdw'
CHEMGRID_OUTPUT_OUTCHEM = 'OUTCHEM'

CHEMGRID_DEFAULT_DYNAMIC_PARAMS = {
    'grid_spacing': CHEMGRID_SPACING,
}
CHEMGRID_PARAM_PATHS = (
    CHEMGRID_CHARGE_PATH,
    CHEMGRID_VDW_PATH,
)
CHEMGRID_PARAMS = """
{pdb}
{charge}
{vdw}
{box}
{{grid_spacing}}
{es_type}
{es_scale}
{unknonw_1}
{bump_distance_polar} {bump_distance_nonpolar}
{output_prefix}
""".format(
    pdb=os.path.join('.', CHEMGRID_INPUT_FILE),
    charge=os.path.join('.', os.path.basename(CHEMGRID_CHARGE_FILE)),
    vdw=os.path.join('.', os.path.basename(CHEMGRID_VDW_FILE)),
    box=CHEMGRID_BOX_FILE,
    grid_spacing=CHEMGRID_SPACING,
    es_type=CHEMGRID_ES_TYPE,
    es_scale=CHEMGRID_ES_SCALE,
    unknonw_1=10,
    bump_distance_polar=CHEMGRID_BUMP_POLAR,
    bump_distance_nonpolar=CHEMGRID_BUMP_NONPOLAR,
    output_prefix=CHEMGRID_OUTPUT_PREFIX
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
ATOM      7  DUG BOX     1    {x0:8.3f}{y1:8.3f}{z1:8.3f}
ATOM      8  DUH BOX     1    {x1:8.3f}{y1:8.3f}{z1:8.3f}
CONECT    1    2    4    5
CONECT    2    1    3    6
CONECT    3    2    4    7
CONECT    4    1    3    8
CONECT    5    1    6    8
CONECT    6    2    5    7
CONECT    7    3    6    8
CONECT    8    4    5    7
""".strip() + '\n'


chemgrid_cmd = sh.Command(os.path.join(CHEMGRID_DIR, 'chemgrid'))


def run_chemgrid(pdb, box, vdw=None, outchem=None, bmp=None, log=None, params=CHEMGRID_DEFAULT_DYNAMIC_PARAMS):
    try:
        exec_dir = tempfile.mkdtemp()
        param_file = os.path.join(exec_dir, CHEMGRID_PARAM_FILE)
        cmd = chemgrid_cmd.bake(param_file)

        with open(param_file, 'w') as f:
            f.write(CHEMGRID_PARAMS.format(**params))

        for param_path in CHEMGRID_PARAM_PATHS:
            os.symlink(
                param_path,
                os.path.join(exec_dir, os.path.basename(param_path))
            )
        os.symlink(
            pdb,
            os.path.join(exec_dir, CHEMGRID_INPUT_FILE)
        )
        os.symlink(
            box,
            os.path.join(exec_dir, CHEMGRID_BOX_FILE)
        )
        logging.info("Running chemgrid {0!s} in {1}".format(cmd, exec_dir))
        cmd(_err_to_out=True, _cwd=exec_dir, _out=log)

        if outchem is not None:
            shutil.move(
                os.path.join(exec_dir, CHEMGRID_OUTPUT_OUTCHEM),
                outchem
            )
        if bmp is not None:
            shutil.move(
                os.path.join(exec_dir, CHEMGRID_OUTPUT_BMP),
                bmp
            )
        if vdw is not None:
            shutil.move(
                os.path.join(exec_dir, CHEMGRID_OUTPUT_VDW),
                vdw
            )
    except Exception:
        raise
    else:
        shutil.rmtree(exec_dir)


def convert_vdwgrid_to_phi(vdw_grid, box, vdw_vdw_phi):
    pass


def get_vanderderwaals_grids(code, chain=None, model=0, alignment_id=None,
                             force=False, log=None,
                             box=None, scale=None):
    protonated = get_charged_pdb(code,
                                 chain=chain,
                                 model=model,
                                 alignment_id=alignment_id,
                                 force=force)
    base, ext = split_ext_gz(protonated)
    box_file = base + '.box'
    vdw = base + '.vdw'
    bmp = base + '.bmp'
    outchem = os.path.join(os.path.dirname(base), 'OUTCHEM')

    if not os.path.exists(box_file) or box is not None:
        with open(box_file, 'w') as f:
            write_vanderwaals_box(f, box)

    if force or not os.path.exists(vdw):
        params = {}
        if scale is not None:
            params['grid_spacing'] = scale
        run_chemgrid(protonated, box_file,
                     vdw=vdw,
                     outchem=outchem,
                     bmp=bmp,
                     log=log,
                     params=params)

    vdw_grids = load_vdw_grids(outchem, vdw)
    return vdw_grids


def write_vanderwaals_box(io, box):
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
