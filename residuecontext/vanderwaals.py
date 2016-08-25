#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six import StringIO

import sh

try:
    from .pdbfiles import split_ext_gz
    from .prepare_structure import get_charged_pdb
except (ValueError, SystemError):
    from pdbfiles import split_ext_gz
    from prepare_structure import get_charged_pdb


def grid_parser(path):
    return phigrid.phi(path, is64=False, byteswap=False)


CHEMGRID_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'chemgrid')

# QNIFFT has an issue with (longer) absolute paths, it seems
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
CHEMGRID_OUTPUT_PREFIX = 'vdw'

CHEMGRID_PARAMS = os.path.join(CHEMGRID_DIR, 'INCHEM')
CHEMGRID_PARAM_PATHS = (
    CHEMGRID_PARAMS,
    CHEMGRID_CHARGE_PATH,
    CHEMGRID_VDW_PATH,
)
CHEMGRID_OVERRIDE_PARAMS = """
{pdb}
{charge}
{vdw}
{box}
{grid_spacing}
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


chemgrid_cmd = sh.Command(os.path.join(CHEMGRID_DIR, 'chemgrid'))


def run_chemgrid(pdb, crg, phi, log=None):
    try:
        exec_dir = tempfile.mkdtemp()
        cmd = chemgrid_cmd.bake()

        for param_path in CHEMGRID_PARAM_PATHS:
            os.symlink(
                param_path,
                os.path.join(exec_dir, os.path.basename(param_path))
            )
        os.symlink(
            pdb,
            os.path.join(exec_dir, CHEMGRID_INPUT_FILE)
        )
        logging.info("Running chemgrid {0!s} in {1}".format(cmd, exec_dir))
        cmd(_err_to_out=True, _cwd=exec_dir, _out=log)
        shutil.move(
            os.path.join(exec_dir, QNIFFT_OUTPUT_CRG),
            crg
        )
        shutil.move(
            os.path.join(exec_dir, QNIFFT_OUTPUT_PHI),
            phi
        )
    except Exception:
        raise
    else:
        shutil.rmtree(exec_dir)


def get_electrostatics_grid(code, chain=None, model=0, alignment_id=None, force=False, log=None):
    protonated = get_charged_pdb(code,
                                 chain=chain,
                                 model=model,
                                 alignment_id=alignment_id,
                                 force=force)
    base, ext = split_ext_gz(protonated)
    charged = base + '.crg'
    phi = base + '.phi'
    if force or not os.path.exists(phi):
        run_qnifft(protonated, charged, phi, log=log)
    phi_grid = grid_parser(phi)
    return phi_grid
