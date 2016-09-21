#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six import StringIO

import sh

try:
    from residuecontext import phigrid
    from residuecontext.config import QNIFFT_DIR
    from residuecontext.pdbfiles import split_ext_gz
    from residuecontext.prepare_structure import get_charged_pdb
except ImportError:
    import phigrid
    from config import QNIFFT_DIR
    from pdbfiles import split_ext_gz
    from prepare_structure import get_charged_pdb

def grid_parser(path):
    return phigrid.phi(path, is64=False, byteswap=False)


QNIFFT_SIZE = 193

# QNIFFT has an issue with (longer) absolute paths, it seems
QNIFFT_INPUT_FILE = 'input.pdb'
QNIFFT_OUTPUT_CRG = 'output.crg'
QNIFFT_OUTPUT_PHI = 'output.phi'
QNIFFT_CHARGE_PATH = os.path.join(QNIFFT_DIR, 'amb.crg.oxt')
QNIFFT_RADIUS_PATH = os.path.join(QNIFFT_DIR, 'default.siz')
QNIFFT_DEFAULT_PARAMS = os.path.join(QNIFFT_DIR, 'delphi.def')
QNIFFT_USER_PARAMS = 'qnifft.params'
QNIFFT_PARAM_PATHS = (
    QNIFFT_DEFAULT_PARAMS,
    QNIFFT_CHARGE_PATH,
    QNIFFT_RADIUS_PATH
)
QNIFFT_OVERRIDE_PARAMS = """
# DOCK37 Parameters
grid={size}
charge={charge}
radius={radius}
pdb_input={pdb}
pdb_output_file={crg}
phi_output_file={phi}
border=10
sizing=border
border_solvent=10.
""".format(
    size=QNIFFT_SIZE,
    charge=os.path.join('.', os.path.basename(QNIFFT_CHARGE_PATH)),
    radius=os.path.join('.', os.path.basename(QNIFFT_RADIUS_PATH)),
    pdb=os.path.join('.', QNIFFT_INPUT_FILE),
    crg=os.path.join('.', QNIFFT_OUTPUT_CRG),
    phi=os.path.join('.', QNIFFT_OUTPUT_PHI)
).strip()


qnifft_cmd = sh.Command(os.path.join(QNIFFT_DIR, 'qnifft22_{0}'.format(QNIFFT_SIZE)))


def run_qnifft(pdb, crg, phi, log=None):
    try:
        exec_dir = tempfile.mkdtemp()
        param_file = os.path.join(exec_dir, QNIFFT_USER_PARAMS)
        with open(param_file, 'w') as params:
            with open(QNIFFT_DEFAULT_PARAMS) as defs:
                params.write(defs.read())
            params.write(
                QNIFFT_OVERRIDE_PARAMS.format(
                    pdb=pdb,
                    crg=crg,
                    phi=phi
                )
            )
        cmd = qnifft_cmd.bake(QNIFFT_USER_PARAMS)

        for param_path in QNIFFT_PARAM_PATHS:
            os.symlink(
                param_path,
                os.path.join(exec_dir, os.path.basename(param_path))
            )
        os.symlink(
            pdb,
            os.path.join(exec_dir, QNIFFT_INPUT_FILE)
        )
        logging.info("Running QNIFFT {0!s} in {1}".format(cmd, exec_dir))
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
