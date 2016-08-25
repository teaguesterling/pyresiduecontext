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
except (ValueError, SystemError):
    from pdbfiles import split_ext_gz


JAVA_RUNTIME = "java"

BIOJAVA_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'biojava')
BIOJAVA_CLASSPATH = os.path.join(BIOJAVA_DIR, "jars", "*")
BIOJAVA_CE_CLASS = "org.biojava.bio.structure.align.ce.CeMain"
BIOJAVA_FATCAT_CLASS = "org.biojava.bio.structure.align.fatcat.FatCat"

ALIGNMENT_TEMP_PDB = "aligned.pdb"
ALIGNMENT_OUTPUT = "alignment.txt"

BIOJAVA_RUNTIME_ARGS = [
    "-Xmx500M",
    "-cp",
    BIOJAVA_CLASSPATH
]

BIOJAVA_COMMON_ARGS = [
    "-autoFetch",
    "false",
    "-printCE",
    "-outputPDB",
    "-outFile",
    ALIGNMENT_TEMP_PDB
]

_biojava_base_cmd = sh.Command(JAVA_RUNTIME, *BIOJAVA_COMMON_ARGS)
ce_cmd = _biojava_base_cmd.bake(BIOJAVA_CE_CLASS)
fatcat_cmd = _biojava_base_cmd.bake(BIOJAVA_FATCAT_CLASS)


def run_biojava_alignment(cmd, code1, code2, alignment, superposed, pdb1=None, pdb2=None, log=None):
    """

    :param cmd:
    :param code1:
    :param code2:
    :param alignment:
    :param superposed:
    :param code1:
    :param code2:
    :param log:
    :return:
    """
    
    try:
        exec_dir = tempfile.mkdtemp()
        cmd = cmd.bake(
            "-pdb1",
            *BIOJAVA_COMMON_ARGS
        )

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
