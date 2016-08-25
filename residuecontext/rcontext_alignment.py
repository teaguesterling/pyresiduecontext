#!/usr/bin/env python
# residuecontext.algorithm.AlgorithmWrapper
# 1AZP.pdb 1C8C.pdb A A working -v --verbhist
# -DPATH=$PATH:../src/external/emd/emd/runEMD
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six import StringIO

import sh


from residuecontext.config import (
    RESIDUE_CONTEXT_DIR,
    JAVA_RUNTIME,
)
from residuecontext.pdbfiles import split_ext_gz, get_pdb_selection


RESIDUE_CONTEXT_CLASSPATH = os.path.join(RESIDUE_CONTEXT_DIR, "jars", "*")
RESIDUE_CONTEXT_EMD_PATH = os.path.join(RESIDUE_CONTEXT_DIR, "emd", "runEMD")
RESIDUE_CONTEXT_ENV = os.environ.copy()
RESIDUE_CONTEXT_ENV['PATH'] += os.pathsep + RESIDUE_CONTEXT_EMD_PATH

ALIGNMENT_TEMP_PDB = "aligned.pdb"
ALIGNMENT_OUTPUT = "alignment.txt"

RESIDUE_CONTEXT_RUNTIME_ARGS = [
    "-Xmx500M",
    "-cp",
    RESIDUE_CONTEXT_CLASSPATH,
    "residuecontext.algorithm.AlgorithmWrapper",
]

RESIDUE_CONTEXT_COMMON_ARGS = [
    '-v',
    '--verbhist',
]

_java_cmd = sh.Command(JAVA_RUNTIME)
_residue_context_cmd = _java_cmd.bake(*RESIDUE_CONTEXT_RUNTIME_ARGS)


def run_residuecontext_alignment(
        code1,
        code2,
        alignment,
        translation,
        transformed,
        pdb1=None,
        pdb2=None,
        log=None
):
    """
    :param cmd:
    :param code1:
    :param code2:
    :param alignment:
    :param transformed:
    :param superposed:
    :param code1:
    :param code2:
    :param log:
    :return:
    """

    try:
        exec_dir = tempfile.mkdtemp()

        # Setup command based on parameters
        cmd = _biojava_base_cmd.bake(
            cmd_class,
            "-pdb1",
            code1,
            "-pdb2",
            code2,
            *BIOJAVA_COMMON_ARGS
        )

        # BioJava Structure is very picky about names
        for (pdb, code) in [(pdb1, code1), (pdb2, code2)]:
            entfile = "ent{}.pdb".format(code.split('.')[0])
            os.symlink(pdb, os.path.join(exec_dir, entfile))

        logging.info("Running biojava: {0!s} in {1}".format(cmd, exec_dir))
        cmd(
            _err_to_out=True,
            _cwd=exec_dir,
            _out=ALIGNMENT_OUTPUT,
            _err=log,
            _env=RESIDUE_CONTEXT_ENV
        )

        pdb2_transformed = get_pdb_selection(
            code2[:4],
            chain=code2[4],
            model=2,
            root=ALIGNMENT_TEMP_PDB,
        )

        shutil.move(
            os.path.join(exec_dir, ALIGNMENT_OUTPUT),
            alignment
        )
        shutil.move(
            os.path.join(exec_dir, ALIGNMENT_TEMP_PDB),
            superposed
        )
        shutil.move(
            pdb2_transformed,
            translation
        )
    except Exception:
        raise
    else:
        shutil.rmtree(exec_dir)
