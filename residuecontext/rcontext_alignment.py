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

try:
    from residuecontext.config import (
        RESIDUE_CONTEXT_DIR,
        JAVA_RUNTIME,
    )
    from residuecontext.pdbfiles import split_ext_gz, get_pdb_selection
except ImportError:
    from config import (
        RESIDUE_CONTEXT_DIR,
        JAVA_RUNTIME,
    )
    from pdbfiles import split_ext_gz, get_pdb_selection


RESIDUE_CONTEXT_CLASSPATH = os.path.join(RESIDUE_CONTEXT_DIR, 'jars', '*')
RESIDUE_CONTEXT_EMD_PATH = os.path.join(RESIDUE_CONTEXT_DIR, "emd", "runEMD")
RESIDUE_CONTEXT_ENV = os.environ.copy()
RESIDUE_CONTEXT_ENV['PATH'] += os.pathsep + RESIDUE_CONTEXT_EMD_PATH

ALIGNMENT_TEMP_PDB = "aligned.pdb"
ALIGNMENT_OUTPUT = "workingalignment.txt"

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
        working_dir = os.path.join(exec_dir, 'working')
        os.mkdir(working_dir)

        ident1 = code1.split()
        ident2 = code2.split()

        pdbid1 = ident1[0]
        pdbid2 = ident2[0]

        if len(ident1) > 1:
            chain1 = ident1[1]
        else:
            chain1 = 'A'

        if len(ident2) > 1:
            chain2 = ident2[1]
        else:
            chain2 = 'A'

        # Setup command based on parameters
        cmd = _residue_context_cmd.bake(
            pdbid1.upper(),
            pdbid2.upper(),
            chain1,
            chain2,
            "working",
            *RESIDUE_CONTEXT_COMMON_ARGS
        )

        # ResidueContext is picky about names
        # And also requires all PDB record types to be 6 characters
        # But biopython writes them out without padding
        for (pdb, code) in [(pdb1, ident1[0]), (pdb2, ident2[0])]:
            entfile = "{}.pdb".format(code.upper())
            with open(pdb) as src, \
                    open(os.path.join(working_dir, entfile), 'w') as dst:
                for line in src:
                    if line == 'END\n':
                        dst.write('END   \n')
                    elif line == 'TER\n':
                        dst.write('TER   \n')
                    else:
                        dst.write(line)

        logging.info("Running residue context: {0!s} in {1}".format(cmd, exec_dir))
        cmd(
            _cwd=exec_dir,
            _err=log,
            _env=RESIDUE_CONTEXT_ENV
        )

        pdb2_transformed = get_pdb_selection(
            pdbid2.upper(),
            chain=chain2,
            root=working_dir,
        )

        shutil.move(
            os.path.join(exec_dir, ALIGNMENT_OUTPUT),
            alignment
        )
        shutil.move(
            os.path.join(working_dir, "{0}.pdb-transform.txt".format(ident2[0].upper())),
            translation
        )
        shutil.move(
            pdb2_transformed,
            transformed
        )
    except Exception as e:
        raise
    else:
        shutil.rmtree(exec_dir)
