#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile

import sh

from residuecontext.config import (
    BIOJAVA_DIR,
    JAVA_RUNTIME,
)
from residuecontext.pdbfiles import (
    split_ext_gz,
    get_pdb_selection,
    extract_ident,
)


BIOJAVA_CLASSPATH = os.path.join(BIOJAVA_DIR, "jars", "*")
BIOJAVA_CE_CLASS = "org.biojava.bio.structure.align.ce.CeMain"
BIOJAVA_FATCAT_CLASS = "org.biojava.bio.structure.align.fatcat.FatCat"

ALIGNMENT_TEMP_PDB = "aligned.pdb"
ALIGNMENT_OUTPUT = "alignment.txt"

BIOJAVA_RUNTIME_ARGS = [
    "-Xmx500M",
    "-cp",
    BIOJAVA_CLASSPATH,
]

BIOJAVA_COMMON_ARGS = [
    "-autoFetch",
    "false",
    "-pdbDirSplit",
    "false",
    "-printFatCat",
    "-outputPDB",
    "-outFile",
    ALIGNMENT_TEMP_PDB
]

_java_cmd = sh.Command(JAVA_RUNTIME)
_biojava_base_cmd = _java_cmd.bake(*BIOJAVA_RUNTIME_ARGS)


def run_biojava_alignment(
        cmd_class,
        code1,
        code2,
        alignment,
        transformed,
        superposed,
        pdb1=None,
        pdb2=None,
        log=None
):
    try:
        exec_dir = tempfile.mkdtemp()

        pdbid1, chain1 = extract_ident(code1)
        pdbid2, chain2 = extract_ident(code2)

        pdb1 = pdb1 or '{0}.pdb'.format(pdbid1)
        pdb2 = pdb2 or '{0}.pdb'.format(pdbid2)

        # Setup command based on parameters
        cmd = _biojava_base_cmd.bake(
            cmd_class,
            "-pdb1",
            code1,
            "-pdb2",
            code2,
            "-pdbFilePath",
            exec_dir,
            *BIOJAVA_COMMON_ARGS
        )

        # BioJava Structure is very picky about names
        for pdb, code in [(pdb1, pdbid1), (pdb2, pdbid2)]:
            entfile = "pdb{0}.ent".format(code.lower())
            os.symlink(pdb, os.path.join(exec_dir, entfile))

        logging.info("Running biojava: {0!s} in {1}".format(cmd, exec_dir))
        alignment_result = cmd(
            _cwd=exec_dir,
            _err=log
        )

        with open(alignment, 'w') as f:
            f.write(alignment_result.stdout)

        superposed_pdb = os.path.join(exec_dir, ALIGNMENT_TEMP_PDB)

        pdb2_transformed = get_pdb_selection(
            pdbid2,
            chain=chain2,
            model=1,
            root=superposed_pdb,
        )

        shutil.move(
            os.path.join(exec_dir, ALIGNMENT_TEMP_PDB),
            superposed
        )
        shutil.move(
            pdb2_transformed,
            transformed
        )
    except Exception as e:
        raise
    else:
        shutil.rmtree(exec_dir)
