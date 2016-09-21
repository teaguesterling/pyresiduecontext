#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six import StringIO

import sh

try:
    from residuecontext.config import REDUCE_DIR
    from residuecontext.pdbfiles import (
        get_pdb_selection,
        split_ext_gz,
    )
except ImportError:
    from config import REDUCE_DIR
    from pdbfiles import (
        get_pdb_selection,
        split_ext_gz,
    )

reduce_cmd = sh.Command(os.path.join(REDUCE_DIR, 'reduce')).bake(
    '-db', os.path.join(REDUCE_DIR, 'reduce_wwPDB_het_dict.txt'),
    '-OH',
    '-HIS',
    '-FLIP',
    '-ALLALT',
    '-ROTNH3',
    '-Keep',
    '-METALBump1.5',
    '-NONMETALBump-5.0'
)


def run_reduce(pdb, crg, log=None):
    cmd = reduce_cmd.bake(pdb)
    try:
        exec_dir = tempfile.mkdtemp()
        logging.info("Running Reduce {0!s} in {1}".format(cmd, exec_dir))
        cmd(_cwd=exec_dir, _out=crg, _err=log)
    except Exception:
        raise
    else:
        shutil.rmtree(exec_dir)


def get_charged_pdb(code, chain=None, model=0, alignment_id=None, force=False, log=None):
    selection_path = get_pdb_selection(code,
                                       chain=chain,
                                       model=model,
                                       alignment_id=alignment_id)
    base, ext = split_ext_gz(selection_path)
    protonated = base + '.H'
    if force or not os.path.exists(protonated):
        run_reduce(selection_path, protonated, log=log)
    return protonated
