#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile

from celery import Celery

from residuecontext.config import ALIGNMENT_JOB_DIR
from residuecontext.pdbfiles import (
    get_pdb_selection,
    extract_ident,
)
from residuecontext.biojava_alignment import (
    BIOJAVA_CE_CLASS,
    BIOJAVA_FATCAT_CLASS,
    run_biojava_alignment,
)
from residuecontext.rcontext_alignment import run_residuecontext_alignment
from residuecontext.electrostatics import get_electrostatics_grid
from residuecontext.vanderwaals import get_vanderderwaals_grid


queue = Celery('tasks',
               backend='db+sqlite:///celerydb.sqlite',
               broker='sqla+sqlite:///celerydb.sqlite')
queue.conf['CELERY_RESULT_DB_SHORT_LIVED_SESSIONS'] = True

@queue.task
def run_alignment_comparisons(number, code1, code2):
    working_dir = os.path.join(ALIGNMENT_JOB_DIR, number)

    pdbid1, chain1 = extract_ident(code1)
    pdbid2, chain2 = extract_ident(code2)

    pdb1src = get_pdb_selection(pdbid1, chain=chain1)
    pdb2src = get_pdb_selection(pdbid2, chain=chain2)

    pdb1 = os.path.join(working_dir, os.path.basename(pdb1src))
    pdb2 = os.path.join(working_dir, os.path.basename(pdb2src))

    shutil.copy(pdb1src, pdb1)
    shutil.copy(pdb2src, pdb2)

    rcontext_dir = os.path.join(working_dir, "rcontext")
    ce_dir = os.path.join(working_dir, "ce")
    fatcat_dir = os.path.join(working_dir, "fatcat")
    dali_dir = os.path.join(working_dir, "dali")

    alignment_basename = 'alignment.txt'
    translated_basename = 'transformed-{0}.pdb'.format(code2)
    superposed_basename = 'superposed.pdb'
    transform_basename = '{0}-transform.txt'

    os.mkdir(ce_dir)
    run_biojava_alignment(
        BIOJAVA_CE_CLASS,
        "{0}.{1}".format(pdbid1, chain1),
        "{0}.{1}".format(pdbid2, chain2),
        alignment=os.path.join(ce_dir, alignment_basename),
        transformed=os.path.join(ce_dir, translated_basename),
        superposed=os.path.join(ce_dir, superposed_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )

    os.mkdir(fatcat_dir)
    run_biojava_alignment(
        BIOJAVA_FATCAT_CLASS,
        "{0}.{1}".format(pdbid1, chain1),
        "{0}.{1}".format(pdbid2, chain2),
        alignment=os.path.join(fatcat_dir, alignment_basename),
        transformed=os.path.join(fatcat_dir, translated_basename),
        superposed=os.path.join(fatcat_dir, superposed_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )

    os.mkdir(rcontext_dir)
    run_residuecontext_alignment(
        "{0}.{1}".format(pdbid1, chain1),
        "{0}.{1}".format(pdbid2, chain2),
        alignment=os.path.join(rcontext_dir, alignment_basename),
        transformed=os.path.join(rcontext_dir, translated_basename),
        translation=os.path.join(rcontext_dir, transform_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )

    os.mkdir(dali_dir)
    #run_dali_alignment(
    #    "{0}.{1}".format(pdbid1, chain1),
    #    "{0}.{1}".format(pdbid2, chain2),
    #    alignment=os.path.join(dali_dir, alignment_basename),
    #    translation=os.path.join(dali_dir, translated_basename),
    #    superposed=os.path.join(dali_dir, superposed_basename),
    #    pdb1=pdb1,
    #    pdb2=pdb2
    #)

