#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
import json

import werkzeug.security

from celery import Celery

from residuecontext.config import (
    ALIGNMENT_JOB_DIR,
    ALIGNMENT_DATA_FILE,
)
from residuecontext.pdbfiles import (
    get_pdb_selection,
    extract_ident,
    ident_to_biojava,
    ident_to_rcontext,
)
from residuecontext.biojava_alignment import (
    BIOJAVA_CE_CLASS,
    BIOJAVA_FATCAT_CLASS,
    run_biojava_alignment,
)
from residuecontext.rcontext_alignment import run_residuecontext_alignment
from residuecontext.dalilite_alignment import run_dalilite_alignment
from residuecontext.electrostatics import get_electrostatics_grid
from residuecontext.vanderwaals import get_vanderderwaals_grids


queue = Celery('tasks',
               backend='db+sqlite:///celerydb.sqlite',
               broker='sqla+sqlite:///celerydb.sqlite')
queue.conf['CELERY_RESULT_DB_SHORT_LIVED_SESSIONS'] = True

@queue.task
def run_alignment_comparisons(number, code1, code2):
    working_dir = werkzeug.security.safe_join(ALIGNMENT_JOB_DIR, str(number))

    if os.path.exists(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)

    pdbid1, chain1 = extract_ident(code1)
    pdbid2, chain2 = extract_ident(code2)

    with open(os.path.join(working_dir, ALIGNMENT_DATA_FILE), 'w') as f:
        json.dump({
            'structure1': {
                'ident': code1,
                'pdb': pdbid1,
                'chain': chain1
            },
            'structure2': {
                'ident': code2,
                'pdb': pdbid2,
                'chain': chain2
            }
        }, f)

    pdb1src = get_pdb_selection(pdbid1, chain=chain1)
    pdb2src = get_pdb_selection(pdbid2, chain=chain2)

    pdb1 = os.path.join(working_dir, "{0}.pdb".format(pdbid1))
    pdb2 = os.path.join(working_dir, "{0}.pdb".format(pdbid2))

    shutil.copy(pdb1src, pdb1)
    shutil.copy(pdb2src, pdb2)

    rcontext_dir = os.path.join(working_dir, "rcontext")
    ce_dir = os.path.join(working_dir, "ce")
    fatcat_dir = os.path.join(working_dir, "fatcat")
    dali_dir = os.path.join(working_dir, "dali")

    alignment_basename = 'alignment.txt'
    translated_basename = 'transformed-{0}.pdb'.format(pdbid2)
    superposed_basename = 'superposed.pdb'
    transform_basename = '{0}.pdb-transform.txt'.format(pdbid2)

    biojava_code1 = ident_to_biojava(pdbid1, chain1)
    biojava_code2 = ident_to_biojava(pdbid2, chain2)

    os.mkdir(ce_dir)
    run_biojava_alignment(
        BIOJAVA_CE_CLASS,
        biojava_code1,
        biojava_code2,
        alignment=os.path.join(ce_dir, alignment_basename),
        transformed=os.path.join(ce_dir, translated_basename),
        superposed=os.path.join(ce_dir, superposed_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )
    os.symlink(pdb1, os.path.join(ce_dir, "{0}.pdb".format(pdbid1)))
    os.symlink(pdb2, os.path.join(ce_dir, "{0}.pdb".format(pdbid2)))

    os.mkdir(fatcat_dir)
    run_biojava_alignment(
        BIOJAVA_FATCAT_CLASS,
        biojava_code1,
        biojava_code2,
        alignment=os.path.join(fatcat_dir, alignment_basename),
        transformed=os.path.join(fatcat_dir, translated_basename),
        superposed=os.path.join(fatcat_dir, superposed_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )
    os.symlink(pdb1, os.path.join(fatcat_dir, "{0}.pdb".format(pdbid1)))
    os.symlink(pdb2, os.path.join(fatcat_dir, "{0}.pdb".format(pdbid2)))

    os.mkdir(rcontext_dir)
    run_residuecontext_alignment(
        ident_to_rcontext(pdbid1, chain1),
        ident_to_rcontext(pdbid2, chain2),
        alignment=os.path.join(rcontext_dir, alignment_basename),
        transformed=os.path.join(rcontext_dir, translated_basename),
        translation=os.path.join(rcontext_dir, transform_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )
    os.symlink(pdb1, os.path.join(rcontext_dir, "{0}.pdb".format(pdbid1)))
    os.symlink(pdb2, os.path.join(rcontext_dir, "{0}.pdb".format(pdbid2)))

    os.mkdir(dali_dir)
    run_dalilite_alignment(
        ident_to_rcontext(pdbid1, chain1),
        ident_to_rcontext(pdbid2, chain2),
        alignment=os.path.join(dali_dir, alignment_basename),
        translation=os.path.join(dali_dir, translated_basename),
        transformed=os.path.join(dali_dir, transform_basename),
        pdb1=pdb1,
        pdb2=pdb2
    )
    os.symlink(pdb1, os.path.join(dali_dir, "{0}.pdb".format(pdbid1)))
    os.symlink(pdb2, os.path.join(dali_dir, "{0}.pdb".format(pdbid2)))

