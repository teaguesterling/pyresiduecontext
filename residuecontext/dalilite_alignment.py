#!/usr/bin/env python
# residuecontext.algorithm.AlgorithmWrapper
# 1AZP.pdb 1C8C.pdb A A working -v --verbhist
# -DPATH=$PATH:../src/external/emd/emd/runEMD
from __future__ import absolute_import, division, print_function

import logging
import os
import shutil
import tempfile
from six.moves.urllib import parse

import numpy as np
from BeautifulSoup import BeautifulSoup

import sh


from residuecontext.config import (
    DALILITE_DIR,
)
from residuecontext.pdbfiles import split_ext_gz, get_pdb_selection


DALILITE_ENV = os.environ.copy()
DALILITE_ENV['PATH'] += os.pathsep + DALILITE_DIR

ALIGNMENT_TEMP_PDB = "aligned.pdb"
ALIGNMENT_OUTPUT = "workingalignment.txt"

_dalilite_cmd = sh.Command(os.path.join(DALILITE_DIR, 'DaliLite'))
_dalilite_convert = _dalilite_cmd.bake('-readbrk')
_dalilite_align = _dalilite_cmd.bake('-align')
_dalilite_format = _dalilite_cmd.bake('-format')


def parse_dali_html(path):
    with open(path) as f:
        content = f.read()

    soup = BeautifulSoup(content)
    anchors = soup.findAll('a')
    query = anchors[0]
    ident1 = query['name']
    pdbid1 = ident1[:4]
    chain1 = ident1[4]

    id_url = parse.urlparse(anchors[4]['href'])
    id_params = parse.parse_qs(id_url.query)

    pdbid2 = id_params['pdbid'][0]
    chain2 = id_params['chainid'][0]

    transform_url = parse.urlparse(anchors[5]['href'])
    transform_params = parse.parse_qs(transform_url.query)

    rotation_raw = [float(x) for x in transform_params['u'][0].split(',')]
    translation_raw = [float(x) for x in transform_params['t'][0].split(',')]
    rotation = np.array(rotation_raw).reshape((3,3))
    translation = np.array(translation_raw)

    pre_tags = soup.findAll('pre')
    raw_alignment = pre_tags[1].string
    alignment_blocks = raw_alignment.split('\n\n')
    alignment = []

    offset1 = 1
    offset2 = 1
    number = 1
    for alignment_block in alignment_blocks:
        block_lines = alignment_block.strip().splitlines()
        if len(block_lines) < 5:
            continue
        alignment_lines = block_lines[1:4]
        alignment_data = zip(*[line[6:51] for line in alignment_lines])
        alignment_next1 = int(alignment_lines[0][66:].strip())
        alignment_next2 = int(alignment_lines[2][66:].strip())

        block = []
        for q, i, s in alignment_data:
            if q == '-' or s == '-':
                number += 1
                n = ' '
            else:
                n = i
            block.append((q, n, s))
        alignment.append((
            (offset1, offset2),
            block,
            (alignment_next1, alignment_next2)
        ))
        offset1 = alignment_next1
        offset2 = alignment_next2

    return {
        'transform': (rotation, translation),
        'ident1': '{0}.{1}'.format(pdbid1, chain1),
        'ident2': '{0}.{1}'.format(pdbid2, chain2),
        'n1': offset1,
        'n2': offset2,
        'alignment': alignment,
    }


def write_faux_dali_alignment_file(io, dali_data):
    print('DaliLight Alignment FatCat Style Alignment File', file=io)
    print('Align {ident1}.pdb {n1} with {ident2}.pdb {n2}'.format(**dali_data), file=io)
    print('', file=io)

    for (start1, start2), (_end1, _end2), block in dali_data['alignment']:
        buffer1 = ['Chain 1: {0: <4} '.format(start1)]
        buffer2 = ['              ']
        buffer3 = ['Chain 2: {0: <4} '.format(start2)]
        for a, b, c in block:
            buffer1.append(a)
            buffer2.append(b)
            buffer3.append(c)

        print(''.join(buffer1), file=io)
        print(''.join(buffer2), file=io)
        print(''.join(buffer3), file=io)
        print('', file=io)

    print('', file=io)


def run_dalilite_alignment(
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

        pdb1 = pdb1 or '{0}.pdb'.format(pdbid1)
        pdb2 = pdb2 or '{0}.pdb'.format(pdbid2)

        formatted = '{0}{1}.html'.format(pdbid1, chain1)
        dccp = '{0}{1}.dccp'.format(pdbid2, chain2)

        logging.info("Preparing PDB inputs in {0}".format(exec_dir))
        _convert_cmd = _dalilite_convert.bake(_cwd=exec_dir, _err=log, _env=DALILITE_ENV)
        for pdb, code in [(pdb1, pdbid1), (pdb2, pdbid2)]:
            entfile = "{0}.pdb".format(code.lower())
            os.symlink(pdb, os.path.join(exec_dir, entfile))
            _convert_cmd(entfile, code)

        # The PDBIDs are intentionally reversed here
        align_cmd = _dalilite_align.bake(
            "{0}{1}".format(pdbid2, chain2),
            "{0}{1}".format(pdbid1, chain1)
        )
        logging.info("Running DaliLite alignment {0!s} in {1}".format(_dalilite_align, exec_dir))
        align_cmd(
            _cwd=exec_dir,
            _err=log,
            _env=DALILITE_ENV
        )

        _format_cmd = _dalilite_format.bake(
            "{0}{1}".format(pdbid1, chain1),
            dccp,
            'list1',
            formatted
        )
        logging.info("Extracting DaliLite Results: {0!s} in {1}".format(_format_cmd, exec_dir))
        _format_cmd(
            _cwd=exec_dir,
            _err=log,
            _env=DALILITE_ENV
        )

        alignment_data = parse_dali_html(os.path.join(exec_dir, formatted))

        with open(os.path.join(exec_dir, ALIGNMENT_OUTPUT), 'w') as f:
            write_faux_dali_alignment_file(f, alignment_data)

        transform_file = os.path.join(exec_dir, "{0}.pdb-transform.txt".format(ident2[0].upper()))
        with open(transform_file) as f:
            rot, trans = alignment_data['transform']
            print('\t'.join(trans), file=f)
            print('\t'.join(rot[0]), file=f)
            print('\t'.join(rot[1]), file=f)
            print('\t'.join(rot[2]), file=f)

        pdb2_transformed = get_pdb_selection(
            pdbid2.upper(),
            chain=chain2,
            root=exec_dir,
        )

        shutil.move(
            os.path.join(exec_dir, ALIGNMENT_OUTPUT),
            alignment
        )
        shutil.move(
            os.path.join(exec_dir, "{0}.pdb-transform.txt".format(ident2[0].upper())),
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


