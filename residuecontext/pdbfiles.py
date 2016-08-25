#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import contextlib
import gzip
import os

import werkzeug.security

from Bio.PDB.PDBParser import PDBParser as _PDBParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBIO import (
    PDBIO,
    Select,
)

import numpy as np

from residuecontext.config import (
    STATIC_PDB_DIR,
    ALIGNMENT_JOB_DIR,
)

STATIC_PDB_LIST = PDBList(pdb=STATIC_PDB_DIR)


def get_pdb_id_from_filename(filename):
    pdbid = os.path.basename(filename).split('.')[0]
    return pdbid


def extract_ident(ident):
    if len(ident) not in (4, 5):
        return None, None
    pdbid = ident[:4]
    chain = ident[4:] or '_'
    return pdbid, chain


def split_ext_gz(path):
    path, ext = os.path.splitext(path)

    if ext == '.gz':
        path, ext = os.path.splitext(path)
        ext += '.gz'

    return path, ext


def chunks_of(it, n, fill=None, wrap=tuple):
    chunk = []
    for el in it:
        chunk.append(el)
        if len(chunk) == n:
            yield wrap(chunk)
            chunk = []
    if len(chunk) > 0:
        if fill is not None:
            while len(chunk) < n:
                chunk.append(fill)
        yield wrap(chunk)


class PDBParser:
    NAME_UNKONWN = "Unkown"

    def __init__(self):
        self._parser = _PDBParser()

    @contextlib.contextmanager
    def open(self, source):
        if isinstance(source, str):
            if source.endswith('.gz'):
                ctx = gzip.open(source, 'r')
            else:
                ctx = open(source, 'r')
        else:
            ctx = source

        yield ctx

    @staticmethod
    def guess_name(source, default=NAME_UNKONWN):
        if isinstance(source, str):
            name = get_pdb_id_from_filename(source)
        elif hasattr(source, 'name'):
            name = get_pdb_id_from_filename(source.name)
        elif hasattr(source, 'url'):
            name = get_pdb_id_from_filename(source.url)
        else:
            name = default

        return name

    def load_structure(self, source, name=None):
        if name is None:
            name = self.guess_name(source)
        with self.open(source) as f:
            structure = self._parser.get_structure(name, f)
        return structure


class ChainSelect(Select):

    def __init__(self, chain=None, model=0):
        super(ChainSelect, self).__init__()
        self.chain = chain
        self.model = model

    def accept_model(self, model):

        if self.model is None or self.model == model.get_id():
            return 1
        else:
            return 0

    def accept_chain(self, chain):
        if self.chain is None or self.chain == chain.get_id():
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        res_id = residue.get_id()
        res_code = res_id[0]
        ins_code = res_id[2]

        if res_code == 'W' or res_code.startswith('H_'):
            return 0
        elif ins_code in ('_', 'A', ' '):
            return 1


pdb_parser = PDBParser()
pdb_writer = PDBIO(use_model_flag=0)  # Reduce doesn't like "model 0"


def get_pdb_selection(code, chain=None, model=0, alignment_id=None, root=None, target=None):
    code = os.path.basename(code)

    if alignment_id is not None:
        root = root or ALIGNMENT_JOB_DIR
        root = werkzeug.security.safe_join(root, str(alignment_id))

    if root is None:
        path = STATIC_PDB_LIST.retrieve_pdb_file(code)
    elif os.path.isdir(root):
        path = werkzeug.security.safe_join(root, code + '.pdb')
    else:
        path = root

    if root is not None:
        transformed_path = werkzeug.security.safe_join(root, 'transformed-' + code + '.pdb')
    else:
        transformed_path = None
    if transformed_path is not None and os.path.exists(transformed_path):
        path = transformed_path
        transforms = None
    elif os.path.exists(path + '-transform.txt'):
        transforms = path + '-transform.txt'
    else:
        transforms = None

    structure = pdb_parser.load_structure(path, name=code)

    if transforms:
        array = np.loadtxt(transforms)
        translation, rotation = array[0], array[1:4]
        structure.transform(rotation, translation)
        path = transformed_path

    selector = ChainSelect(chain=chain, model=model)

    if model is not None or chain is not None:
        pdb_base, ext = split_ext_gz(path)
        selection_path = pdb_base

        if model is None:
            selection_path += '_'
        else:
            selection_path += str(model)

        if chain is None:
            selection_path += '_'
        else:
            selection_path += str(chain)

        selection_path += '.pdb'
    else:
        selection_path = path

    if target is not None:
        if os.path.isdir(target):
            selection_path = os.path.join(target, os.path.basename(selection_path))
        else:
            selection_path = target

    if not os.path.exists(selection_path):
        pdb_writer.set_structure(structure)
        pdb_writer.save(selection_path, select=selector)

    return selection_path


def load_alignment_file(io):
    chunks = []
    buf = []

    title = None

    for line in io:
        line = line.rstrip()

        if title is None:
            title = line
        elif len(buf) == 0 and line.startswith('Structure 1'):
            buf.append(line)
        elif len(buf) == 1 and line[:28].strip() == '':
            buf.append(line)
        elif len(buf) == 2 and line.startswith('Structure 2'):
            buf.append(line)
            chunks.append(buf)
            buf = []
        elif len(chunks) > 0 and len(buf) == 0 and len(line) > 0:
            break

    title_parts = title.split()
    structure1_id = title_parts[4] + title_parts[6]
    structure2_id = title_parts[8] + title_parts[10]

    alignment = {
        structure1_id: [],
        structure2_id: []
    }

    for structure1_line, blocks, structure2_line in chunks:
        line_start1 = int(structure1_line[19:27].strip())
        line_start2 = int(structure2_line[19:27].strip())

        block_nums = [int(num) for num in chunks_of(blocks[28:129], 2, wrap=' '.join)]
        line1_rescodes = structure1_line[28:129].split()
        line2_rescodes = structure2_line[28:129].split()

        joined_block = zip(
            block_nums,
            enumerate(line1_rescodes, start=line_start1),
            enumerate(line2_rescodes, start=line_start2)
        )

        for block_num, (resnum1, rescode1), (resnum2, rescode2) in joined_block:
            if block_num > len(alignment[structure1_id]):
                alignment[structure1_id].append({
                    'block_num': block_num,
                    'indices': [],
                    'sequence': [],
                })

            if block_num > len(alignment[structure2_id]):
                alignment[structure2_id].append({
                    'block_num': block_num,
                    'indices': [],
                    'sequence': [],
                })

            alignment[structure1_id][block_num-1]['indices'].append(resnum1)
            alignment[structure1_id][block_num-1]['sequence'].append(rescode1)

            alignment[structure2_id][block_num-1]['indices'].append(resnum2)
            alignment[structure2_id][block_num-1]['sequence'].append(rescode2)

    return alignment


def read_alignment_file(alignment_id, kind=None, root=ALIGNMENT_JOB_DIR):
    alignment_dir = werkzeug.security.safe_join(root, str(alignment_id))
    if kind is not None:
        alignment_dir = werkzeug.security.safe_join(alignment_dir, kind)
    alignment_file = os.path.join(alignment_dir, 'alignment.txt')

    with open(alignment_file) as f:
        blocks = load_alignment_file(f)

    return blocks
