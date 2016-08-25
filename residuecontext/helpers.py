#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

from collections import defaultdict
import logging
import os
import tempfile

import sh


from residuecontext.builder import (
    ResidueContextParams,
    ContextBuilder,
)
from residuecontext.pdbfiles import get_pdb_selection


def create_context_from_pdb_code(code, chain):
    pdb = get_pdb_selection(code, chain)
    params = ResidueContextParams(
        source=pdb,
        chain=chain,
        name=code
    )
    builder = ContextBuilder(params)
    context = builder.run()
    return context


def chain_context_to_jsonable(context, sparse=True, inverse=True):
    data = {
        'name': context.name,
        'bin_data': dict(
            (k, v.tolist() if hasattr(v, 'tolist') else v)
            for k, v in context.bin_data.items()
        ),
        'contexts': [],
        'sparse': sparse,
        'inverse': inverse,
        'templates': {
            'vertices_per_sphere': context.surface['vertices_per_sphere'],
            'faces_per_sphere': context.surface['faces_per_sphere'],  # Should be same as bins
            'vertices': [[round(i, 3) for i in v] for v in context.surface['vertices']],
            'faces': context.surface['faces'].tolist()
        }
    }

    if inverse:
        inverses = defaultdict(lambda: defaultdict(list))

    for index, ctx in enumerate(context.contexts):
        if sparse:
            histogram = dict(
                (index, int(value))
                for index, value in enumerate(ctx.histogram)
                if value > 0
            )
        else:
            histogram = ctx.histogram.tolist()

        this_context = {
            'index': int(ctx.index),
            'name': ctx.name,
            'serial': int(ctx.serial),
            'resnum': int(ctx.resnum),
            'aacode': ctx.aacode,
            'size': float(ctx.size),
            'coords': [
                round(float(ctx.x), 3),
                round(float(ctx.y), 3),
                round(float(ctx.z), 3)
            ],
            'histogram': histogram
        }

        if inverse:
            this_context['inverse'] = inverses[index]
            for inverse_index, bin_number in enumerate(ctx.inverted.tolist()):
                if bin_number is not None:
                    inverses[inverse_index][bin_number].append(index)

        data['contexts'].append(this_context)

    return data



def grid_context_to_jsonable(context, sparse=True, inverse=True):
    data = {
        'name': context.name,
        'bin_data': {
            k: v.tolist() if hasattr(v, 'tolist') else v
            for k, v in context.bin_data.items()
        },
        'contexts': [],
        'sparse': sparse,
        'inverse': inverse,
        'templates': {
            'vertices_per_sphere': context.surface['vertices_per_sphere'],
            'faces_per_sphere': context.surface['faces_per_sphere'],  # Should be same as bins
            'vertices': [[round(i, 3) for i in v] for v in context.surface['vertices']],
            'faces': context.surface['faces'].tolist()
        }
    }

    if inverse:
        inverses = defaultdict(lambda: defaultdict(list))

    for index, ctx in enumerate(context.contexts):
        if sparse:
            histogram = [
                (index, int(value))
                for index, value in enumerate(ctx.histogram)
                if value > 0
            ]
        else:
            histogram = ctx.histogram.tolist()

        this_context = {
            'index': int(ctx.index),
            'name': ctx.name,
            'serial': int(ctx.serial),
            'aacode': ctx.aacode,
            'size': float(ctx.size),
            'coords': [
                round(float(ctx.x), 3),
                round(float(ctx.y), 3),
                round(float(ctx.z), 3)
            ],
            'histogram': histogram
        }

        if inverse:
            this_context['inverse'] = inverses[index]
            for inverse_index, bin_number in enumerate(ctx.inverted.tolist()):
                if bin_number is not None:
                    inverses[inverse_index][bin_number].append(index)

        data['contexts'].append(this_context)

    return data

