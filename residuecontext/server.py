#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import time

from residuecontext.electrostatics import get_electrostatics_grid

try:
    from .builder import (
        ResidueContextParams,
        ContextBuilder,
        create_grid_histogram,
    )
    from .helpers import (
        get_pdb_selection,
        create_context_from_pdb_code,
        chain_context_to_jsonable,
    )
    from .pdbfiles import (
        get_pdb_selection,
        read_alignment_file,
    )
except (ValueError, SystemError):
    from builder import (
        ResidueContextParams,
        ContextBuilder,
        create_grid_histogram,
    )
    from helpers import (
        get_pdb_selection,
        create_context_from_pdb_code,
        chain_context_to_jsonable,
    )
    from pdbfiles import (
        get_pdb_selection,
        read_alignment_file,
    )

from flask import (
    abort,
    Flask,
    g,
    jsonify,
    request,
    render_template,
    send_file,
)


app = Flask(__name__)


def after_this_request(f):
    if not hasattr(g, 'after_request_callbacks'):
        g.after_request_callbacks = []
    g.after_request_callbacks.append(f)
    return f


def extract_ident(ident):
    if len(ident) not in (4, 5):
        abort(400)
    pdbid = ident[:4]
    chain = ident[4:] or '_'
    return pdbid, chain


@app.route('/contexts/sterics/<ident>.json', defaults={'alignment_id': None})
@app.route('/alignments/<alignment_id>/contexts/sterics/<ident>.json')
def get_context(ident, alignment_id=None):
    pdbid, chain = extract_ident(ident)
    if chain == '_':
        chain = None

    app.logger.info("Starting download of {}".format(pdbid))
    start = time.time()
    data = get_pdb_selection(pdbid, chain, alignment_id=alignment_id)
    download_time = time.time() - start
    app.logger.info("Download took {:.4f}sec".format(download_time))

    if data is None:
        raise abort(404)

    start = time.time()
    app.logger.info("Starting residue context generation of {}".format(ident))
    params = ResidueContextParams(
        source=data,
        chain=chain,
        name=pdbid
    )
    builder = ContextBuilder(params)
    context = builder.run()
    generate_time = time.time() - start
    app.logger.info("Generation took {:.4f}sec".format(generate_time))

    app.logger.info("Starting production of JSON results for {}".format(ident))
    start = time.time()
    result = chain_context_to_jsonable(context)
    render_time = time.time() - start
    app.logger.info("Production took {:.4f}sec".format(generate_time))

    result['timing'] = {
        'download': round(download_time, 4),
        'generate': round(generate_time, 4),
        'render': round(render_time, 4),
    }

    @after_this_request
    def report_timing():
        transmit_time = time.time() - start
        app.logger.info("Transmission took {:.4f}sec".format(transmit_time))

    app.logger.info("Starting transmission results for {}".format(ident))
    start = time.time()
    return jsonify(result)


@app.route('/contexts/phi/<ident>.json', defaults={'alignment_id': None})
@app.route('/alignments/<alignment_id>/contexts/phi/<ident>.json')
def get_phi_context(ident, alignment_id=None):
    pdbid, chain = extract_ident(ident)
    if chain == '_':
        chain = None
    coords = map(float, request.args['coords'].split())
    app.logger.info("Starting generation of {}".format(pdbid))
    start = time.time()
    phi = get_electrostatics_grid(pdbid, chain=chain, alignment_id=alignment_id)
    download_time = time.time() - start
    app.logger.info("Generation took {:.4f}sec".format(download_time))

    if phi is None:
        abort(404)

    start = time.time()
    app.logger.info("Starting grid context generation of {}".format(ident))
    grid_histograms = create_grid_histogram(phi, coords, extents=(-7, 7))
    generate_time = time.time() - start
    app.logger.info("Histogramming took {:.4f}sec".format(generate_time))

    app.logger.info("Starting production of JSON results for {}".format(ident))
    start = time.time()
    result = {
        'histograms': grid_histograms.tolist()
    }
    render_time = time.time() - start
    app.logger.info("Production took {:.4f}sec".format(generate_time))

    result['timing'] = {
        'download': round(download_time, 4),
        'generate': round(generate_time, 4),
        'render': round(render_time, 4),
    }

    @after_this_request
    def report_timing():
        transmit_time = time.time() - start
        app.logger.info("Transmission took {:.4f}sec".format(transmit_time))

    app.logger.info("Starting transmission results for {}".format(ident))
    start = time.time()
    return jsonify(result)


@app.route('/pdb/<pdb>.pdb', defaults={'alignment_id': None})
@app.route('/alignments/<alignment_id>/pdb/<pdb>.pdb')
def get_pdb(pdb, alignment_id=None):
    return send_file(get_pdb_selection(pdb, chain=None, model=None, alignment_id=alignment_id))


@app.route('/alignments/<alignment_id>/alignment.json')
def get_alignment_blocks(alignment_id):
    alignment = read_alignment_file(alignment_id)
    return jsonify(**alignment)


@app.route('/')
def index():
    return render_template('client.html')


if __name__ == "__main__":
    import logging
    app.logger.setLevel(logging.INFO)
    app.run(debug=True)