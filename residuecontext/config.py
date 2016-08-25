import os

ALIGNMENT_JOB_DIR = os.path.join(os.path.dirname(__file__), 'static', 'job-data')
ALIGNMENT_DATA_FILE = 'alignment-params.json'
STATIC_PDB_DIR = os.path.join(os.path.dirname(__file__), 'static', 'pdb')

REDUCE_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'reduce')
QNIFFT_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'qnifft')
BIOJAVA_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'biojava')
RESIDUE_CONTEXT_DIR = os.path.join(os.path.dirname(__file__), 'scripts', 'rcontext')

JAVA_RUNTIME = "java"
