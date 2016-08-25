def get_parser(config=None):
    parser = PDBParser()
    return parser


def init_env(env):
    env.pdb_parser = get_parser(env.config)


