from .config import ResidueContextConfig


class ContextSpec:
    def __init__(self, env, pdb, chain):
        self.env = env
        self.pdb = pdb
        self.chain = chain
        self.model = 0


class Environment:
    def __init__(self, config=None):
        self.config = config or ResidueContextConfig()
        self.pdb_parser = None

