class PhEvalRunner:
    def __init__(self):
        """PhEvalRunner"""

    def prepare(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError

    def post_process(self):
        raise NotImplementedError
