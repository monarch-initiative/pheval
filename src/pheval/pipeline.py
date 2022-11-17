import click


class PhEvalRunner:
    def __init__(
        self,
        inputdir: click.Path,
        testdatadir: click.Path,
        tmpdir: click.Path,
        outputdir: click.Path,
        config: click.Path,
    ):
        self.inputdir = inputdir
        self.testdatadir = testdatadir
        self.tmpdir = tmpdir
        self.outputdir = outputdir
        self.config = config

    def prepare(self):
        raise NotImplementedError()

    def run(self):
        raise NotImplementedError()

    def post_process(self):
        raise NotImplementedError()


class DefaultPhEvalRunner(PhEvalRunner):
    def __init__(
        self,
        inputdir: click.Path,
        testdatadir: click.Path,
        tmpdir: click.Path,
        outputdir: click.Path,
        config: click.Path,
    ):
        pass

    def prepare(self):
        print("preparing")

    def run(self):
        print("running")

    def post_process(self):
        print("post processing")
