import click

from pheval.pipeline import PhEvalRunner


class ExomiserPhEvalRunner(PhEvalRunner):
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
