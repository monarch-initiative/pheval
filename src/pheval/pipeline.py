"""pipeline"""
from dataclasses import dataclass

import click


@dataclass
class PhEvalRunner:
    """_summary_

    Raises:
        NotImplementedError: _description_
        NotImplementedError: _description_
        NotImplementedError: _description_
    """

    inputdir: click.Path
    testdatadir: click.Path
    tmpdir: click.Path
    outputdir: click.Path
    config: click.Path

    """_summary_"""

    def prepare(self) -> str:
        """prepare
        Raises:
            NotImplementedError:
        """
        raise NotImplementedError

    def run(self):
        """run
        Raises:
            NotImplementedError:
        """
        raise NotImplementedError

    def post_process(self):
        """post_process
        Raises:
            NotImplementedError:
        """
        raise NotImplementedError


@dataclass
class DefaultPhEvalRunner(PhEvalRunner):
    """DefaultPhEvalRunner
    Args:
        PhEvalRunner (_type_): _description_
    """

    inputdir: click.Path
    testdatadir: click.Path
    tmpdir: click.Path
    outputdir: click.Path
    config: click.Path

    def prepare(self) -> str:
        print("preparing")
        return "default"

    def run(self):
        print("running")

    def post_process(self):
        print("post processing")
