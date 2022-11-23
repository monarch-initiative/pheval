"""Runners Module"""
from abc import ABC, abstractmethod
from dataclasses import dataclass

import click


@dataclass
class PhEvalRunner(ABC):
    """PhEvalRunner Class"""

    inputdir: click.Path
    testdatadir: click.Path
    tmpdir: click.Path
    outputdir: click.Path
    config: click.Path

    @abstractmethod
    def prepare(self) -> str:
        """prepare"""

    @abstractmethod
    def run(self):
        """run"""

    @abstractmethod
    def post_process(self):
        """post_process"""


class DefaultPhEvalRunner(PhEvalRunner):
    """DefaultPhEvalRunner

    Args:
        PhEvalRunner (PhEvalRunner): Abstract PhEvalRunnerClass
    """

    inputdir: click.Path
    testdatadir: click.Path
    tmpdir: click.Path
    outputdir: click.Path
    config: click.Path

    def prepare(self):
        print("preparing")

    def run(self):
        print("running")

    def post_process(self):
        print("post processing")
