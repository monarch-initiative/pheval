"""Runners Module"""
from abc import ABC, abstractmethod
from dataclasses import dataclass

import click


@dataclass
class PhEvalRunner(ABC):
    """PhEvalRunner Class"""

    input_dir: click.Path
    testdata_dir: click.Path
    tmp_dir: click.Path
    output_dir: click.Path
    config_file: click.Path

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

    input_dir: click.Path
    testdata_dir: click.Path
    tmp_dir: click.Path
    output_dir: click.Path
    config_file: click.Path

    def prepare(self):
        print("preparing")

    def run(self):
        print("running")

    def post_process(self):
        print("post processing")
