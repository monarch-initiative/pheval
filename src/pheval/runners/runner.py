"""Runners Module"""
import timeit
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

from pheval.config_parser import parse_input_dir_config
from pheval.utils.constants import (
    PHEVAL_GENE_RESULTS_DIR,
    PHEVAL_VARIANT_RESULTS_DIR,
    RUNNER_INPUT_COMMANDS_DIR,
    TOOL_RESULTS_DIR,
)
from pheval.construct_metadata import create_run_metadata


@dataclass
class PhEvalRunner(ABC):
    """PhEvalRunner Class"""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str
    directory_path = None
    input_dir_config = None
    start_time = timeit.default_timer()

    def __post_init__(self):
        self.input_dir_config = parse_input_dir_config(self.input_dir)

    def _get_tool(self):
        return self.input_dir_config.tool

    def _get_phenotype_only(self):
        return self.input_dir_config.phenotype_only

    @property
    def runner_version_dir(self):
        return Path(self.output_dir).joinpath(f"{self._get_tool()}-{self.version}")

    @runner_version_dir.setter
    def runner_version_dir(self, directory_path):
        self.directory_path = Path(directory_path)

    @property
    def runner_input_commands_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self._get_tool()}-{self.version}/{RUNNER_INPUT_COMMANDS_DIR}"
        )

    @runner_input_commands_dir.setter
    def runner_input_commands_dir(self, directory_path):
        self.directory_path = Path(directory_path)

    @property
    def corpus_variant_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self._get_tool()}-{self.version}/{Path(self.input_dir).name}-"
            f"{Path(self.testdata_dir).parent.name}-"
            f"{Path(self.testdata_dir).name}"
        )

    @corpus_variant_dir.setter
    def corpus_variant_dir(self, directory_path):
        self.directory_path = Path(directory_path)

    @property
    def runner_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self._get_tool()}-{self.version}/{Path(self.input_dir).name}-"
            f"{Path(self.testdata_dir).parent.name}-"
            f"{Path(self.testdata_dir).name}/{TOOL_RESULTS_DIR}"
        )

    @runner_results_dir.setter
    def runner_results_dir(self, directory_path):
        self.directory_path = Path(directory_path)

    @property
    def pheval_gene_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self._get_tool()}-{self.version}/{Path(self.input_dir).name}-"
            f"{Path(self.testdata_dir).parent.name}-"
            f"{Path(self.testdata_dir).name}/"
            f"{PHEVAL_GENE_RESULTS_DIR}"
        )

    @pheval_gene_results_dir.setter
    def pheval_gene_results_dir(self, directory_path):
        self.directory_path = Path(directory_path)

    @property
    def pheval_variant_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self._get_tool()}-{self.version}/{Path(self.input_dir).name}-"
            f"{Path(self.testdata_dir).parent.name}-"
            f"{Path(self.testdata_dir).name}/"
            f"{PHEVAL_VARIANT_RESULTS_DIR}"
        )

    @pheval_variant_results_dir.setter
    def pheval_variant_results_dir(self, directory_path):
        self.directory_path = Path(directory_path)

    def build_output_directory_structure(self):
        """build output directory structure"""
        self.runner_version_dir.mkdir(exist_ok=True, parents=True)
        self.runner_input_commands_dir.mkdir(exist_ok=True, parents=True)
        self.corpus_variant_dir.mkdir(parents=True, exist_ok=True)
        self.runner_results_dir.mkdir(parents=True, exist_ok=True)
        self.pheval_gene_results_dir.mkdir(parents=True, exist_ok=True)
        if not self._get_phenotype_only():
            self.pheval_variant_results_dir.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def prepare(self) -> str:
        """prepare"""

    @abstractmethod
    def run(self):
        """run"""

    @abstractmethod
    def post_process(self):
        """post_process"""

    def create_metadata(self):
        """create metadata for run"""
        return create_run_metadata(
            input_dir_config=self.input_dir_config,
            input_dir=self.input_dir,
            tool_version=self.version,
            corpus=self.testdata_dir,
            start_time=self.start_time,
            corpus_variant_path=self.corpus_variant_dir,
        )


class DefaultPhEvalRunner(PhEvalRunner):
    """DefaultPhEvalRunner

    Args:
        PhEvalRunner (PhEvalRunner): Abstract PhEvalRunnerClass
    """

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        print("preparing")

    def run(self):
        print("running")

    def post_process(self):
        print("post processing")
