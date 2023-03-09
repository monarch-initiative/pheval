from pathlib import Path

from pheval.utils.constants import PHEVAL_GENE_RESULTS_DIR, PHEVAL_VARIANT_RESULTS_DIR


class OutputDirectoryStructure:
    def __init__(
            self, output_dir: Path, input_dir: Path, testdata_dir: Path, tool: str, version: str
    ):
        self.output_dir = output_dir
        self.input_dir = input_dir
        self.testdata_dir = testdata_dir
        self.tool = tool
        self.version = version
        self.directory_path = None

    @property
    def runner_version_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}-{self.version}"
        )

    @runner_version_dir.setter
    def runner_version_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def runner_input_commands_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}-{self.version}/{self.tool}_batch_files"
        )

    @runner_input_commands_dir.setter
    def runner_input_commands_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def corpus_variant_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}-{self.version}/{Path(self.input_dir.name)}-{Path(self.testdata_dir).name}-"
            f"{Path(self.testdata_dir).parents[0]}"
        )

    @corpus_variant_dir.setter
    def corpus_variant_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def runner_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}-{self.version}/{Path(self.input_dir.name)}-{Path(self.testdata_dir).name}-"
            f"{Path(self.testdata_dir).parents[0]}/results"
        )

    @runner_results_dir.setter
    def runner_results_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def pheval_gene_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}-{self.version}/{Path(self.input_dir.name)}-{Path(self.testdata_dir).name}-"
            f"{Path(self.testdata_dir).parents[0]}/"
            f"{PHEVAL_GENE_RESULTS_DIR}"
        )

    @pheval_gene_results_dir.setter
    def pheval_gene_results_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def pheval_variant_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}-{self.version}/{Path(self.input_dir.name)}-{Path(self.testdata_dir).name}-"
            f"{Path(self.testdata_dir).parents[0]}/"
            f"{PHEVAL_VARIANT_RESULTS_DIR}"
        )

    @pheval_variant_results_dir.setter
    def pheval_variant_results_dir(self, directory_path):
        self.directory_path = directory_path

    def build_directory_structure(self, phenotype_only: bool):
        self.runner_version_dir.mkdir(exist_ok=True, parents=True)
        self.runner_input_commands_dir.mkdir(exist_ok=True, parents=True)
        self.corpus_variant_dir.mkdir(parents=True, exist_ok=True)
        self.runner_results_dir.mkdir(parents=True, exist_ok=True)
        self.pheval_gene_results_dir.mkdir(parents=True, exist_ok=True)
        if not phenotype_only:
            self.pheval_variant_results_dir.mkdir(parents=True, exist_ok=True)
