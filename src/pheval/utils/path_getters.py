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
        self.version = version.replace(".", "_")
        self.directory_path = None

    @property
    def tool_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}_{self.version}_{Path(self.input_dir).name}"
        )

    @tool_dir.setter
    def tool_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def batch_file_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}_{self.version}_{Path(self.input_dir).name}/{self.tool}_batch_files"
        )

    @batch_file_dir.setter
    def batch_file_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def testdata_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}_{self.version}_{Path(self.input_dir).name}/{Path(self.testdata_dir).name}_results"
        )

    @testdata_results_dir.setter
    def testdata_results_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def tool_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}_{self.version}_{Path(self.input_dir).name}/"
            f"{Path(self.testdata_dir).name}_results/{self.tool}_results"
        )

    @tool_results_dir.setter
    def tool_results_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def pheval_gene_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}_{self.version}_{Path(self.input_dir).name}/{Path(self.testdata_dir).name}_results/"
            f"{PHEVAL_GENE_RESULTS_DIR}"
        )

    @pheval_gene_results_dir.setter
    def pheval_gene_results_dir(self, directory_path):
        self.directory_path = directory_path

    @property
    def pheval_variant_results_dir(self):
        return Path(self.output_dir).joinpath(
            f"{self.tool}_{self.version}_{Path(self.input_dir).name}/{Path(self.testdata_dir).name}_results/"
            f"{PHEVAL_VARIANT_RESULTS_DIR}"
        )

    @pheval_variant_results_dir.setter
    def pheval_variant_results_dir(self, directory_path):
        self.directory_path = directory_path

    def build_directory_structure(self, phenotype_only: bool):
        self.tool_dir.mkdir(exist_ok=True, parents=True)
        self.batch_file_dir.mkdir(exist_ok=True, parents=True)
        self.testdata_results_dir.mkdir(parents=True, exist_ok=True)
        self.tool_results_dir.mkdir(parents=True, exist_ok=True)
        self.pheval_gene_results_dir.mkdir(parents=True, exist_ok=True)
        if not phenotype_only:
            self.pheval_variant_results_dir.mkdir(parents=True, exist_ok=True)

