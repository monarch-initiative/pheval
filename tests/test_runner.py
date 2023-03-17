import unittest
from pathlib import Path

from pheval.runners.runner import DefaultPhEvalRunner


class TestDefaultPhEvalRunner(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.pheval_runner = DefaultPhEvalRunner(
            input_dir="./tests/input_dir/configs/default/",
            testdata_dir="/corpora/corpus1/default",
            output_dir="./exomiser-13.2.0/default-corpus1-default",
            version="1.0.0",
            config_file=None,
            tmp_dir=None,
        )

    def test__get_tool(self):
        self.assertEqual(self.pheval_runner._get_tool(), "defaultrunner")

    def test__get_phenotype_only(self):
        self.assertEqual(self.pheval_runner._get_phenotype_only(), False)

    def test_runner_input_commands_dir(self):
        self.assertEqual(
            self.pheval_runner.tool_input_commands_dir,
            Path("/results/defaultrunner-1.0.0/tool_input_commands"),
        )

    def test_corpus_variant_dir(self):
        self.assertEqual(
            self.pheval_runner.corpus_variant_dir,
            Path("/results/defaultrunner-1.0.0/default-corpus1-default"),
        )

    def test_runner_results_dir(self):
        self.assertEqual(
            self.pheval_runner.runner_results_dir,
            Path("/results/defaultrunner-1.0.0/default-corpus1-default/results"),
        )

    def test_pheval_gene_results_dir(self):
        self.assertEqual(
            self.pheval_runner.pheval_gene_results_dir,
            Path("/results/defaultrunner-1.0.0/default-corpus1-default/pheval_gene_results"),
        )

    def test_pheval_variant_results_dir(self):
        self.assertEqual(
            self.pheval_runner.pheval_variant_results_dir,
            Path("/results/defaultrunner-1.0.0/default-corpus1-default/pheval_variant_results"),
        )
