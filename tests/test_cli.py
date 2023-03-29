"""CLI Test """
import logging
import shutil
import tempfile
import unittest
from pathlib import Path

from click.testing import CliRunner

from pheval.cli_pheval import run
from pheval.cli_pheval_utils import semsim_comparison


class TestCommandLineInterface(unittest.TestCase):
    """
    Tests all command-line subcommands
    """

    def setUp(self) -> None:
        runner = CliRunner(mix_stderr=False)
        self.runner = runner
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_cli_runner(self):
        """test_cli_runner"""
        result = self.runner.invoke(
            run,
            [
                "-i",
                "./tests/input_dir/configs/default",
                "-t",
                "./testdata/phenopackets/lirical",
                "-r",
                "defaultphevalrunner",
                "-o",
                self.test_dir,
                "-v",
                "1.0.0",
            ],
        )
        err = result.stderr
        self.assertEqual(None, result.exception)
        logging.info("ERR=%s", err)
        exit_code = result.exit_code
        self.assertEqual(0, exit_code)
        self.assertTrue(Path(self.test_dir).joinpath("pheval_gene_results").exists())
        self.assertTrue(Path(self.test_dir).joinpath("pheval_variant_results").exists())
        self.assertTrue(Path(self.test_dir).joinpath("raw_results").exists())
        self.assertTrue(Path(self.test_dir).joinpath("tool_input_commands").exists())

    def test_semsim_heatmap(self):
        """test_semsim_heatmap"""
        semsim_left = "./testdata/semsim/hp-mp.semsim.tsv"
        semsim_right = "./testdata/semsim/hp-mp2.semsim.tsv"
        result = self.runner.invoke(
            semsim_comparison,
            [
                "--semsim-left",
                semsim_left,
                "--semsim-right",
                semsim_right,
                "-s",
                "jaccard_similarity",
                "-a",
                "heatmap",
            ],
        )
        err = result.stderr
        self.assertEqual(None, result.exception)
        logging.info("ERR=%s", err)
        self.assertEqual(0, result.exit_code)

    def test_semsim_heatmap_invalid_col(self):
        """test_semsim_heatmap"""
        semsim_left = "./testdata/semsim/hp-mp.semsim.tsv"
        semsim_right = "./testdata/semsim/hp-mp2.semsim.tsv"
        result = self.runner.invoke(
            semsim_comparison,
            [
                "--semsim-left",
                semsim_left,
                "--semsim-right",
                semsim_right,
                "-s",
                "invalid_col",
                "-a",
                "heatmap",
            ],
        )
        self.assertIn("Invalid value for '--score-column'", str(result.stderr))
        logging.info("ERR=%s", result.exception)
        self.assertEqual(2, result.exit_code)

    def test_semsim_heatmap_invalid_file(self):
        """test_semsim_heatmap"""
        semsim_left = "./testdata/semsim/hp-mpx.semsim.tsv"
        semsim_right = "./testdata/semsim/hp-mp2.semsim.tsv"
        result = self.runner.invoke(
            semsim_comparison,
            [
                "--semsim-left",
                semsim_left,
                "--semsim-right",
                semsim_right,
                "-s",
                "jaccard_similarity",
                "-a",
                "heatmap",
            ],
        )
        self.assertEqual(f"File {semsim_left} not found", str(result.exception))
        logging.info("ERR=%s", result.exception)
        self.assertEqual(1, result.exit_code)

    def test_semsim_heatmap_invalid_equal_file(self):
        """test_semsim_heatmap"""
        semsim_left = "./testdata/semsim/hp-mp.semsim.tsv"
        semsim_right = "./testdata/semsim/hp-mp.semsim.tsv"
        result = self.runner.invoke(
            semsim_comparison,
            [
                "--semsim-left",
                semsim_left,
                "--semsim-right",
                semsim_right,
                "-s",
                "jaccard_similarity",
                "-a",
                "heatmap",
            ],
        )
        errmsg = "Semantic similarity profiles are equal. Make sure you have selected different files to analyze"
        self.assertEqual(errmsg, str(result.exception))
        logging.info("ERR=%s", result.exception)
        self.assertEqual(1, result.exit_code)
