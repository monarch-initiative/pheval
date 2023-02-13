"""CLI Test """
import logging
import unittest

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

    def test_scramble_semsim(self):
        """test_scramble_semsim"""
        result = self.runner.invoke(
            run, ["-i", "./", "-t", "./", "-r", "defaultphevalrunner", "-o", "./"]
        )
        err = result.stderr
        self.assertEqual(None, result.exception)
        logging.info("ERR=%s", err)
        exit_code = result.exit_code
        self.assertEqual(0, exit_code)

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
