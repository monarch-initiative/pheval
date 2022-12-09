"""CLI Test """
import logging
import os
import unittest
from pathlib import Path

import pytest
from click.testing import CliRunner

from pheval.cli_pheval_utils import compare_semsim


class TestCommandLineInterface(unittest.TestCase):
    """
    Tests all command-line subcommands
    """

    def setUp(self) -> None:
        Path("./res.tsv").unlink(missing_ok=True)
        runner = CliRunner(mix_stderr=False)
        self.runner = runner

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        """inject_fixtures"""
        self._caplog = caplog  # pylint: disable=attribute-defined-outside-init

    def test_compare_semsim(self):
        """test_compare_semsim"""
        params = [
            "-l",
            "./tests/fixtures/left.tsv",
            "-r",
            "./tests/fixtures/right.tsv",
            "-o",
            "res.tsv",
            "-j",
            "GENE_SYMBOL",
            "-c",
            "MOUSE_PHENO_SCORE",
            "-c",
            "FISH_PHENO_SCORE",
        ]
        result = self.runner.invoke(
            compare_semsim,
            params,
        )
        exit_code = result.exit_code
        self.assertEqual(0, exit_code)
        self.assertTrue(os.path.exists("./res.tsv"), "Result file must exists")

    def test_compare_semsim_inexistent_column(self):
        """test_compare_semsim"""
        params = [
            "-l",
            "./tests/fixtures/left.tsv",
            "-r",
            "./tests/fixtures/right.tsv",
            "-o",
            "res.tsv",
            "-j",
            "GENE_SYMBOL",
            "-c",
            "MOUSE_PHENO_SCORE",
            "-c",
            "MOUSE_PHENO_SCOR",
        ]
        result = self.runner.invoke(
            compare_semsim,
            params,
        )
        err_log = self._caplog.records[0].message
        self.assertTrue("must be present in both left and right files" in err_log)
        logging.info("ERR=%s", err_log)
        exit_code = result.exit_code
        self.assertEqual(0, exit_code)
        self.assertFalse(os.path.exists("./res.tsv"), "Result file mustn't exists")

    def test_compare_semsim_inexistent_input(self):
        """test_compare_semsim"""
        params = [
            "-l",
            "l.tsv",
            "-r",
            "./tests/fixtures/right.tsv",
            "-o",
            "res.tsv",
            "-j",
            "GENE_SYMBOL",
            "-c",
            "MOUSE_PHENO_SCORE",
            "-c",
            "MOUSE_PHENO_SCOR",
        ]
        result = self.runner.invoke(
            compare_semsim,
            params,
        )
        err_log = self._caplog.records[0].message
        self.assertTrue("No such file or directory" in err_log)
        logging.info("ERR=%s", err_log)
        exit_code = result.exit_code
        self.assertEqual(0, exit_code)
        self.assertFalse(os.path.exists("./res.tsv"))
