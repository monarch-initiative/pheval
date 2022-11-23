"""CLI Test """
import logging
import unittest

from click.testing import CliRunner

from pheval.cli_pheval import run


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
