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
        input_arg = "x"
        result = self.runner.invoke(run, ["-i", input_arg])
        err = result.stderr
        logging.info(f"ERR={err}")
        self.assertEqual(1, result.exit_code)
