import logging
import unittest

from click.testing import CliRunner

from pheval.cli_pheval import run


class TestCommandLineInterface(unittest.TestCase):
    """
    Tests all command-line subcommands
    """

    # @click.option(
    #     "--inputdir",
    #     "-i",
    #     metavar="INPUTDIR",
    #     required=True,
    #     help="The input directory (relative path: e.g exomiser-13.11)",
    # )
    # @click.option(
    #     "--testdatadir", "-t", metavar="TESTDATA", required=True, help="Test data dir"
    # )
    # @click.option("--runner", "-r", metavar="RUNNER", required=True, help="Runner")
    # @click.option("--tmpdir", "-m", metavar="TMPDIR", required=True, help="Tmp dir")
    # @click.option("--outputdir", "-o", metavar="OUTPUTDIR", required=True, help="Ouput dir")
    # @click.option("--config", "-c", metavar="CONFIG", required=True, help="Config file")

    def setUp(self) -> None:
        runner = CliRunner(mix_stderr=False)
        self.runner = runner

    def test_scramble_semsim(self):
        result = self.runner.invoke(run, ["-i", "./", "-t", "./", "-r", "exomiser", "-o", "./"])
        err = result.stderr
        self.assertEqual("", err)
        logging.info(f"ERR={err}")
        exit_code = result.exit_code
        self.assertEqual(0, exit_code)
