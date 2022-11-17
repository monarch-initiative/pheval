"""
Monarch Initiative
"""

import click

from .pipeline import DefaultPhEvalRunner
from .runners.exomiser_pheval_runner import ExomiserPhEvalRunner


@click.command()
@click.option(
    "--inputdir",
    "-i",
    metavar="INPUTDIR",
    required=True,
    help="The input directory (relative path: e.g exomiser-13.11)",
)
@click.option(
    "--testdatadir", "-t", metavar="TESTDATA", required=True, help="Test data dir"
)
@click.option("--runner", "-r", metavar="RUNNER", required=True, help="Runner")
@click.option("--tmpdir", "-m", metavar="TMPDIR", required=False, help="Tmp dir")
@click.option("--outputdir", "-o", metavar="OUTPUTDIR", required=True, help="Ouput dir")
@click.option("--config", "-c", metavar="CONFIG", required=False, help="Config file")
def run(inputdir, testdatadir, runner, tmpdir, outputdir, config) -> None:
    """Run method
    Args:
        inputdir str: Input file
        testdatadir str: Test data dir
        runner str: Runner
        tmpdir str: Tmp dir
        outputdir str: Output dir
        config str: Config file
    """
    if runner == "exomiser":
        runner = ExomiserPhEvalRunner(inputdir, testdatadir, tmpdir, outputdir, config)
    else:
        runner = DefaultPhEvalRunner(inputdir, testdatadir, tmpdir, outputdir, config)

    runner.prepare()
    runner.run()
    runner.post_process()
