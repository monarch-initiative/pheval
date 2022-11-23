"""
Monarch Initiative
"""

import click

from pheval.implementations import get_implementation_resolver


@click.command()
@click.option(
    "--inputdir",
    "-i",
    metavar="INPUTDIR",
    required=True,
    help="The input directory (relative path: e.g exomiser-13.11)",
)
@click.option(
    "--testdatadir",
    "-t",
    metavar="TESTDATA",
    required=True,
    help="The input directory (relative path: e.g ./data)",
)
@click.option(
    "--runner",
    "-r",
    metavar="RUNNER",
    required=True,
    help="Runner implementation (e.g exomiser-13.11)",
)
@click.option(
    "--tmpdir",
    "-m",
    metavar="TMPDIR",
    required=False,
    help="The path of the temporary directory (optional)",
)
@click.option(
    "--outputdir",
    "-o",
    metavar="OUTPUTDIR",
    required=True,
    help="The path of the output directory",
)
@click.option(
    "--config",
    "-c",
    metavar="CONFIG",
    required=False,
    help="The path of the configuration file (optional e.g config.yaml)",
)
def run(inputdir, testdatadir, runner, tmpdir, outputdir, config) -> None:
    """PhEval Runner Command Line Interface

    Args:
        inputdir (Click.Path): The input directory (relative path: e.g exomiser-13.11)
        testdatadir (Click.Path): The input directory (relative path: e.g ./data
        runner (str): Runner implementation (e.g exomiser-13.11)
        tmpdir (Click.Path): The path of the temporary directory (optional)
        outputdir (Click.Path): The path of the output directory
        config (Click.Path): The path of the configuration file (optional e.g config.yaml)
    """
    runner_class = get_implementation_resolver().lookup(runner)
    runner_instance = runner_class(inputdir, testdatadir, tmpdir, outputdir, config)
    runner_instance.prepare()
    runner_instance.run()
    runner_instance.post_process()
