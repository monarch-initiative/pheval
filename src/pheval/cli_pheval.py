"""
Monarch Initiative
"""
from pathlib import Path

import click

from pheval.implementations import get_implementation_resolver
from pheval.utils.file_utils import write_metadata


@click.command()
@click.option(
    "--input-dir",
    "-i",
    metavar="INPUTDIR",
    required=True,
    help="The input directory (relative path: e.g exomiser-13.11)",
    type=Path,
)
@click.option(
    "--testdata-dir",
    "-t",
    metavar="TESTDATA",
    required=True,
    help="The input directory (relative path: e.g ./data)",
    type=Path,
)
@click.option(
    "--runner",
    "-r",
    metavar="RUNNER",
    required=True,
    help="Runner implementation (e.g exomiser-13.11)",
)
@click.option(
    "--tmp-dir",
    "-m",
    metavar="TMPDIR",
    required=False,
    help="The path of the temporary directory (optional)",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    metavar="OUTPUTDIR",
    required=True,
    help="The path of the output directory",
    type=Path,
)
@click.option(
    "--config",
    "-c",
    metavar="CONFIG",
    required=False,
    help="The path of the configuration file (optional e.g config.yaml)",
    type=Path,
)
@click.option(
    "--version",
    "-v",
    required=False,
    help="Version of the tool implementation.",
    type=str,
)
def run(
    input_dir: Path,
    testdata_dir: Path,
    runner: str,
    tmp_dir: Path,
    output_dir: Path,
    config: Path,
    version: str,
) -> None:
    """PhEval Runner Command Line Interface
    Args:
        input_dir (Path): The input directory (relative path: e.g exomiser-13.11)
        testdata_dir (Path): The input directory (relative path: e.g ./data
        runner (str): Runner implementation (e.g exomiser-13.11)
        tmp_dir (Path): The path of the temporary directory (optional)
        output_dir (Path): The path of the output directory
        config (Path): The path of the configuration file (optional e.g., config.yaml)
        version (str): The version of the tool implementation
    """
    runner_class = get_implementation_resolver().lookup(runner)
    runner_instance = runner_class(input_dir, testdata_dir, tmp_dir, output_dir, config, version)
    runner_instance.build_output_directory_structure()
    runner_instance.prepare()
    runner_instance.run()
    runner_instance.post_process()
    run_metadata = runner_instance.construct_meta_data()
    write_metadata(output_dir, run_metadata)
