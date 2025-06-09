"""
Monarch Initiative
"""

import time
from pathlib import Path

import click

from pheval.implementations import get_implementation_resolver
from pheval.utils.file_utils import write_metadata
from pheval.utils.logger import get_logger
from pheval.utils.utils import download_hgnc_data, download_mondo_mapping

logger = get_logger()


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
    logger.info(f"Executing {runner}.")
    start_time = time.perf_counter()
    runner_class = get_implementation_resolver().lookup(runner)
    runner_instance = runner_class(input_dir, testdata_dir, tmp_dir, output_dir, config, version)
    runner_instance.build_output_directory_structure()
    logger.info("Executing prepare phase.")
    runner_instance.prepare()
    logger.info("Executing run phase.")
    runner_instance.run()
    logger.info("Executing post-processing phase.")
    runner_instance.post_process()
    run_metadata = runner_instance.construct_meta_data()
    logger.info(f"Writing metadata for run to {output_dir}.")
    write_metadata(output_dir, run_metadata)
    logger.info(f"Run completed! Total time: {time.perf_counter() - start_time:.2f} seconds.")


@click.command()
def update():
    """
    Download the latest MONDO and HGNC mapping files.

    This command fetches the most recent versions of:

    * The MONDO SSSOM mapping file from the Monarch Initiative

    * The HGNC complete gene set from the HGNC download site

    These files are saved to the `resources/` directory and will
    overwrite any existing versions. This ensures that PhEval has
    access to the most up-to-date identifier mappings for disease
    and gene normalisation.
    """
    download_mondo_mapping()
    download_hgnc_data()
