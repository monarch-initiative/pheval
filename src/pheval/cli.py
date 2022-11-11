import logging

import click

from .cli_pheval import run
from .cli_pheval_utils import scramble_phenopacket, scramble_semsim
from .containers import Container

info_log = logging.getLogger("info")


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet")
def main(verbose: int, quiet: bool) -> None:
    """main CLI method for pheval

    Args:
        verbose (int): _description_
        quiet (bool): _description_
    """
    if verbose >= 2:
        info_log.setLevel(level=logging.DEBUG)
    elif verbose == 1:
        info_log.setLevel(level=logging.INFO)
    else:
        info_log.setLevel(level=logging.WARNING)
    if quiet:
        info_log.setLevel(level=logging.ERROR)


container = Container()
container.wire(packages=["pheval"])


@click.group()
@click.option(
    "--runner",
    "-r",
    metavar="RUNNER",
    type=click.Choice(["exomiser", "foo"]),
    required=True,
    help="Runner",
)
def pheval(runner):
    container.config.runner.type.from_value(runner)


@click.group()
def pheval_utils():
    pass


pheval.add_command(run)

pheval_utils.add_command(scramble_semsim)
pheval_utils.add_command(scramble_phenopacket)
