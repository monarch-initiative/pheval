"""PhEval CLI Module """
import logging

import click

from .cli_pheval import run
from .cli_pheval_utils import compare_semsim
from .cli_pheval_utils import scramble_phenopacket as sp
from .cli_pheval_utils import scramble_semsim

info_log = logging.getLogger("info")


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet")
def main(verbose=1, quiet=False) -> None:
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


@click.group()
def pheval():
    """pheval"""


@click.group()
def pheval_utils():
    """pheval_utils"""


pheval.add_command(run)

pheval_utils.add_command(scramble_semsim)
pheval_utils.add_command(sp)
pheval_utils.add_command(compare_semsim)

if __name__ == "__main__":
    main()
