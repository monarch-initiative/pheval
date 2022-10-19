import click
import logging
from .db import *
import subprocess
import os


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet")
def main(verbose: int, quiet: bool) -> None:
    """main CLI method for pheval

    Args:
        verbose (int): _description_
        quiet (bool): _description_
    """
    logger = logging.getLogger()
    if verbose >= 2:
        logger.setLevel(level=logging.DEBUG)
    elif verbose == 1:
        logger.setLevel(level=logging.INFO)
    else:
        logger.setLevel(level=logging.WARNING)
    if quiet:
        logger.setLevel(level=logging.ERROR)


# def doscramble(table: str, scramble_factor: float):
#     scramble_table(table, scramble_factor)


@main.command()
@click.option("-T", "--table", help="Number of greetings.")
@click.option("-S", "--scramble_factor", default=0.5, help="Number of greetings.")
def run(table: str, scramble_factor: float):
    res = subprocess.check_call(
        f"{os.path.dirname(__file__)}/run.sh",
        stdout=sys.stdout,
        stderr=subprocess.STDOUT,
    )
    scramble_table(table, scramble_factor)


if __name__ == "__main__":
    main()
