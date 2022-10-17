import click
import logging
from .db import *


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


@main.command()
@click.option("-T", "--table", help="Number of greetings.")
@click.option("-S", "--scramble_factor", default=0.5, help="Number of greetings.")
def doscramble(table: str, scramble_factor: float):
    scramble_table(table, scramble_factor)


if __name__ == "__main__":
    main()
