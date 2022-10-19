import click
from .db import *
import subprocess
import os
import logging

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


def run_bash():
    return subprocess.check_call(
        f"{os.path.dirname(__file__)}/run.sh",
        stdout=sys.stdout,
        stderr=subprocess.STDOUT,
    )


@main.command()
@click.option("-T", "--table", help="Table Name", required=True)
@click.option("-S", "--scramble_factor", default=0.5, help="Scramble Factor")
def run(table: str, scramble_factor: float):
    run_bash()
    scramble_table(table, scramble_factor)
    rename_table(f"{table}_scramble", table)
    run_bash()
    new_table_name = f"{table}_scramble"
    rename_table(table, new_table_name)
    info_log.log("Done")


if __name__ == "__main__":
    main()
