import click
import pheval.db as db
import sys
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


def run_bash(table=None):
    print(f"{os.path.dirname(__file__)}/run.sh")
    if table is None:
        return subprocess.check_call(
            [f"{os.path.dirname(__file__)}/run.sh"],
            stdout=sys.stdout,
            stderr=subprocess.STDOUT,
        )
    return subprocess.check_call(
        [f"{os.path.dirname(__file__)}/run.sh", "-t", table],
        stdout=sys.stdout,
        stderr=subprocess.STDOUT,
    )


def dump(table):
    return subprocess.check_call(
        [f"{os.path.dirname(__file__)}/run.sh", "-d", table],
        stdout=sys.stdout,
        stderr=subprocess.STDOUT,
    )


valid_tables = ["HP_HP_MAPPINGS", "HP_MP_MAPPINGS", "HP_ZP_MAPPINGS"]


def validate_table(table):
    valid_tables.append("ALL")
    if table.upper() not in valid_tables:
        raise Exception(f"{table} is invalid")


def iterate_table(table, scramble_factor):
    if table.upper() != "ALL":
        return pipeline(table, scramble_factor, first=True)
    for idx, t in enumerate(valid_tables):
        print(t)
        pipeline(t, scramble_factor, first=idx == 0)


@main.command()
@click.option("-T", "--table", help="Table Name", required=True)
@click.option("-S", "--scramble_factor", default=0.5, help="Scramble Factor")
def run(table: str, scramble_factor: float):
    validate_table(table)
    iterate_table(table, scramble_factor)
    logging.info("Done")


def pipeline(table, scramble_factor, first=True):
    if first:
        run_bash()
    info_log.info("P1")
    db.clean_aux_table(table)
    dump(table)
    db.rename_table(f"{table}", f"{table}_ORIGINAL")
    info_log.info("P2")
    db.scramble_table(table, scramble_factor)
    info_log.info("P3")
    db.rename_table(f"{table}_SCRAMBLE", table)
    run_bash(table)
    info_log.info("P4")
    db.rename_table(table, f"{table}_SCRAMBLE")
    info_log.info("P5")
    db.rename_table(f"{table}_ORIGINAL", f"{table}")


if __name__ == "__main__":
    main()
