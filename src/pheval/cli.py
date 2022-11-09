import click
import logging
from .cli_pheval import run
from .cli_pheval_utils import scramble_semsim, scramble_phenopacket


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

@click.group()
def pheval():
    pass

@click.group()
def pheval_utils():
    pass

pheval.add_command(run)

pheval_utils.add_command(scramble_semsim)
pheval_utils.add_command(scramble_phenopacket)

if __name__ == "__main__":
    main()
