from collections import defaultdict
from pathlib import Path

import click

from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import (
    PhenopacketRebuilder,
    PhenopacketUpdater,
    create_hgnc_dict,
    phenopacket_reader,
)


def update_outdated_gene_context(phenopacket: Path, gene_identifier: str, hgnc_data: defaultdict):
    phenopacket_contents = phenopacket_reader(phenopacket)
    updated_phenopacket_contents = PhenopacketUpdater(
        phenopacket, phenopacket_contents, hgnc_data, gene_identifier
    ).update_phenopacket_interpretations()
    PhenopacketRebuilder(updated_phenopacket_contents).write_phenopacket(phenopacket)


@click.command()
@click.option(
    "--phenopacket",
    "-p",
    metavar="FILE",
    required=True,
    help="Path to phenopacket for updating.",
    type=Path,
)
@click.option(
    "--gene-identifier",
    "-g",
    required=False,
    default="ensembl_id",
    show_default=True,
    help="Gene identifier to add to phenopacket",
    type=click.Choice(["ensembl_id", "entrez_id", "hgnc_id"]),
)
def update_phenopacket(phenopacket: Path, gene_identifier: str):
    """Update gene symbols and identifiers for a phenopacket."""
    hgnc_data = create_hgnc_dict()
    update_outdated_gene_context(phenopacket, gene_identifier, hgnc_data)


@click.command()
@click.option(
    "--phenopacket-dir",
    "-p",
    metavar="PATH",
    required=True,
    help="Path to phenopacket directory for updating.",
    type=Path,
)
@click.option(
    "--gene-identifier",
    "-g",
    required=False,
    default="ensembl_id",
    show_default=True,
    help="Gene identifier to add to phenopacket",
    type=click.Choice(["ensembl_id", "entrez_id", "hgnc_id"]),
)
def update_phenopackets(phenopacket_dir: Path, gene_identifier: str):
    """Update gene symbols and identifiers for phenopackets."""
    hgnc_data = create_hgnc_dict()
    for phenopacket in all_files(phenopacket_dir):
        update_outdated_gene_context(phenopacket, gene_identifier, hgnc_data)
