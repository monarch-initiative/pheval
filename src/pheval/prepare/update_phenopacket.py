from collections import defaultdict
from pathlib import Path

import click

from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import (
    GeneIdentifierUpdater,
    PhenopacketRebuilder,
    PhenopacketUtil,
    create_hgnc_dict,
    phenopacket_reader,
    write_phenopacket,
)


def update_outdated_gene_context(
    phenopacket_path: Path, gene_identifier: str, hgnc_data: defaultdict
):
    """Updates the gene context of the phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    interpretations = PhenopacketUtil(phenopacket).interpretations()
    updated_interpretations = GeneIdentifierUpdater(
        hgnc_data, gene_identifier
    ).update_genomic_interpretations_gene_identifier(interpretations)

    return PhenopacketRebuilder(phenopacket).update_interpretations(updated_interpretations)


def update_phenopacket(gene_identifier: str, phenopacket_path: Path, output_dir: Path):
    """Updates the gene context within the interpretations for a phenopacket."""
    hgnc_data = create_hgnc_dict()
    updated_phenopacket = update_outdated_gene_context(phenopacket_path, gene_identifier, hgnc_data)
    write_phenopacket(updated_phenopacket, output_dir.joinpath(phenopacket_path.name))


def update_phenopackets(gene_identifier: str, phenopacket_dir: Path, output_dir: Path):
    """Updates the gene context within the interpretations for phenopackets."""
    hgnc_data = create_hgnc_dict()
    for phenopacket_path in all_files(phenopacket_dir):
        updated_phenopacket = update_outdated_gene_context(
            phenopacket_path, gene_identifier, hgnc_data
        )
        write_phenopacket(updated_phenopacket, output_dir.joinpath(phenopacket_path.name))


@click.command("update-phenopacket")
@click.option(
    "--phenopacket-path",
    "-p",
    metavar="FILE",
    required=True,
    help="Path to phenopacket for updating.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    metavar="PATH",
    required=True,
    help="Path to write phenopacket.",
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
def update_phenopacket_command(phenopacket_path: Path, output_dir: Path, gene_identifier: str):
    """Update gene symbols and identifiers for a phenopacket."""
    update_phenopacket(gene_identifier, phenopacket_path, output_dir)


@click.command("update-phenopackets")
@click.option(
    "--phenopacket-dir",
    "-p",
    metavar="PATH",
    required=True,
    help="Path to phenopacket directory for updating.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    metavar="PATH",
    required=True,
    help="Path to write phenopacket.",
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
def update_phenopackets_command(phenopacket_dir: Path, output_dir: Path, gene_identifier: str):
    """Update gene symbols and identifiers for phenopackets."""
    update_phenopackets(gene_identifier, phenopacket_dir, output_dir)
