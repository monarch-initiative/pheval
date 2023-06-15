from collections import defaultdict
from pathlib import Path

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
        hgnc_data=hgnc_data, gene_identifier=gene_identifier
    ).update_genomic_interpretations_gene_identifier(interpretations)

    return PhenopacketRebuilder(phenopacket).update_interpretations(updated_interpretations)


def create_updated_phenopacket(gene_identifier: str, phenopacket_path: Path, output_dir: Path):
    """Updates the gene context within the interpretations for a phenopacket."""
    hgnc_data = create_hgnc_dict()
    updated_phenopacket = update_outdated_gene_context(phenopacket_path, gene_identifier, hgnc_data)
    write_phenopacket(updated_phenopacket, output_dir.joinpath(phenopacket_path.name))


def create_updated_phenopackets(gene_identifier: str, phenopacket_dir: Path, output_dir: Path):
    """Updates the gene context within the interpretations for phenopackets."""
    hgnc_data = create_hgnc_dict()
    for phenopacket_path in all_files(phenopacket_dir):
        updated_phenopacket = update_outdated_gene_context(
            phenopacket_path, gene_identifier, hgnc_data
        )
        write_phenopacket(updated_phenopacket, output_dir.joinpath(phenopacket_path.name))


def update_phenopackets(
    gene_identifier: str, phenopacket_path: Path, phenopacket_dir: Path, output_dir: Path
):
    """Update the gene identifiers in either a single phenopacket or a directory of phenopackets."""
    output_dir.mkdir(exist_ok=True)
    if phenopacket_path is not None:
        create_updated_phenopacket(gene_identifier, phenopacket_path, output_dir)
    elif phenopacket_dir is not None:
        create_updated_phenopackets(gene_identifier, phenopacket_dir, output_dir)
