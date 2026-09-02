import shutil
import time
from pathlib import Path

import polars as pl

from pheval.prepare.create_spiked_vcf import create_spiked_vcf
from pheval.prepare.update_phenopacket import create_updated_phenopacket
from pheval.utils.file_utils import files_with_suffix
from pheval.utils.logger import get_logger
from pheval.utils.phenopacket_utils import (
    PhenopacketUtil,
    create_gene_identifier_map,
    phenopacket_reader,
)

logger = get_logger()


def _exclude_phenopacket(
    phenopacket_util: PhenopacketUtil,
    phenopacket_name: str,
    variant_analysis: bool,
    gene_analysis: bool,
    disease_analysis: bool,
) -> bool:
    """
    Check whether a Phenopacket fails a completeness check and should be dropped from the corpus.

    Logs the reason for exclusion as a side effect.

    Args:
        phenopacket_util (PhenopacketUtil): Wrapper around the Phenopacket being checked.
        phenopacket_name (str): Name of the Phenopacket file, used in log messages.
        variant_analysis (bool): If True, check for complete variant records.
        gene_analysis (bool): If True, check for complete gene records.
        disease_analysis (bool): If True, check for complete disease records.

    Returns:
        bool: True if the Phenopacket should be excluded from the corpus.
    """
    if not phenopacket_util.observed_phenotypic_features():
        logger.warning(f"Removed {phenopacket_name} from the corpus due to no observed phenotypic features.")
        return True
    if variant_analysis:
        if phenopacket_util.check_incomplete_variant_record():
            logger.warning(f"Removed {phenopacket_name} from the corpus due to missing variant fields.")
            return True
        elif phenopacket_util.check_variant_alleles():
            logger.warning(
                f"Removed {phenopacket_name} from the corpus due to identical reference and alternate allele fields."
            )
            return True
    if gene_analysis and phenopacket_util.check_incomplete_gene_record():
        logger.warning(f"Removed {phenopacket_name} from the corpus due to missing gene fields.")
        return True
    if disease_analysis and phenopacket_util.check_incomplete_disease_record():
        logger.warning(f"Removed {phenopacket_name} from the corpus due to missing disease fields.")
        return True
    return False


def _write_phenopacket(
    phenopacket_path: Path,
    output_dir: Path,
    gene_identifier: str,
    identifier_map: pl.DataFrame,
) -> None:
    """
    Write a Phenopacket to the output directory, updating its gene identifiers first if requested.

    Args:
        phenopacket_path (Path): Path to the source Phenopacket file.
        output_dir (Path): The directory to save the prepared Phenopacket to.
        gene_identifier (str): Identifier to update gene identifiers to, if applicable.
        identifier_map (pl.DataFrame): Gene identifier lookup table.
    """
    output_phenopacket_path = output_dir.joinpath(f"phenopackets/{phenopacket_path.name}")
    if gene_identifier:
        logger.info(f"Updating gene identifiers to {gene_identifier} for {phenopacket_path.name}")
        create_updated_phenopacket(
            gene_identifier,
            phenopacket_path,
            output_dir.joinpath("phenopackets"),
            identifier_map,
        )
    # if not updating phenopacket gene identifiers then copy phenopacket as is to output directory
    elif phenopacket_path != output_phenopacket_path:
        shutil.copy(phenopacket_path, output_phenopacket_path)


def prepare_corpus(
    phenopacket_dir: Path,
    variant_analysis: bool,
    gene_analysis: bool,
    disease_analysis: bool,
    gene_identifier: str,
    hg19_template_vcf: Path,
    hg38_template_vcf: Path,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
    output_dir: Path,
) -> None:
    """
    Prepare a corpus of Phenopackets for analysis, optionally checking for complete variant records and updating
    gene identifiers.

    Args:
        phenopacket_dir (Path): The path to the directory containing Phenopackets.
        variant_analysis (bool): If True, check for complete variant records in the Phenopackets.
        gene_analysis (bool): If True, check for complete gene records in the Phenopackets.
        disease_analysis (bool): If True, check for complete disease records in the Phenopackets.
        gene_identifier (str): Identifier for updating gene identifiers, if applicable.
        hg19_template_vcf (Path): Path to the hg19 template VCF file (optional), to spike variants into
        VCFs for variant-based analysis at least one of hg19_template_vcf or hg38_template_vcf is required.
        hg38_template_vcf (Path): Path to the hg38 template VCF file (optional), to spike variants into
        VCFs for variant-based analysis at least one of hg19_template_vcf or hg38_template_vcf is required.
        hg19_vcf_dir (Path): Path to the directory containing hg19 template VCF files (optional).
        hg38_vcf_dir (Path): Path to the directory containing hg38 template VCF files (optional).
        output_dir (Path): The directory to save the prepared Phenopackets and, optionally, VCF files.
    Notes:
        To spike variants into VCFs for variant-based analysis at least one of hg19_template_vcf, hg38_template_vcf,
        hg19_vcf_dir or hg38_vcf_dir is required.
    """
    start_time = time.perf_counter()
    logger.info(f"Preparing corpus for {phenopacket_dir}")
    output_dir.joinpath("phenopackets").mkdir(exist_ok=True, parents=True)
    logger.info(f" Created output directory: {output_dir.joinpath('phenopackets')}")
    identifier_map = create_gene_identifier_map()
    for phenopacket_path in files_with_suffix(phenopacket_dir, ".json"):
        phenopacket_util = PhenopacketUtil(phenopacket_reader(phenopacket_path))
        if _exclude_phenopacket(
            phenopacket_util, phenopacket_path.name, variant_analysis, gene_analysis, disease_analysis
        ):
            continue
        logger.info(f"{phenopacket_path.name} OK!")
        if hg19_template_vcf or hg38_template_vcf:
            output_dir.joinpath("vcf").mkdir(exist_ok=True)
            logger.info(f" Created output directory: {output_dir.joinpath('vcf')}")
            logger.info(f"Spiking VCF for {phenopacket_path}.")
            create_spiked_vcf(
                output_dir.joinpath("vcf"),
                phenopacket_path,
                hg19_template_vcf,
                hg38_template_vcf,
                hg19_vcf_dir,
                hg38_vcf_dir,
            )
        _write_phenopacket(phenopacket_path, output_dir, gene_identifier, identifier_map)
    logger.info(
        f"Finished preparing corpus for {phenopacket_dir}. Total time: {time.perf_counter() - start_time:.2f} seconds."
    )
