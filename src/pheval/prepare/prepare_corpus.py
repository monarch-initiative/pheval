import shutil
import time
from pathlib import Path

from pheval.prepare.create_spiked_vcf import create_spiked_vcf
from pheval.prepare.update_phenopacket import create_updated_phenopacket
from pheval.utils.file_utils import all_files
from pheval.utils.logger import get_logger
from pheval.utils.phenopacket_utils import (
    PhenopacketUtil,
    create_gene_identifier_map,
    phenopacket_reader,
)

logger = get_logger()


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
    for phenopacket_path in all_files(phenopacket_dir):
        phenopacket_util = PhenopacketUtil(phenopacket_reader(phenopacket_path))
        if not phenopacket_util.observed_phenotypic_features():
            logger.warning(f"Removed {phenopacket_path.name} from the corpus due to no observed phenotypic features.")
            continue
        if variant_analysis:
            if phenopacket_util.check_incomplete_variant_record():
                logger.warning(f"Removed {phenopacket_path.name} from the corpus due to missing variant fields.")
                continue
            elif phenopacket_util.check_variant_alleles():
                logger.warning(
                    f"Removed {phenopacket_path.name} from the corpus due to identical "
                    "reference and alternate allele fields."
                )
        if gene_analysis:
            if phenopacket_util.check_incomplete_gene_record():
                logger.warning(f"Removed {phenopacket_path.name} from the corpus due to missing gene fields.")
                continue
        if disease_analysis:
            if phenopacket_util.check_incomplete_disease_record():
                logger.warning(f"Removed {phenopacket_path.name} from the corpus due to missing disease fields.")
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
        if gene_identifier:
            logger.info(f"Updating gene identifiers to {gene_identifier} for {phenopacket_dir}")
            create_updated_phenopacket(
                gene_identifier,
                phenopacket_path,
                output_dir.joinpath("phenopackets"),
                identifier_map,
            )
        else:
            # if not updating phenopacket gene identifiers then copy phenopacket as is to output directory
            (
                shutil.copy(phenopacket_path, output_dir.joinpath(f"phenopackets/{phenopacket_path.name}"))
                if phenopacket_path != output_dir.joinpath(f"phenopackets/{phenopacket_path.name}")
                else None
            )
    logger.info(
        f"Finished preparing corpus for {phenopacket_dir}. Total time: {time.perf_counter() - start_time:.2f} seconds."
    )
