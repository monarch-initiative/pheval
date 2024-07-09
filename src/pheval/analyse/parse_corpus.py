from pathlib import Path
from typing import List

from pheval.utils.phenopacket_utils import GenomicVariant, ProbandCausativeGene, phenopacket_reader, PhenopacketUtil, \
    ProbandDisease


def _obtain_causative_diseases(phenopacket_path: Path) -> List[ProbandDisease]:
    """
    Obtain known diseases from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.

    Returns:
       List[ProbandDisease]: A list of known diseases associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnoses()


def _obtain_causative_variants(phenopacket_path: Path) -> List[GenomicVariant]:
    """
    Obtain known variants from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.

    Returns:
       List[GenomicVariant]: A list of known variants associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_variants()


def _obtain_causative_genes(phenopacket_path: Path) -> List[ProbandCausativeGene]:
    """
    Obtain known genes from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.
    Returns:
       List[ProbandCausativeGene]: A list of known genes associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_genes()
