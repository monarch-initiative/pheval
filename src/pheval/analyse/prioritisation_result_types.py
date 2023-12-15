from dataclasses import dataclass
from pathlib import Path

from pheval.utils.phenopacket_utils import GenomicVariant, ProbandDisease


@dataclass
class GenePrioritisationResult:
    """
    Store rank data for causative genes.

    Attributes:
        phenopacket_path (Path): Path to the phenopacket.
        gene (str): The causative gene.
        rank (int): The assigned rank for the gene. Defaults to 0.
    """

    phenopacket_path: Path
    gene: str
    rank: int = 0


@dataclass
class VariantPrioritisationResult:
    """
    Store rank data for variants.

    Attributes:
        phenopacket_path (Path): Path to the phenopacket.
        variant (GenomicVariant): The genomic variant.
        rank (int): The assigned rank for the variant. Defaults to 0.
    """

    phenopacket_path: Path
    variant: GenomicVariant
    rank: int = 0


@dataclass
class DiseasePrioritisationResult:
    """
    Store rank data for known diseases.

    Attributes:
        phenopacket_path (Path): Path to the phenopacket.
        disease (ProbandDisease): The proband disease.
        rank (int): The assigned rank for the disease. Defaults to 0.
    """

    phenopacket_path: Path
    disease: ProbandDisease
    rank: int = 0
