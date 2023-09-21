from dataclasses import dataclass
from pathlib import Path

from pheval.utils.phenopacket_utils import GenomicVariant, ProbandDisease


@dataclass
class GenePrioritisationResult:
    """Store rank data for causative genes."""

    phenopacket_path: Path
    gene: str
    rank: int = 0


@dataclass
class VariantPrioritisationResult:
    """Store rank data for causative variants."""

    phenopacket_path: Path
    variant: GenomicVariant
    rank: int = 0


@dataclass
class DiseasePrioritisationResult:
    """Store rank data for known diseases."""

    phenopacket_path: Path
    disease: ProbandDisease
    rank: int = 0
