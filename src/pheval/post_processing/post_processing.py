from dataclasses import dataclass


@dataclass
class PhEvalGeneResult:
    """Minimal data required from tool-specific output for gene prioritisation."""
    gene_symbol: str
    gene_identifier: str
    score: float


@dataclass
class RankedPhEvalGeneResult:
    """PhEval gene result with corresponding rank."""
    pheval_gene_result: PhEvalGeneResult
    rank: int


@dataclass
class PhEvalVariantResult:
    """Minimal data required from tool-specific output for variant prioritisation."""
    chromosome: str
    start: int
    stop: int
    ref: str
    alt: str


@dataclass
class RankedPhEvalVariantResult:
    """PhEval variant result with corresponding rank."""
    pheval_variant_result: PhEvalVariantResult
    rank: int
