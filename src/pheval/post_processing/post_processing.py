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

    def as_dict(self):
        return {'gene_symbol': self.pheval_gene_result.gene_symbol,
                'gene_identifier': self.pheval_gene_result.gene_identifier,
                'score': self.pheval_gene_result.score,
                'rank': self.rank}


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
