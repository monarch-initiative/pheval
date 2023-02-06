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
        """Return PhEval gene result as dictionary."""
        return {'gene_symbol': self.pheval_gene_result.gene_symbol,
                'gene_identifier': self.pheval_gene_result.gene_identifier,
                'score': self.pheval_gene_result.score,
                'rank': self.rank}


@dataclass
class PhEvalVariantResult:
    """Minimal data required from tool-specific output for variant prioritisation."""
    chromosome: str
    start: int
    end: int
    ref: str
    alt: str
    score: float


@dataclass
class RankedPhEvalVariantResult:
    """PhEval variant result with corresponding rank."""
    pheval_variant_result: PhEvalVariantResult
    rank: int

    def as_dict(self):
        """Return PhEval variant result as dictionary."""
        return {'chromosome': self.pheval_variant_result.chromosome,
                'start': self.pheval_variant_result.start,
                'end': self.pheval_variant_result.end,
                'ref': self.pheval_variant_result.ref,
                'alt': self.pheval_variant_result.alt,
                'score': self.pheval_variant_result.score,
                'rank': self.rank}
