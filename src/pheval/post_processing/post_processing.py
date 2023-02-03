from dataclasses import dataclass


@dataclass
class PhEvalGeneResult:
    gene_symbol: str
    gene_identifier: str
    score: float


@dataclass
class RankedPhEvalGeneResult:
    pheval_gene_result: PhEvalGeneResult
    rank: int


@dataclass
class PhEvalVariantResult:
    chromosome: str
    start: int
    stop: int
    ref: str
    alt: str


@dataclass
class RankedPhEvalVariantResult:
    pheval_variant_result: PhEvalVariantResult
    rank: int
