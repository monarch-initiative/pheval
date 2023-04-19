import operator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import pandas as pd


def calculate_end_pos(variant_start: int, variant_ref: str) -> int:
    """Calculate the end position for a variant."""
    return variant_start + len(variant_ref) - 1


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
        return {
            "gene_symbol": self.pheval_gene_result.gene_symbol,
            "gene_identifier": self.pheval_gene_result.gene_identifier,
            "score": self.pheval_gene_result.score,
            "rank": self.rank,
        }


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
        return {
            "chromosome": self.pheval_variant_result.chromosome,
            "start": self.pheval_variant_result.start,
            "end": self.pheval_variant_result.end,
            "ref": self.pheval_variant_result.ref,
            "alt": self.pheval_variant_result.alt,
            "score": self.pheval_variant_result.score,
            "rank": self.rank,
        }


class SortOrder(Enum):
    ASCENDING = 1
    DESCENDING = 2


class ResultSorter:
    def __init__(
        self, pheval_results: [PhEvalGeneResult] or [PhEvalVariantResult], sort_order: SortOrder
    ):
        self.pheval_results = pheval_results
        self.sort_order = sort_order

    def _sort_by_decreasing_score(self) -> [PhEvalGeneResult] or [PhEvalVariantResult]:
        """Sort results in descending order."""
        return sorted(self.pheval_results, key=operator.attrgetter("score"), reverse=True)

    def _sort_by_increasing_score(self) -> [PhEvalGeneResult] or [PhEvalVariantResult]:
        """Sort results in ascending order."""
        return sorted(self.pheval_results, key=operator.attrgetter("score"), reverse=False)

    def sort_pheval_results(self) -> [PhEvalGeneResult] or [PhEvalVariantResult]:
        """Sort results with best score first."""
        return (
            self._sort_by_increasing_score()
            if self.sort_order == SortOrder.ASCENDING
            else self._sort_by_decreasing_score()
        )


class ScoreRanker:
    rank: int = 0
    current_score: float = float("inf")
    count: int = 0

    def __init__(self, sort_order: SortOrder):
        self.sort_order = sort_order

    def _check_rank_order(self, round_score: float) -> None:
        """Check the results are correctly ordered."""
        if self.sort_order == SortOrder.ASCENDING and round_score < self.current_score != float(
            "inf"
        ):
            raise ValueError("Results are not correctly sorted!")
        elif self.sort_order == SortOrder.DESCENDING and round_score > self.current_score != float(
            "inf"
        ):
            raise ValueError("Results are not correctly sorted!")

    def rank_scores(self, round_score: float) -> int:
        """Add ranks to a result, equal scores are given the same rank e.g., 1,1,3."""
        self._check_rank_order(round_score)
        self.count += 1
        if self.current_score == round_score:
            return self.rank
        self.current_score = round_score
        self.rank = self.count
        return self.rank


def _rank_pheval_result(
    pheval_result: [PhEvalGeneResult] or [PhEvalVariantResult], sort_order: SortOrder
) -> [RankedPhEvalGeneResult] or [RankedPhEvalVariantResult]:
    """Ranks either a PhEval gene or variant result post-processed from a tool specific output.
    Deals with ex aequo scores"""
    score_ranker = ScoreRanker(sort_order)
    ranked_result = []
    for result in pheval_result:
        ranked_result.append(
            RankedPhEvalGeneResult(
                pheval_gene_result=result, rank=score_ranker.rank_scores(result.score)
            )
        ) if type(result) == PhEvalGeneResult else ranked_result.append(
            RankedPhEvalVariantResult(
                pheval_variant_result=result, rank=score_ranker.rank_scores(result.score)
            )
        )
    return ranked_result


def _return_sort_order(sort_order_str: str) -> SortOrder:
    """Return the SortOrder Enum from string derived from config."""
    try:
        return SortOrder[sort_order_str.upper()]
    except KeyError:
        raise ValueError("Incompatible ordering method specified.")


def _create_pheval_result(
    pheval_result: [PhEvalGeneResult] or [PhEvalVariantResult], sort_order_str: str
) -> [RankedPhEvalGeneResult] or [RankedPhEvalVariantResult]:
    """Create PhEval gene/variant result with corresponding ranks."""
    sort_order = _return_sort_order(sort_order_str)
    sorted_pheval_result = ResultSorter(pheval_result, sort_order).sort_pheval_results()
    return _rank_pheval_result(sorted_pheval_result, sort_order)


def _write_pheval_gene_result(
    ranked_pheval_result: [RankedPhEvalGeneResult], output_dir: Path, tool_result_path: Path
) -> None:
    """Write ranked PhEval gene result to tsv."""
    ranked_result = pd.DataFrame([x.as_dict() for x in ranked_pheval_result])
    pheval_gene_output = ranked_result.loc[:, ["rank", "score", "gene_symbol", "gene_identifier"]]
    pheval_gene_output.to_csv(
        output_dir.joinpath(
            "pheval_gene_results/" + tool_result_path.stem + "-pheval_gene_result.tsv"
        ),
        sep="\t",
        index=False,
    )


def _write_pheval_variant_result(
    ranked_pheval_result: [RankedPhEvalVariantResult], output_dir: Path, tool_result_path: Path
) -> None:
    """Write ranked PhEval variant result to tsv."""
    ranked_result = pd.DataFrame([x.as_dict() for x in ranked_pheval_result])
    pheval_variant_output = ranked_result.loc[
        :, ["rank", "score", "chromosome", "start", "end", "ref", "alt"]
    ]
    pheval_variant_output.to_csv(
        output_dir.joinpath(
            "pheval_variant_results/" + tool_result_path.stem + "-pheval_variant_result.tsv"
        ),
        sep="\t",
        index=False,
    )


def generate_pheval_result(
    pheval_result: [PhEvalGeneResult] or [PhEvalVariantResult],
    sort_order_str: str,
    output_dir: Path,
    tool_result_path: Path,
):
    """Generate either a PhEval variant or PhEval gene tsv result."""
    ranked_pheval_result = _create_pheval_result(pheval_result, sort_order_str)
    _write_pheval_variant_result(ranked_pheval_result, output_dir, tool_result_path) if all(
        isinstance(result, RankedPhEvalVariantResult) for result in ranked_pheval_result
    ) else _write_pheval_gene_result(ranked_pheval_result, output_dir, tool_result_path)
