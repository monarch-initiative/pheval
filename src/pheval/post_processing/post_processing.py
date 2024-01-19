import logging
import operator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import pandas as pd

info_log = logging.getLogger("info")


def calculate_end_pos(variant_start: int, variant_ref: str) -> int:
    """Calculate the end position for a variant
    Args:
        variant_start (int): The start position of the variant
        variant_ref (str): The reference allele of the variant

    Returns:
        int: The end position of the variant
    """
    return variant_start + len(variant_ref) - 1


@dataclass
class PhEvalResult:
    """Base class for PhEval results."""


@dataclass
class PhEvalGeneResult(PhEvalResult):
    """Minimal data required from tool-specific output for gene prioritisation result
    Args:
        gene_symbol (str): The gene symbol for the result entry
        gene_identifier (str): The ENSEMBL gene identifier for the result entry
        score (float): The score for the gene result entry
    Notes:
        While we recommend providing the gene identifier in the ENSEMBL namespace,
        any matching format used in Phenopacket interpretations is acceptable for result matching purposes
        in the analysis.
    """

    gene_symbol: str
    gene_identifier: str
    score: float


@dataclass
class RankedPhEvalGeneResult(PhEvalGeneResult):
    """PhEval gene result with corresponding rank
    Args:
        rank (int): The rank for the result entry
    """

    rank: int

    @staticmethod
    def from_gene_result(pheval_gene_result: PhEvalGeneResult, rank: int):
        """Return RankedPhEvalGeneResult from a PhEvalGeneResult and rank
        Args:
            pheval_gene_result (PhEvalGeneResult): The gene result entry
            rank (int): The corresponding rank for the result entry

        Returns:
            RankedPhEvalGeneResult: The result as a RankedPhEvalGeneResult
        """
        return RankedPhEvalGeneResult(
            gene_symbol=pheval_gene_result.gene_symbol,
            gene_identifier=pheval_gene_result.gene_identifier,
            score=pheval_gene_result.score,
            rank=rank,
        )


@dataclass
class PhEvalVariantResult(PhEvalResult):
    """Minimal data required from tool-specific output for variant prioritisation
    Args:
        chromosome (str): The chromosome position of the variant recommended to be provided in the following format.
        This includes numerical designations from 1 to 22 representing autosomal chromosomes,
        as well as the sex chromosomes X and Y, and the mitochondrial chromosome MT.
        start (int): The start position of the variant
        end (int): The end position of the variant
        ref (str): The reference allele of the variant
        alt (str): The alternate allele of the variant
        score (float): The score for the variant result entry
    Notes:
        While we recommend providing the variant's chromosome in the specified format,
        any matching format used in Phenopacket interpretations is acceptable for result matching purposes
        in the analysis.
    """

    chromosome: str
    start: int
    end: int
    ref: str
    alt: str
    score: float


@dataclass
class RankedPhEvalVariantResult(PhEvalVariantResult):
    """PhEval variant result with corresponding rank
    Args:
        rank (int): The rank for the result entry
    """

    rank: int

    @staticmethod
    def from_variant_result(pheval_variant_result: PhEvalVariantResult, rank: int):
        """Return RankedPhEvalVariantResult from a PhEvalVariantResult and rank
        Args:
            pheval_variant_result (PhEvalVariantResult): The variant result entry
            rank (int): The corresponding rank for the result entry

        Returns:
            RankedPhEvalVariantResult: The result as a RankedPhEvalVariantResult
        """
        return RankedPhEvalVariantResult(
            chromosome=pheval_variant_result.chromosome,
            start=pheval_variant_result.start,
            end=pheval_variant_result.end,
            ref=pheval_variant_result.ref,
            alt=pheval_variant_result.alt,
            score=pheval_variant_result.score,
            rank=rank,
        )


@dataclass
class PhEvalDiseaseResult(PhEvalResult):
    """Minimal data required from tool-specific output for disease prioritisation
    Args:
        disease_name (str): Disease name for the result entry
        disease_identifier (str): Identifier for the disease result entry in the OMIM namespace
        score (str): Score for the disease result entry
    Notes:
        While we recommend providing the disease identifier in the OMIM namespace,
        any matching format used in Phenopacket interpretations is acceptable for result matching purposes
        in the analysis.
    """

    disease_name: str
    disease_identifier: str
    score: float


@dataclass
class RankedPhEvalDiseaseResult(PhEvalDiseaseResult):
    """PhEval disease result with corresponding rank
    Args:
        rank (int): The rank for the result entry
    """

    rank: int

    @staticmethod
    def from_disease_result(pheval_disease_result: PhEvalDiseaseResult, rank: int):
        """Return RankedPhEvalDiseaseResult from a PhEvalDiseaseResult and rank
        Args:
            pheval_disease_result (PhEvalDiseaseResult): The disease result entry
            rank (int): The corresponding rank for the result entry

        Returns:
            RankedPhEvalDiseaseResult: The result as a RankedPhEvalDiseaseResult
        """
        return RankedPhEvalDiseaseResult(
            disease_name=pheval_disease_result.disease_name,
            disease_identifier=pheval_disease_result.disease_identifier,
            score=pheval_disease_result.score,
            rank=rank,
        )


class SortOrder(Enum):
    """Enumeration representing sorting orders."""

    ASCENDING = 1
    """Ascending sort order."""
    DESCENDING = 2
    """Descending sort order."""


class ResultSorter:
    """Class for sorting PhEvalResult instances based on a given sort order."""

    def __init__(self, pheval_results: [PhEvalResult], sort_order: SortOrder):
        """
        Initialise ResultSorter

        Args:
            pheval_results ([PhEvalResult]): List of PhEvalResult instances to be sorted
            sort_order (SortOrder): Sorting order to be applied
        """
        self.pheval_results = pheval_results
        self.sort_order = sort_order

    def _sort_by_decreasing_score(self) -> [PhEvalResult]:
        """
        Sort results in descending order based on the score

        Returns:
            [PhEvalResult]: Sorted list of PhEvalResult instances.
        """
        return sorted(self.pheval_results, key=operator.attrgetter("score"), reverse=True)

    def _sort_by_increasing_score(self) -> [PhEvalResult]:
        """
        Sort results in ascending order based on the score

        Returns:
            [PhEvalResult]: Sorted list of PhEvalResult instances.
        """
        return sorted(self.pheval_results, key=operator.attrgetter("score"), reverse=False)

    def sort_pheval_results(self) -> [PhEvalResult]:
        """
        Sort results based on the specified sort order.

        Returns:
            [PhEvalResult]: Sorted list of PhEvalResult instances.
        """
        return (
            self._sort_by_increasing_score()
            if self.sort_order == SortOrder.ASCENDING
            else self._sort_by_decreasing_score()
        )


class ScoreRanker:
    """
    Class for ranking scores based on a given sort order

    Attributes:
       rank (int): Represents the current rank, initialised with 0
       current_score (float): Represents the current score, initialised with positive infinity (float("inf"))
       count (int): Used for counting, initialised with 0
    """

    rank: int = 0
    current_score: float = float("inf")
    count: int = 0

    def __init__(self, sort_order: SortOrder):
        """
        Initialise ScoreRanker

        Args:
            sort_order (SortOrder): Sorting order to be applied
        """
        self.sort_order = sort_order

    def _check_rank_order(self, round_score: float) -> None:
        """
        Check if the results are correctly ordered

        Args:
            round_score (float): Score to be checked against the current score

        Raises:
            ValueError: If results are not correctly sorted.
        """
        if self.sort_order == SortOrder.ASCENDING and round_score < self.current_score != float(
            "inf"
        ):
            raise ValueError("Results are not correctly sorted!")
        elif self.sort_order == SortOrder.DESCENDING and round_score > self.current_score != float(
            "inf"
        ):
            raise ValueError("Results are not correctly sorted!")

    def rank_scores(self, round_score: float) -> int:
        """
        Add ranks to a result; equal scores are given the same rank, e.g., 1, 1, 3

        Args:
            round_score (float): Score to be ranked

        Returns:
            int: Rank assigned to the score
        """
        self._check_rank_order(round_score)
        self.count += 1
        if self.current_score == round_score:
            return self.rank
        self.current_score = round_score
        self.rank = self.count
        return self.rank


def _rank_pheval_result(pheval_result: [PhEvalResult], sort_order: SortOrder) -> [PhEvalResult]:
    """
    Rank PhEval results post-processed from tool-specific output, managing tied scores (ex aequo)

    Args:
        pheval_result ([PhEvalResult]): PhEval results obtained from tool-specific output
        sort_order (SortOrder): Sorting order based on which ranking is performed

    Returns:
        List[PhEvalResult]: Ranked PhEval results with tied scores managed

    Raises:
        ValueError: If an incompatible PhEval result type is encountered
    """
    score_ranker = ScoreRanker(sort_order)
    ranked_result = []
    for result in pheval_result:
        if type(result) == PhEvalGeneResult:
            ranked_result.append(
                RankedPhEvalGeneResult.from_gene_result(
                    result, score_ranker.rank_scores(result.score)
                )
            )
        elif type(result) == PhEvalVariantResult:
            ranked_result.append(
                RankedPhEvalVariantResult.from_variant_result(
                    result, score_ranker.rank_scores(result.score)
                )
            )
        elif type(result) == PhEvalDiseaseResult:
            ranked_result.append(
                RankedPhEvalDiseaseResult.from_disease_result(
                    result, score_ranker.rank_scores(result.score)
                )
            )
        else:
            raise ValueError("Incompatible PhEval result type.")
    return ranked_result


def _return_sort_order(sort_order_str: str) -> SortOrder:
    """
    Convert a string derived from the config file into SortOrder Enum

    Args:
        sort_order_str (str): String representation of the sorting order

    Returns:
        SortOrder: Enum representing the specified sorting order

    Raises:
        ValueError: If an incompatible or unknown sorting method is provided
    """
    try:
        return SortOrder[sort_order_str.upper()]
    except KeyError:
        raise ValueError("Incompatible ordering method specified.")


def _create_pheval_result(pheval_result: [PhEvalResult], sort_order_str: str) -> [PhEvalResult]:
    """
    Create PhEval results with corresponding ranks based on the specified sorting order.

    Args:
        pheval_result ([PhEvalResult]): List of PhEvalResult instances to be processed.
        sort_order_str (str): String representation of the desired sorting order.

    Returns:
        List[PhEvalResult]: PhEval results with ranks assigned.
    """
    sort_order = _return_sort_order(sort_order_str)
    sorted_pheval_result = ResultSorter(pheval_result, sort_order).sort_pheval_results()
    return _rank_pheval_result(sorted_pheval_result, sort_order)


def _write_pheval_gene_result(
    ranked_pheval_result: [PhEvalResult], output_dir: Path, tool_result_path: Path
) -> None:
    """
    Write ranked PhEval gene results to a TSV file

    Args:
        ranked_pheval_result ([PhEvalResult]): List of ranked PhEval gene results
        output_dir (Path): Path to the output directory
        tool_result_path (Path): Path to the tool-specific result file
    """
    ranked_result = pd.DataFrame([data.__dict__ for data in ranked_pheval_result])
    pheval_gene_output = ranked_result.loc[:, ["rank", "score", "gene_symbol", "gene_identifier"]]
    pheval_gene_output.to_csv(
        output_dir.joinpath(
            "pheval_gene_results/" + tool_result_path.stem + "-pheval_gene_result.tsv"
        ),
        sep="\t",
        index=False,
    )


def _write_pheval_variant_result(
    ranked_pheval_result: [PhEvalResult], output_dir: Path, tool_result_path: Path
) -> None:
    """
    Write ranked PhEval variant results to a TSV file

    Args:
        ranked_pheval_result ([PhEvalResult]): List of ranked PhEval gene results
        output_dir (Path): Path to the output directory
        tool_result_path (Path): Path to the tool-specific result file
    """
    ranked_result = pd.DataFrame([data.__dict__ for data in ranked_pheval_result])
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


def _write_pheval_disease_result(
    ranked_pheval_result: [RankedPhEvalDiseaseResult], output_dir: Path, tool_result_path: Path
) -> None:
    """
    Write ranked PhEval disease results to a TSV file

    Args:
        ranked_pheval_result ([PhEvalResult]): List of ranked PhEval gene results
        output_dir (Path): Path to the output directory
        tool_result_path (Path): Path to the tool-specific result file
    """
    ranked_result = pd.DataFrame([data.__dict__ for data in ranked_pheval_result])
    pheval_disease_output = ranked_result.loc[
        :, ["rank", "score", "disease_name", "disease_identifier"]
    ]
    pheval_disease_output.to_csv(
        output_dir.joinpath(
            "pheval_disease_results/" + tool_result_path.stem + "-pheval_disease_result.tsv"
        ),
        sep="\t",
        index=False,
    )


def generate_pheval_result(
    pheval_result: [PhEvalResult],
    sort_order_str: str,
    output_dir: Path,
    tool_result_path: Path,
) -> None:
    """
    Generate PhEval variant, gene or disease TSV result based on input results.

    Args:
        pheval_result ([PhEvalResult]): List of PhEvalResult instances to be processed.
        sort_order_str (str): String representation of the desired sorting order.
        output_dir (Path): Path to the output directory.
        tool_result_path (Path): Path to the tool-specific result file.

    Raises:
        ValueError: If the results are not all the same type or an error occurs during file writing.
    """
    if not pheval_result:
        info_log.warning(f"No results found for {tool_result_path.name}")
        return
    ranked_pheval_result = _create_pheval_result(pheval_result, sort_order_str)
    if all(isinstance(result, RankedPhEvalGeneResult) for result in ranked_pheval_result):
        _write_pheval_gene_result(ranked_pheval_result, output_dir, tool_result_path)
    elif all(isinstance(result, RankedPhEvalVariantResult) for result in ranked_pheval_result):
        _write_pheval_variant_result(ranked_pheval_result, output_dir, tool_result_path)
    elif all(isinstance(result, RankedPhEvalDiseaseResult) for result in ranked_pheval_result):
        _write_pheval_disease_result(ranked_pheval_result, output_dir, tool_result_path)
    else:
        raise ValueError("Results are not all of the same type.")
