from collections import defaultdict
from pathlib import Path
from typing import List

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import VariantPrioritisationResult
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalVariantResult
from pheval.utils.file_utils import (
    all_files,
    files_with_suffix,
    obtain_phenopacket_path_from_pheval_result,
)
from pheval.utils.phenopacket_utils import GenomicVariant, PhenopacketUtil, phenopacket_reader


class AssessVariantPrioritisation:
    """Class for assessing variant prioritisation based on thresholds and scoring orders."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_variant_results: List[RankedPhEvalVariantResult],
        threshold: float,
        score_order: str,
        proband_causative_variants: List[GenomicVariant],
    ):
        """
        Initialise AssessVariantPrioritisation class

        Args:
            phenopacket_path (Path): Path to the phenopacket file
            results_dir (Path): Path to the results directory
            standardised_variant_results (List[RankedPhEvalVariantResult]): List of ranked PhEval variant results
            threshold (float): Threshold for scores
            score_order (str): Score order for results, either ascending or descending
            proband_causative_variants (List[GenomicVariant]): List of proband variants

        """
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_variant_results = standardised_variant_results
        self.threshold = threshold
        self.score_order = score_order
        self.proband_causative_variants = proband_causative_variants

    def _record_variant_prioritisation_match(
        self,
        result_entry: RankedPhEvalVariantResult,
        rank_stats: RankStats,
    ) -> VariantPrioritisationResult:
        """
        Record the variant prioritisation rank if found within the results
        Args:
            result_entry (RankedPhEvalVariantResult): Ranked PhEval variant result entry
            rank_stats (RankStats): RankStats class instance
        Returns:
            VariantPrioritisationResult: Recorded correct variant prioritisation rank result
        """
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return VariantPrioritisationResult(
            self.phenopacket_path,
            GenomicVariant(
                chrom=result_entry.chromosome,
                pos=result_entry.start,
                ref=result_entry.ref,
                alt=result_entry.alt,
            ),
            rank,
        )

    def _assess_variant_with_threshold_ascending_order(
        self, result_entry: RankedPhEvalVariantResult, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """
        Record the variant prioritisation rank if it meets the ascending order threshold.

        This method checks if the variant prioritisation rank meets the ascending order threshold.
        If the score of the result entry is less than the threshold, it records the variant rank.

        Args:
            result_entry (RankedPhEvalVariantResult): Ranked PhEval variant result entry
            rank_stats (RankStats): RankStats class instance

        Returns:
            VariantPrioritisationResult: Recorded correct variant prioritisation rank result
        """
        if float(self.threshold) > float(result_entry.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _assess_variant_with_threshold(
        self, result_entry: RankedPhEvalVariantResult, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """
        Record the variant prioritisation rank if it meets the score threshold.

        This method checks if the variant prioritisation rank meets the score threshold.
        If the score of the result entry is greater than the threshold, it records the variant rank.

        Args:
            result_entry (RankedPhEvalVariantResult): Ranked PhEval variant result entry
            rank_stats (RankStats): RankStats class instance

        Returns:
            VariantPrioritisationResult: Recorded correct variant prioritisation rank result
        """
        if float(self.threshold) < float(result_entry.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _record_matched_variant(
        self, rank_stats: RankStats, standardised_variant_result: RankedPhEvalVariantResult
    ) -> VariantPrioritisationResult:
        """
        Return the variant rank result - handling the specification of a threshold.

        This method determines and returns the variant rank result based on the specified threshold
        and score order. If the threshold is 0.0, it records the variant rank directly.
        Otherwise, it assesses the variant with the threshold based on the score order.

        Args:
            rank_stats (RankStats): RankStats class instance
            standardised_variant_result (RankedPhEvalVariantResult): Ranked PhEval variant result entry

        Returns:
            VariantPrioritisationResult: Recorded correct variant prioritisation rank result
        """
        if float(self.threshold) == 0.0:
            return self._record_variant_prioritisation_match(
                standardised_variant_result, rank_stats
            )
        else:
            return (
                self._assess_variant_with_threshold(standardised_variant_result, rank_stats)
                if self.score_order != "ascending"
                else self._assess_variant_with_threshold_ascending_order(
                    standardised_variant_result, rank_stats
                )
            )

    def assess_variant_prioritisation(
        self,
        rank_stats: RankStats,
        rank_records: defaultdict,
        binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess variant prioritisation.

        This method assesses the prioritisation of variants based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            rank_stats (RankStats): RankStats class instance
            rank_records (defaultdict): A defaultdict to store the correct ranked results.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        for variant in self.proband_causative_variants:
            rank_stats.total += 1
            variant_match = VariantPrioritisationResult(self.phenopacket_path, variant)
            for result in self.standardised_variant_results:
                result_variant = GenomicVariant(
                    chrom=result.chromosome,
                    pos=result.start,
                    ref=result.ref,
                    alt=result.alt,
                )
                if variant == result_variant:
                    variant_match = self._record_matched_variant(rank_stats, result)
                    (
                        relevant_ranks.append(variant_match.rank)
                        if variant_match
                        else relevant_ranks.append(0)
                    )
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                (
                    VariantPrioritisationResult(self.phenopacket_path, variant)
                    if variant_match is None
                    else variant_match
                ),
                rank_records,
            ).record_rank()
        rank_stats.relevant_result_ranks.append(relevant_ranks)
        binary_classification_stats.add_classification(
            self.standardised_variant_results, relevant_ranks
        )


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


def assess_phenopacket_variant_prioritisation(
    standardised_variant_result: Path,
    score_order: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    variant_rank_stats: RankStats,
    variant_rank_comparison: defaultdict,
    variant_binary_classification_stats: BinaryClassificationStats,
) -> None:
    """
    Assess variant prioritisation for a Phenopacket by comparing PhEval standardised variant results
    against the recorded causative variants for a proband in the Phenopacket.

    Args:
        standardised_variant_result (Path): Path to the PhEval standardised variant result file.
        score_order (str): The order in which scores are arranged, either ascending or descending.
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories.
        threshold (float): Threshold for assessment.
        variant_rank_stats (RankStats): RankStats class instance.
        variant_rank_comparison (defaultdict): Default dictionary for variant rank comparisons.
        variant_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
    """
    phenopacket_path = obtain_phenopacket_path_from_pheval_result(
        standardised_variant_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    proband_causative_variants = _obtain_causative_variants(phenopacket_path)
    pheval_variant_result = read_standardised_result(standardised_variant_result)
    AssessVariantPrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_variant_results/"),
        parse_pheval_result(RankedPhEvalVariantResult, pheval_variant_result),
        threshold,
        score_order,
        proband_causative_variants,
    ).assess_variant_prioritisation(
        variant_rank_stats, variant_rank_comparison, variant_binary_classification_stats
    )


def benchmark_variant_prioritisation(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    variant_rank_comparison: defaultdict,
):
    """
    Benchmark a directory based on variant prioritisation results.

    Args:
        results_directory_and_input (TrackInputOutputDirectories): Input and output directories.
        score_order (str): The order in which scores are arranged.
        threshold (float): Threshold for assessment.
        variant_rank_comparison (defaultdict): Default dictionary for variant rank comparisons.

    Returns:
        BenchmarkRunResults: An object containing benchmarking results for variant prioritisation,
        including ranks and rank statistics for the benchmarked directory.
    """
    variant_rank_stats = RankStats()
    variant_binary_classification_stats = BinaryClassificationStats()
    for standardised_result in files_with_suffix(
        results_directory_and_input.results_dir.joinpath("pheval_variant_results/"),
        ".tsv",
    ):
        assess_phenopacket_variant_prioritisation(
            standardised_result,
            score_order,
            results_directory_and_input,
            threshold,
            variant_rank_stats,
            variant_rank_comparison,
            variant_binary_classification_stats,
        )
    return BenchmarkRunResults(
        results_dir=results_directory_and_input.results_dir,
        ranks=variant_rank_comparison,
        rank_stats=variant_rank_stats,
        binary_classification_stats=variant_binary_classification_stats,
    )
