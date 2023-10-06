from collections import defaultdict
from pathlib import Path

import pandas as pd

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import VariantPrioritisationResult
from pheval.analyse.rank_stats import RankStats, RankStatsWriter
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalVariantResult
from pheval.utils.file_utils import all_files, files_with_suffix, obtain_closest_file_name
from pheval.utils.phenopacket_utils import GenomicVariant, PhenopacketUtil, phenopacket_reader


class AssessVariantPrioritisation:
    """Assess variant prioritisation."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_variant_results: [RankedPhEvalVariantResult],
        threshold: float,
        score_order: str,
        proband_causative_variants: [GenomicVariant],
    ):
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
        """Record the variant prioritisation rank if found within results."""
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
        """Record the variant prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _assess_variant_with_threshold(
        self, result_entry: pd.Series, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """Record the variant prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _record_matched_variant(
        self, rank_stats: RankStats, standardised_variant_result: pd.Series
    ) -> VariantPrioritisationResult:
        """Return the variant rank result - dealing with the specification of a threshold."""
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
        self, rank_stats: RankStats, rank_records: defaultdict
    ) -> None:
        """Assess variant prioritisation."""
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
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                VariantPrioritisationResult(self.phenopacket_path, variant)
                if variant_match is None
                else variant_match,
                rank_records,
            ).record_rank()


def _obtain_causative_variants(phenopacket_path: Path) -> [GenomicVariant]:
    """Obtain causative variants from a phenopacket."""
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
) -> None:
    """Assess variant prioritisation for a phenopacket"""
    phenopacket_path = obtain_closest_file_name(
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
    ).assess_variant_prioritisation(variant_rank_stats, variant_rank_comparison)


def benchmark_variant_prioritisation(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    variant_rank_comparison: defaultdict,
):
    """Benchmark a directory based on variant prioritisation results."""
    variant_rank_stats = RankStats()
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
        )
    return BenchmarkRunResults(
        results_dir=results_directory_and_input.results_dir,
        ranks=variant_rank_comparison,
        rank_stats=variant_rank_stats,
    )
