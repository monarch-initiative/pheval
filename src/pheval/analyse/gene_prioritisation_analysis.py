from collections import defaultdict
from pathlib import Path

import pandas as pd

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import GenePrioritisationResult
from pheval.analyse.rank_stats import RankStats, RankStatsWriter
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalGeneResult
from pheval.utils.file_utils import all_files, files_with_suffix, obtain_closest_file_name
from pheval.utils.phenopacket_utils import PhenopacketUtil, ProbandCausativeGene, phenopacket_reader


class AssessGenePrioritisation:
    """Assess gene prioritisation."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_gene_results: [RankedPhEvalGeneResult],
        threshold: float,
        score_order: str,
        proband_causative_genes: [ProbandCausativeGene],
    ):
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_gene_results = standardised_gene_results
        self.threshold = threshold
        self.score_order = score_order
        self.proband_causative_genes = proband_causative_genes

    def _record_gene_prioritisation_match(
        self,
        gene: ProbandCausativeGene,
        result_entry: RankedPhEvalGeneResult,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if found within results."""
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol, rank)

    def _assess_gene_with_threshold_ascending_order(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _assess_gene_with_threshold(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _record_matched_gene(
        self, gene: ProbandCausativeGene, rank_stats: RankStats, standardised_gene_result: pd.Series
    ) -> GenePrioritisationResult:
        """Return the gene rank result - dealing with the specification of a threshold."""
        if float(self.threshold) == 0.0:
            return self._record_gene_prioritisation_match(
                gene, standardised_gene_result, rank_stats
            )
        else:
            return (
                self._assess_gene_with_threshold(standardised_gene_result, gene, rank_stats)
                if self.score_order != "ascending"
                else self._assess_gene_with_threshold_ascending_order(
                    standardised_gene_result, gene, rank_stats
                )
            )

    def assess_gene_prioritisation(self, rank_stats: RankStats, rank_records: defaultdict) -> None:
        """Assess gene prioritisation."""
        for gene in self.proband_causative_genes:
            rank_stats.total += 1
            gene_match = GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol)
            for standardised_gene_result in self.standardised_gene_results:
                if (
                    gene.gene_identifier == standardised_gene_result.gene_identifier
                    or gene.gene_symbol == standardised_gene_result.gene_symbol
                ):
                    gene_match = self._record_matched_gene(
                        gene, rank_stats, standardised_gene_result
                    )
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol)
                if gene_match is None
                else gene_match,
                rank_records,
            ).record_rank()


def _obtain_causative_genes(phenopacket_path: Path) -> [ProbandCausativeGene]:
    """Obtain causative genes from a phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_genes()


def assess_phenopacket_gene_prioritisation(
    standardised_gene_result: Path,
    score_order: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    gene_rank_stats: RankStats,
    gene_rank_comparison: defaultdict,
) -> None:
    """Assess gene prioritisation for a phenopacket."""
    phenopacket_path = obtain_closest_file_name(
        standardised_gene_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    pheval_gene_result = read_standardised_result(standardised_gene_result)
    proband_causative_genes = _obtain_causative_genes(phenopacket_path)
    AssessGenePrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_gene_results/"),
        parse_pheval_result(RankedPhEvalGeneResult, pheval_gene_result),
        threshold,
        score_order,
        proband_causative_genes,
    ).assess_gene_prioritisation(gene_rank_stats, gene_rank_comparison)


def benchmark_gene_prioritisation(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    gene_rank_comparison: defaultdict,
    gene_stats_writer: RankStatsWriter,
) -> BenchmarkRunResults:
    """Benchmark a directory based on gene prioritisation results."""
    gene_rank_stats = RankStats()
    for standardised_result in files_with_suffix(
        results_directory_and_input.results_dir.joinpath("pheval_gene_results/"), ".tsv"
    ):
        assess_phenopacket_gene_prioritisation(
            standardised_result,
            score_order,
            results_directory_and_input,
            threshold,
            gene_rank_stats,
            gene_rank_comparison,
        )
    gene_stats_writer.write_row(results_directory_and_input.results_dir, gene_rank_stats)
    return BenchmarkRunResults(
        results_dir=results_directory_and_input.results_dir,
        ranks=gene_rank_comparison,
        rank_stats=gene_rank_stats,
    )
