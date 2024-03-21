from collections import defaultdict
from pathlib import Path
from typing import List

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import GenePrioritisationResult
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalGeneResult
from pheval.utils.file_utils import (
    all_files,
    files_with_suffix,
    obtain_phenopacket_path_from_pheval_result,
)
from pheval.utils.phenopacket_utils import PhenopacketUtil, ProbandCausativeGene, phenopacket_reader


class AssessGenePrioritisation:
    """Class for assessing gene prioritisation based on thresholds and scoring orders."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_gene_results: List[RankedPhEvalGeneResult],
        threshold: float,
        score_order: str,
        proband_causative_genes: List[ProbandCausativeGene],
    ):
        """
        Initialise AssessGenePrioritisation class.

        Args:
            phenopacket_path (Path): Path to the phenopacket file
            results_dir (Path): Path to the results directory
            standardised_gene_results (List[RankedPhEvalGeneResult]): List of ranked PhEval gene results
            threshold (float): Threshold for scores
            score_order (str): Score order for results, either ascending or descending
            proband_causative_genes (List[ProbandCausativeGene]): List of proband causative genes
        """
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
        """
        Record the gene prioritisation rank if found within the results

        Args:
            gene (ProbandCausativeGene): Diagnosed proband gene
            result_entry (RankedPhEvalGeneResult): Ranked PhEval gene result entry
            rank_stats (RankStats): RankStats class instance

        Returns:
            GenePrioritisationResult: Recorded correct gene prioritisation rank result
        """
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol, rank)

    def _assess_gene_with_threshold_ascending_order(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """
        Record the gene prioritisation rank if it meets the ascending order threshold.

        This method checks if the gene prioritisation rank meets the ascending order threshold.
        If the score of the result entry is less than the threshold, it records the gene rank.

        Args:
            result_entry (RankedPhEvalGeneResult): Ranked PhEval gene result entry
            gene (ProbandCausativeGene): Diagnosed proband gene
            rank_stats (RankStats): RankStats class instance
        Returns:
            GenePrioritisationResult: Recorded correct gene prioritisation rank result
        """
        if float(self.threshold) > float(result_entry.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _assess_gene_with_threshold(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """
        Record the gene prioritisation rank if it meets the score threshold.
        This method checks if the gene prioritisation rank meets the score threshold.
        If the score of the result entry is greater than the threshold, it records the gene rank.

        Args:
            result_entry (RankedPhEvalResult): Ranked PhEval gene result entry
            gene (ProbandCausativeGene): Diagnosed proband gene
            rank_stats (RankStats): RankStats class instance

        Returns:
            GenePrioritisationResult: Recorded correct gene prioritisation rank result
        """
        if float(self.threshold) < float(result_entry.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _record_matched_gene(
        self,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
        standardised_gene_result: RankedPhEvalGeneResult,
    ) -> GenePrioritisationResult:
        """
        Return the gene rank result - handling the specification of a threshold.
        This method determines and returns the gene rank result based on the specified threshold
        and score order. If the threshold is 0.0, it records the gene rank directly.
        Otherwise, it assesses the gene with the threshold based on the score order.
        Args:
            gene (ProbandCausativeGene): Diagnosed proband gene
            rank_stats (RankStats): RankStats class instance
            standardised_gene_result (RankedPhEvalGeneResult): Ranked PhEval gene result entry
        Returns:
            GenePrioritisationResult: Recorded correct gene prioritisation rank result
        """
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

    def assess_gene_prioritisation(
        self,
        rank_stats: RankStats,
        rank_records: defaultdict,
        binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess gene prioritisation.
        This method assesses the prioritisation of genes based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            rank_stats (RankStats): RankStats class instance
            rank_records (defaultdict): A defaultdict to store the correct ranked results.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
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
                    (
                        relevant_ranks.append(gene_match.rank)
                        if gene_match
                        else relevant_ranks.append(0)
                    )
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                (
                    GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol)
                    if gene_match is None
                    else gene_match
                ),
                rank_records,
            ).record_rank()
        rank_stats.relevant_result_ranks.append(relevant_ranks)
        binary_classification_stats.add_classification(
            pheval_results=self.standardised_gene_results, relevant_ranks=relevant_ranks
        )


def _obtain_causative_genes(phenopacket_path: Path) -> List[ProbandCausativeGene]:
    """
    Obtain known genes from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.
    Returns:
       List[ProbandCausativeGene]: A list of known genes associated with the proband,
       extracted from the Phenopacket.
    """
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
    gene_binary_classification_stats: BinaryClassificationStats,
) -> None:
    """
    Assess gene prioritisation for a Phenopacket by comparing PhEval standardised gene results
    against the recorded causative genes for a proband in the Phenopacket.

    Args:
        standardised_gene_result (Path): Path to the PhEval standardised gene result file.
        score_order (str): The order in which scores are arranged, either ascending or descending.
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories.
        threshold (float): Threshold for assessment.
        gene_rank_stats (RankStats): RankStats class instance.
        gene_rank_comparison (defaultdict): Default dictionary for gene rank comparisons.
        gene_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
    """
    phenopacket_path = obtain_phenopacket_path_from_pheval_result(
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
    ).assess_gene_prioritisation(
        gene_rank_stats, gene_rank_comparison, gene_binary_classification_stats
    )


def benchmark_gene_prioritisation(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    gene_rank_comparison: defaultdict,
) -> BenchmarkRunResults:
    """
    Benchmark a directory based on gene prioritisation results.
     Args:
         results_directory_and_input (TrackInputOutputDirectories): Input and output directories.
         score_order (str): The order in which scores are arranged.
         threshold (float): Threshold for assessment.
         gene_rank_comparison (defaultdict): Default dictionary for gene rank comparisons.
     Returns:
         BenchmarkRunResults: An object containing benchmarking results for gene prioritisation,
         including ranks and rank statistics for the benchmarked directory.
    """
    gene_rank_stats = RankStats()
    gene_binary_classification_stats = BinaryClassificationStats()
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
            gene_binary_classification_stats,
        )
    return BenchmarkRunResults(
        results_dir=results_directory_and_input.results_dir,
        ranks=gene_rank_comparison,
        rank_stats=gene_rank_stats,
        binary_classification_stats=gene_binary_classification_stats,
    )
