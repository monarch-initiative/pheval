from collections import defaultdict
from pathlib import Path
from typing import List

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import DiseasePrioritisationResult
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalDiseaseResult
from pheval.utils.file_utils import (
    all_files,
    files_with_suffix,
    obtain_phenopacket_path_from_pheval_result,
)
from pheval.utils.phenopacket_utils import PhenopacketUtil, ProbandDisease, phenopacket_reader


class AssessDiseasePrioritisation:
    """Class for assessing disease prioritisation based on thresholds and scoring orders."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_disease_results: List[RankedPhEvalDiseaseResult],
        threshold: float,
        score_order: str,
        proband_diseases: List[ProbandDisease],
    ):
        """
        Initialise AssessDiseasePrioritisation class

        Args:
            phenopacket_path (Path): Path to the phenopacket file
            results_dir (Path): Path to the results directory
            standardised_disease_results (List[RankedPhEvalDiseaseResult]): List of ranked PhEval disease results
            threshold (float): Threshold for scores
            score_order (str): Score order for results, either ascending or descending
            proband_diseases (List[ProbandDisease]): List of proband diseases

        """
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_disease_results = standardised_disease_results
        self.threshold = threshold
        self.score_order = score_order
        self.proband_diseases = proband_diseases

    def _record_disease_prioritisation_match(
        self,
        disease: ProbandDisease,
        result_entry: RankedPhEvalDiseaseResult,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """
        Record the disease prioritisation rank if found within the results
        Args:
            disease (ProbandDisease): Diagnosed proband disease
            result_entry (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry
            rank_stats (RankStats): RankStats class instance
        Returns:
            DiseasePrioritisationResult: Recorded correct disease prioritisation rank result
        """
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return DiseasePrioritisationResult(self.phenopacket_path, disease, rank)

    def _assess_disease_with_threshold_ascending_order(
        self,
        result_entry: RankedPhEvalDiseaseResult,
        disease: ProbandDisease,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """
        Record the disease prioritisation rank if it meets the ascending order threshold.

        This method checks if the disease prioritisation rank meets the ascending order threshold.
        If the score of the result entry is less than the threshold, it records the disease rank.

        Args:
            result_entry (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry
            disease (ProbandDisease): Diagnosed proband disease
            rank_stats (RankStats): RankStats class instance

        Returns:
            DiseasePrioritisationResult: Recorded correct disease prioritisation rank result
        """
        if float(self.threshold) > float(result_entry.score):
            return self._record_disease_prioritisation_match(disease, result_entry, rank_stats)

    def _assess_disease_with_threshold(
        self,
        result_entry: RankedPhEvalDiseaseResult,
        disease: ProbandDisease,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """
        Record the disease prioritisation rank if it meets the score threshold.

        This method checks if the disease prioritisation rank meets the score threshold.
        If the score of the result entry is greater than the threshold, it records the disease rank.

        Args:
            result_entry (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry
            disease (ProbandDisease): Diagnosed proband disease
            rank_stats (RankStats): RankStats class instance

        Returns:
            DiseasePrioritisationResult: Recorded correct disease prioritisation rank result
        """
        if float(self.threshold) < float(result_entry.score):
            return self._record_disease_prioritisation_match(disease, result_entry, rank_stats)

    def _record_matched_disease(
        self,
        disease: ProbandDisease,
        rank_stats: RankStats,
        standardised_disease_result: RankedPhEvalDiseaseResult,
    ) -> DiseasePrioritisationResult:
        """
        Return the disease rank result - handling the specification of a threshold.

        This method determines and returns the disease rank result based on the specified threshold
        and score order. If the threshold is 0.0, it records the disease rank directly.
        Otherwise, it assesses the disease with the threshold based on the score order.

        Args:
            disease (ProbandDisease): Diagnosed proband disease
            rank_stats (RankStats): RankStats class instance
            standardised_disease_result (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry

        Returns:
            DiseasePrioritisationResult: Recorded correct disease prioritisation rank result
        """
        if float(self.threshold) == 0.0:
            return self._record_disease_prioritisation_match(
                disease, standardised_disease_result, rank_stats
            )
        else:
            return (
                self._assess_disease_with_threshold(
                    standardised_disease_result, disease, rank_stats
                )
                if self.score_order != "ascending"
                else self._assess_disease_with_threshold_ascending_order(
                    standardised_disease_result, disease, rank_stats
                )
            )

    def assess_disease_prioritisation(
        self,
        rank_stats: RankStats,
        rank_records: defaultdict,
        binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess disease prioritisation.

        This method assesses the prioritisation of diseases based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            rank_stats (RankStats): RankStats class instance
            rank_records (defaultdict): A defaultdict to store the correct ranked results.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        for disease in self.proband_diseases:
            rank_stats.total += 1
            disease_match = DiseasePrioritisationResult(self.phenopacket_path, disease)
            for standardised_disease_result in self.standardised_disease_results:
                if (
                    disease.disease_identifier == standardised_disease_result.disease_identifier
                    or disease.disease_name == standardised_disease_result.disease_name
                ):
                    disease_match = self._record_matched_disease(
                        disease, rank_stats, standardised_disease_result
                    )
                    (
                        relevant_ranks.append(disease_match.rank)
                        if disease_match
                        else relevant_ranks.append(0)
                    )
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                (
                    DiseasePrioritisationResult(self.phenopacket_path, disease)
                    if disease_match is None
                    else disease_match
                ),
                rank_records,
            ).record_rank()
        rank_stats.relevant_result_ranks.append(relevant_ranks)
        binary_classification_stats.add_classification(
            self.standardised_disease_results, relevant_ranks
        )


def _obtain_causative_diseases(phenopacket_path: Path) -> List[ProbandDisease]:
    """
    Obtain known diseases from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.

    Returns:
       List[ProbandDisease]: A list of known diseases associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnoses()


def assess_phenopacket_disease_prioritisation(
    standardised_disease_result: Path,
    score_order: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    disease_rank_stats: RankStats,
    disease_rank_comparison: defaultdict,
    disease_binary_classification_stats: BinaryClassificationStats,
) -> None:
    """
    Assess disease prioritisation for a Phenopacket by comparing PhEval standardised disease results
    against the recorded causative diseases for a proband in the Phenopacket.

    Args:
        standardised_disease_result (Path): Path to the PhEval standardised disease result file.
        score_order (str): The order in which scores are arranged, either ascending or descending.
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories.
        threshold (float): Threshold for assessment.
        disease_rank_stats (RankStats): RankStats class instance.
        disease_rank_comparison (defaultdict): Default dictionary for disease rank comparisons.
        disease_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
    """
    phenopacket_path = obtain_phenopacket_path_from_pheval_result(
        standardised_disease_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    pheval_disease_result = read_standardised_result(standardised_disease_result)
    proband_diseases = _obtain_causative_diseases(phenopacket_path)
    AssessDiseasePrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_disease_results/"),
        parse_pheval_result(RankedPhEvalDiseaseResult, pheval_disease_result),
        threshold,
        score_order,
        proband_diseases,
    ).assess_disease_prioritisation(
        disease_rank_stats, disease_rank_comparison, disease_binary_classification_stats
    )


def benchmark_disease_prioritisation(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    disease_rank_comparison: defaultdict,
):
    """
    Benchmark a directory based on disease prioritisation results.

    Args:
        results_directory_and_input (TrackInputOutputDirectories): Input and output directories.
        score_order (str): The order in which scores are arranged.
        threshold (float): Threshold for assessment.
        disease_rank_comparison (defaultdict): Default dictionary for disease rank comparisons.

    Returns:
        BenchmarkRunResults: An object containing benchmarking results for disease prioritisation,
        including ranks and rank statistics for the benchmarked directory.
    """
    disease_rank_stats = RankStats()
    disease_binary_classification_stats = BinaryClassificationStats()
    for standardised_result in files_with_suffix(
        results_directory_and_input.results_dir.joinpath("pheval_disease_results/"),
        ".tsv",
    ):
        assess_phenopacket_disease_prioritisation(
            standardised_result,
            score_order,
            results_directory_and_input,
            threshold,
            disease_rank_stats,
            disease_rank_comparison,
            disease_binary_classification_stats,
        )
    return BenchmarkRunResults(
        results_dir=results_directory_and_input.results_dir,
        ranks=disease_rank_comparison,
        rank_stats=disease_rank_stats,
        binary_classification_stats=disease_binary_classification_stats,
    )
