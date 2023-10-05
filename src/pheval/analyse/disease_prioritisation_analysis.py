from collections import defaultdict
from pathlib import Path

from pheval.analyse.benchmarking_data import BenchmarkRun
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import DiseasePrioritisationResult
from pheval.analyse.rank_stats import RankStats

from pheval.analyse.rank_stats_writer import RankStatsWriter
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalDiseaseResult
from pheval.utils.file_utils import all_files, files_with_suffix, obtain_closest_file_name
from pheval.utils.phenopacket_utils import PhenopacketUtil, ProbandDisease, phenopacket_reader


class AssessDiseasePrioritisation:
    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_disease_results: [RankedPhEvalDiseaseResult],
        threshold: float,
        score_order: str,
        proband_diseases: [ProbandDisease],
    ):
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
        """Record the disease prioritisation rank if found within results."""
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return DiseasePrioritisationResult(self.phenopacket_path, disease, rank)

    def _assess_disease_with_threshold_ascending_order(
        self,
        result_entry: RankedPhEvalDiseaseResult,
        disease: ProbandDisease,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """Record the disease prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.score):
            return self._record_disease_prioritisation_match(disease, result_entry, rank_stats)

    def _assess_disease_with_threshold(
        self,
        result_entry: RankedPhEvalDiseaseResult,
        disease: ProbandDisease,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """Record the disease prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.score):
            return self._record_disease_prioritisation_match(disease, result_entry, rank_stats)

    def _record_matched_disease(
        self,
        disease: ProbandDisease,
        rank_stats: RankStats,
        standardised_disease_result: RankedPhEvalDiseaseResult,
    ) -> DiseasePrioritisationResult:
        """Return the disease rank result - dealing with the specification of a threshold."""
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
        self, rank_stats: RankStats, rank_records: defaultdict
    ) -> None:
        """Assess disease prioritisation."""
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
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                DiseasePrioritisationResult(self.phenopacket_path, disease)
                if disease_match is None
                else disease_match,
                rank_records,
            ).record_rank()


def _obtain_causative_diseases(phenopacket_path: Path) -> [ProbandDisease]:
    """Obtain known diseases from a phenopacket."""
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
) -> None:
    """Assess disease prioritisation for a phenopacket."""
    phenopacket_path = obtain_closest_file_name(
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
    ).assess_disease_prioritisation(disease_rank_stats, disease_rank_comparison)


def benchmark_disease_prioritisation(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    disease_rank_comparison: defaultdict,
    disease_stats_writer: RankStatsWriter,
):
    """Benchmark a directory based on disease prioritisation results."""
    disease_rank_stats = RankStats()
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
        )
    disease_stats_writer.write_row(results_directory_and_input.results_dir, disease_rank_stats)
    return BenchmarkRun(
        results_dir=results_directory_and_input.results_dir,
        ranks=disease_rank_comparison,
        rank_stats=disease_rank_stats,
    )
