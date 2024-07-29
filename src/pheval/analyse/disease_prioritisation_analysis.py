from pathlib import Path
from typing import List

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.get_connection import DBConnector
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalDiseaseResult
from pheval.utils.file_utils import all_files


class AssessDiseasePrioritisation:
    """Class for assessing disease prioritisation based on thresholds and scoring orders."""

    def __init__(
            self,
            db_connection: DBConnector,
            table_name: str,
            results_dir: Path,
            threshold: float,
            score_order: str,
    ):
        """
        Initialise AssessDiseasePrioritisation class

        Args:
            db_connection (DBConnector): Database connection
            table_name (str): Table name
            results_dir (Path): Path to the results directory
            threshold (float): Threshold for scores
            score_order (str): Score order for results, either ascending or descending

        """
        self.results_dir = results_dir
        self.threshold = threshold
        self.score_order = score_order
        self.conn = db_connection.conn
        self.column = str(self.results_dir.parents[0])
        self.table_name = table_name
        db_connection.add_column_integer_default(table_name=table_name, column=self.column, default=0)

    def _assess_disease_with_threshold_ascending_order(
            self,
            result_entry: RankedPhEvalDiseaseResult,
    ) -> int:
        """
        Record the disease prioritisation rank if it meets the ascending order threshold.

        This method checks if the disease prioritisation rank meets the ascending order threshold.
        If the score of the result entry is less than the threshold, it records the disease rank.

        Args:
            result_entry (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry

        Returns:
            int: Recorded disease prioritisation rank
        """
        if float(self.threshold) > float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _assess_disease_with_threshold(
            self,
            result_entry: RankedPhEvalDiseaseResult,
    ) -> int:
        """
        Record the disease prioritisation rank if it meets the score threshold.

        This method checks if the disease prioritisation rank meets the score threshold.
        If the score of the result entry is greater than the threshold, it records the disease rank.

        Args:
            result_entry (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry

        Returns:
            int: Recorded disease prioritisation rank
        """
        if float(self.threshold) < float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _record_matched_disease(
            self,
            standardised_disease_result: RankedPhEvalDiseaseResult,
    ) -> int:
        """
        Return the disease rank result - handling the specification of a threshold.

        This method determines and returns the disease rank result based on the specified threshold
        and score order. If the threshold is 0.0, it records the disease rank directly.
        Otherwise, it assesses the disease with the threshold based on the score order.

        Args:
            standardised_disease_result (RankedPhEvalDiseaseResult): Ranked PhEval disease result entry

        Returns:
            int: Recorded disease prioritisation rank
        """
        if float(self.threshold) == 0.0:
            return standardised_disease_result.rank
        else:
            return (
                self._assess_disease_with_threshold(standardised_disease_result)
                if self.score_order != "ascending"
                else self._assess_disease_with_threshold_ascending_order(
                    standardised_disease_result,
                )
            )

    def assess_disease_prioritisation(
            self,
            standardised_disease_results: List[RankedPhEvalDiseaseResult],
            phenopacket_path: Path,
            binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess disease prioritisation.

        This method assesses the prioritisation of diseases based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            standardised_disease_results (List[RankedPhEvalDiseaseResult]): List of standardised disease results.
            phenopacket_path (Path): Path to the phenopacket.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        df = self.conn.execute(
            f"""SELECT * FROM {self.table_name} WHERE phenopacket = '{phenopacket_path.name}'""").fetchdf()
        for i, row in df.iterrows():
            generated_matches = list(result for result in standardised_disease_results if
                                     row["disease_name"] == result.disease_name or row[
                                         "disease_identifier"] == result.disease_identifier)
            if len(generated_matches) > 0:
                disease_match = self._record_matched_disease(generated_matches[0])
                relevant_ranks.append(disease_match)
                primary_key = f"{phenopacket_path.name}-{row['disease_identifier']}"
                self.conn.execute(
                    f'UPDATE {self.table_name} SET "{self.column}" = ? WHERE identifier = ?',
                    (disease_match, primary_key),
                )
        binary_classification_stats.add_classification(
            standardised_disease_results, relevant_ranks
        )


def assess_phenopacket_disease_prioritisation(
        phenopacket_path: Path,
        results_dir_and_input: TrackInputOutputDirectories,
        disease_binary_classification_stats: BinaryClassificationStats,
        disease_benchmarker: AssessDiseasePrioritisation
) -> None:
    """
    Assess disease prioritisation for a Phenopacket by comparing PhEval standardised disease results
    against the recorded causative diseases for a proband in the Phenopacket.

    Args:
        phenopacket_path (Path): Path to the Phenopacket.
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories.
        disease_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        disease_benchmarker (AssessDiseasePrioritisation): AssessDiseasePrioritisation class instance.
    """
    standardised_disease_result = results_dir_and_input.results_dir.joinpath(
        f"pheval_disease_results/{phenopacket_path.stem}-pheval_disease_result.tsv"
    )
    pheval_disease_result = read_standardised_result(standardised_disease_result)
    disease_benchmarker.assess_disease_prioritisation(
        parse_pheval_result(RankedPhEvalDiseaseResult, pheval_disease_result),
        phenopacket_path,
        disease_binary_classification_stats)


def benchmark_disease_prioritisation(
        results_directory_and_input: TrackInputOutputDirectories,
        score_order: str,
        threshold: float,
):
    """
    Benchmark a directory based on disease prioritisation results.

    Args:
        results_directory_and_input (TrackInputOutputDirectories): Input and output directories.
        score_order (str): The order in which scores are arranged.
        threshold (float): Threshold for assessment.

    Returns:
        BenchmarkRunResults: An object containing benchmarking results for disease prioritisation,
        including ranks and rank statistics for the benchmarked directory.
    """
    disease_binary_classification_stats = BinaryClassificationStats()
    db_connection = DBConnector()
    disease_benchmarker = AssessDiseasePrioritisation(db_connection,
                                                      f"{results_directory_and_input.phenopacket_dir.parents[0].name}_disease",
                                                      results_directory_and_input.results_dir.joinpath(
                                                          "pheval_disease_results/"),
                                                      threshold,
                                                      score_order)
    for phenopacket_path in all_files(results_directory_and_input.phenopacket_dir):
        assess_phenopacket_disease_prioritisation(
            phenopacket_path,
            results_directory_and_input,
            disease_binary_classification_stats,
            disease_benchmarker
        )
    db_connection.close()
    disease_rank_stats = RankStats()
    disease_rank_stats.add_ranks(
        table_name=f'{results_directory_and_input.phenopacket_dir.parents[0].name}_disease',
        column_name=str(results_directory_and_input.results_dir))
    return BenchmarkRunResults(
        rank_stats=disease_rank_stats,
        results_dir=results_directory_and_input.results_dir,
        binary_classification_stats=disease_binary_classification_stats,
        phenopacket_dir=results_directory_and_input.phenopacket_dir,
    )
