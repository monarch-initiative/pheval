from collections import defaultdict
from pathlib import Path
from typing import List

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.get_connection import DBConnector
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.post_processing.post_processing import RankedPhEvalVariantResult
from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import GenomicVariant


class AssessVariantPrioritisation:
    """Class for assessing variant prioritisation based on thresholds and scoring orders."""

    def __init__(
            self,
            db_connection: DBConnector,
            table_name: str,
            results_dir: Path,
            threshold: float,
            score_order: str,
    ):
        """
        Initialise AssessVariantPrioritisation class

        Args:
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

    def _assess_variant_with_threshold_ascending_order(
            self, result_entry: RankedPhEvalVariantResult
    ) -> int:
        """
        Record the variant prioritisation rank if it meets the ascending order threshold.

        This method checks if the variant prioritisation rank meets the ascending order threshold.
        If the score of the result entry is less than the threshold, it records the variant rank.

        Args:
            result_entry (RankedPhEvalVariantResult): Ranked PhEval variant result entry

        Returns:
            int: Recorded variant prioritisation rank
        """
        if float(self.threshold) > float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _assess_variant_with_threshold(
            self, result_entry: RankedPhEvalVariantResult
    ) -> int:
        """
        Record the variant prioritisation rank if it meets the score threshold.

        This method checks if the variant prioritisation rank meets the score threshold.
        If the score of the result entry is greater than the threshold, it records the variant rank.

        Args:
            result_entry (RankedPhEvalVariantResult): Ranked PhEval variant result entry

        Returns:
            int: Recorded variant prioritisation rank
        """
        if float(self.threshold) < float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _record_matched_variant(
            self, standardised_variant_result: RankedPhEvalVariantResult
    ) -> int:
        """
        Return the variant rank result - handling the specification of a threshold.

        This method determines and returns the variant rank result based on the specified threshold
        and score order. If the threshold is 0.0, it records the variant rank directly.
        Otherwise, it assesses the variant with the threshold based on the score order.

        Args:
            standardised_variant_result (RankedPhEvalVariantResult): Ranked PhEval variant result entry

        Returns:
            int: Recorded variant prioritisation rank
        """
        if float(self.threshold) == 0.0:
            return standardised_variant_result.rank
        else:
            return (
                self._assess_variant_with_threshold(standardised_variant_result)
                if self.score_order != "ascending"
                else self._assess_variant_with_threshold_ascending_order(
                    standardised_variant_result,
                )
            )

    def assess_variant_prioritisation(
            self,
            standardised_variant_results: List[RankedPhEvalVariantResult],
            phenopacket_path: Path,
            binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess variant prioritisation.

        This method assesses the prioritisation of variants based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            standardised_variant_results (List[RankedPhEvalVariantResult]): List of standardised variant results.
            phenopacket_path (Path): Path to the phenopacket.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        df = self.conn.execute(
            f"""SELECT * FROM {self.table_name} WHERE phenopacket = '{phenopacket_path.name}'""").fetchdf()
        for i, row in df.iterrows():
            causative_variant = GenomicVariant(chrom=row["chrom"],
                                               pos=int(row["pos"]),
                                               ref=row["ref"],
                                               alt=row["alt"], )
            generated_matches = list(result for result in standardised_variant_results if
                                     causative_variant == GenomicVariant(chrom=result.chromosome,
                                                                         pos=result.start,
                                                                         alt=result.alt,
                                                                         ref=result.ref, ))
            if len(generated_matches) > 0:
                variant_match = self._record_matched_variant(generated_matches[0])
                relevant_ranks.append(variant_match)
                primary_key = (f"{phenopacket_path.name}-{causative_variant.chrom}-{causative_variant.pos}-"
                               f"{causative_variant.ref}-{causative_variant.alt}")
                self.conn.execute(
                    f'UPDATE {self.table_name} SET "{self.column}" = ? WHERE identifier = ?',
                    (variant_match, primary_key),
                )

        binary_classification_stats.add_classification(
            standardised_variant_results, relevant_ranks
        )


def assess_phenopacket_variant_prioritisation(
        phenopacket_path: Path,
        results_dir_and_input: TrackInputOutputDirectories,
        variant_binary_classification_stats: BinaryClassificationStats,
        variant_benchmarker: AssessVariantPrioritisation
) -> None:
    """
    Assess variant prioritisation for a Phenopacket by comparing PhEval standardised variant results
    against the recorded causative variants for a proband in the Phenopacket.

    Args:
        phenopacket_path (Path): Path to the Phenopacket.
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories.
        variant_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        variant_benchmarker (AssessVariantPrioritisation): AssessVariantPrioritisation class instance.
    """
    standardised_variant_result = results_dir_and_input.results_dir.joinpath(
        f"pheval_variant_results/{phenopacket_path.stem}-pheval_variant_result.tsv"
    )
    pheval_variant_result = read_standardised_result(standardised_variant_result)
    variant_benchmarker.assess_variant_prioritisation(
        parse_pheval_result(RankedPhEvalVariantResult, pheval_variant_result),
        phenopacket_path,
        variant_binary_classification_stats)


def benchmark_variant_prioritisation(
        results_directory_and_input: TrackInputOutputDirectories,
        score_order: str,
        threshold: float,
):
    """
    Benchmark a directory based on variant prioritisation results.

    Args:
        results_directory_and_input (TrackInputOutputDirectories): Input and output directories.
        score_order (str): The order in which scores are arranged.
        threshold (float): Threshold for assessment.

    Returns:
        BenchmarkRunResults: An object containing benchmarking results for variant prioritisation,
        including ranks and rank statistics for the benchmarked directory.
    """
    variant_binary_classification_stats = BinaryClassificationStats()
    db_connection = DBConnector()
    variant_benchmarker = AssessVariantPrioritisation(db_connection,
                                                      f"{results_directory_and_input.phenopacket_dir.parents[0].name}"
                                                      f"_variant",
                                                      results_directory_and_input.results_dir.joinpath(
                                                          "pheval_variant_results/"),
                                                      threshold,
                                                      score_order,
                                                      )
    for phenopacket_path in all_files(results_directory_and_input.phenopacket_dir):
        assess_phenopacket_variant_prioritisation(
            phenopacket_path,
            results_directory_and_input,
            variant_binary_classification_stats,
            variant_benchmarker
        )
    variant_rank_stats = RankStats()
    variant_rank_stats.add_ranks(
        table_name=f'{results_directory_and_input.phenopacket_dir.parents[0].name}_variant',
        column_name=str(results_directory_and_input.results_dir))
    return BenchmarkRunResults(
        results_dir=results_directory_and_input.results_dir,
        rank_stats=variant_rank_stats,
        binary_classification_stats=variant_binary_classification_stats,
        phenopacket_dir=results_directory_and_input.phenopacket_dir,
    )
