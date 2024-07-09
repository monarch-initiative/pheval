import ast
import re
from pathlib import Path
from typing import List, Union


from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.prioritisation_result_types import GenePrioritisationResult
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.analyse.get_connection import DBConnector
from pheval.post_processing.post_processing import RankedPhEvalGeneResult
from pheval.utils.file_utils import all_files


class AssessGenePrioritisation:
    """Class for assessing gene prioritisation based on thresholds and scoring orders."""

    def __init__(
            self,
            db_connection: DBConnector,
            table_name: str,
            results_dir: Path,
            threshold: float,
            score_order: str,
    ):
        """
        Initialise AssessGenePrioritisation class.

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
        db_connection.add_column(table_name=table_name, column=self.column, default=0)

    def _assess_gene_with_threshold_ascending_order(
            self,
            result_entry: RankedPhEvalGeneResult,
    ) -> int:
        """
        Record the gene prioritisation rank if it meets the ascending order threshold.
        Args:
            result_entry (RankedPhEvalGeneResult): Ranked PhEval gene result entry
        Returns:
            int: Recorded gene prioritisation rank.
        """
        if float(self.threshold) > float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _assess_gene_with_threshold(
            self,
            result_entry: RankedPhEvalGeneResult,
    ) -> int:
        """
        Record the gene prioritisation rank if it meets the score threshold.
        Args:
            result_entry (RankedPhEvalResult): Ranked PhEval gene result entry

        Returns:
            int: Recorded correct gene prioritisation rank.
        """
        if float(self.threshold) < float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _record_matched_gene(
            self,
            standardised_gene_result: RankedPhEvalGeneResult,
    ) -> int:
        """
        Return the gene rank result - handling the specification of a threshold.
        This method determines and returns the gene rank result based on the specified threshold
        and score order. If the threshold is 0.0, it records the gene rank directly.
        Otherwise, it assesses the gene with the threshold based on the score order.
        Args:
            standardised_gene_result (RankedPhEvalGeneResult): Ranked PhEval gene result entry
        Returns:
            GenePrioritisationResult: Recorded correct gene prioritisation rank result
        """
        if float(self.threshold) == 0.0:
            return standardised_gene_result.rank
        else:
            return (
                self._assess_gene_with_threshold(standardised_gene_result)
                if self.score_order != "ascending"
                else self._assess_gene_with_threshold_ascending_order(
                    standardised_gene_result,
                )
            )

    @staticmethod
    def _check_string_representation(entity: str) -> Union[List[str], str]:
        """
        Check if the input string is a representation of a list and returns the list if true, otherwise the string.

        Args:
            entity (str): The input entity to check.

        Returns:
            Union[List[str], str]: A list if the input string is a list representation, otherwise
            the original string.
        """
        list_pattern = re.compile(r"^\[\s*(?:[^\[\],\s]+(?:\s*,\s*[^\[\],\s]+)*)?\s*\]$")
        if list_pattern.match(str(entity)):
            return ast.literal_eval(entity)
        else:
            return entity

    def assess_gene_prioritisation(
            self,
            standardised_gene_results: List[RankedPhEvalGeneResult],
            phenopacket_path: Path,
            binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess gene prioritisation.
        This method assesses the prioritisation of genes based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            standardised_gene_results (List[RankedPhEvalGeneResult]) List of standardised gene results.
            phenopacket_path (Path): Path to the Phenopacket.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        df = self.conn.execute(
            f"""SELECT * FROM {self.table_name} WHERE phenopacket = '{phenopacket_path.name}'""").fetchdf()
        for i, row in df.iterrows():
            generated_matches = list(
                result for result in standardised_gene_results
                if (
                        isinstance(self._check_string_representation(result.gene_identifier), list)
                        and row["gene_identifier"] in self._check_string_representation(result.gene_identifier)
                        or isinstance(self._check_string_representation(result.gene_identifier), str)
                        and row["gene_identifier"] == self._check_string_representation(result.gene_identifier)
                        or isinstance(self._check_string_representation(result.gene_symbol), list)
                        and row["gene_symbol"] in self._check_string_representation(result.gene_symbol)
                        or isinstance(self._check_string_representation(result.gene_symbol), str)
                        and row["gene_symbol"] == self._check_string_representation(result.gene_symbol)
                )
            )
            if len(generated_matches) > 0:
                gene_match = self._record_matched_gene(generated_matches[0])
                relevant_ranks.append(gene_match)
                primary_key = f"{phenopacket_path.name}-{row['gene_symbol']}"
                self.conn.execute(
                    f'UPDATE {self.table_name} SET "{self.column}" = ? WHERE identifier = ?',
                    (gene_match, primary_key),
                )
        binary_classification_stats.add_classification(
            pheval_results=standardised_gene_results, relevant_ranks=relevant_ranks
        )


def assess_phenopacket_gene_prioritisation(
        phenopacket_path: Path,
        results_dir_and_input: TrackInputOutputDirectories,
        gene_binary_classification_stats: BinaryClassificationStats,
        gene_benchmarker: AssessGenePrioritisation,
) -> None:
    """
    Assess gene prioritisation for a Phenopacket by comparing PhEval standardised gene results
    against the recorded causative genes for a proband in the Phenopacket.

    Args:
        phenopacket_path (Path): Path to the Phenopacket.
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories.
        gene_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        gene_benchmarker (AssessGenePrioritisation): AssessGenePrioritisation class instance.
    """
    standardised_gene_result = results_dir_and_input.results_dir.joinpath(
        f"pheval_gene_results/{phenopacket_path.stem}-pheval_gene_result.tsv"
    )
    pheval_gene_result = read_standardised_result(standardised_gene_result)
    gene_benchmarker.assess_gene_prioritisation(
        parse_pheval_result(RankedPhEvalGeneResult, pheval_gene_result),
        phenopacket_path,
        gene_binary_classification_stats
    )


def benchmark_gene_prioritisation(
        results_directory_and_input: TrackInputOutputDirectories,
        score_order: str,
        threshold: float,
) -> BenchmarkRunResults:
    """
    Benchmark a directory based on gene prioritisation results.
     Args:
         results_directory_and_input (TrackInputOutputDirectories): Input and output directories.
         score_order (str): The order in which scores are arranged.
         threshold (float): Threshold for assessment.
     Returns:
         BenchmarkRunResults: An object containing benchmarking results for gene prioritisation,
         including ranks and rank statistics for the benchmarked directory.
    """
    gene_binary_classification_stats = BinaryClassificationStats()
    db_connection = DBConnector()
    gene_benchmarker = AssessGenePrioritisation(db_connection,
                                                f"{results_directory_and_input.phenopacket_dir.parents[0].name}"
                                                f"_gene",
                                                results_directory_and_input.results_dir.joinpath(
                                                    "pheval_gene_results/"),
                                                threshold,
                                                score_order
                                                )
    for phenopacket_path in all_files(results_directory_and_input.phenopacket_dir):
        assess_phenopacket_gene_prioritisation(
            phenopacket_path,
            results_directory_and_input,
            gene_binary_classification_stats,
            gene_benchmarker
        )
    db_connection.close()
    gene_rank_stats = RankStats()
    gene_rank_stats.add_ranks(
        table_name=f'{results_directory_and_input.phenopacket_dir.parents[0].name}_gene',
        column_name=str(results_directory_and_input.results_dir))
    return BenchmarkRunResults(
        rank_stats=gene_rank_stats,
        results_dir=results_directory_and_input.results_dir,
        binary_classification_stats=gene_binary_classification_stats,
        phenopacket_dir=results_directory_and_input.phenopacket_dir,
    )
