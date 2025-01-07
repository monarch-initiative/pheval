from pathlib import Path

from pheval.analyse.assess_prioritisation_base import AssessPrioritisationBase
from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import RunConfig
from pheval.post_processing.post_processing import RankedPhEvalDiseaseResult
from pheval.utils.file_utils import all_files


class AssessDiseasePrioritisation(AssessPrioritisationBase):
    """Class for assessing disease prioritisation based on thresholds and scoring orders."""

    def assess_disease_prioritisation(
        self,
        standardised_disease_result_path: Path,
        phenopacket_path: Path,
        binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess disease prioritisation.

        This method assesses the prioritisation of diseases based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            standardised_disease_result_path (Path): Path to the standardised disease TSV result.
            phenopacket_path (Path): Path to the phenopacket.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        df = self.conn.execute(
            f"SELECT * FROM {self.table_name} WHERE phenopacket = ? ",
            (phenopacket_path.name,),
        ).fetchdf()
        for _i, row in df.iterrows():
            result = (
                self.conn.execute(
                    (
                        f"SELECT * FROM '{standardised_disease_result_path}' "
                        f"WHERE contains_entity_function(CAST(COALESCE(disease_identifier, '') AS VARCHAR),"
                        f" '{row['disease_identifier']}') "
                        f"OR contains_entity_function(CAST(COALESCE(disease_name, '') AS VARCHAR), "
                        f"'{row['disease_name']}')"
                    )
                    if standardised_disease_result_path.exists()
                    and standardised_disease_result_path.stat().st_size > 0
                    else "SELECT NULL WHERE FALSE"
                )
                .fetchdf()
                .to_dict(orient="records")
            )

            if len(result) > 0:
                disease_match = self._record_matched_entity(RankedPhEvalDiseaseResult(**result[0]))
                relevant_ranks.append(disease_match)
                primary_key = f"{phenopacket_path.name}-{row['disease_identifier']}"
                self.conn.execute(
                    f'UPDATE {self.table_name} SET "{self.column}" = ? WHERE identifier = ?',
                    (disease_match, primary_key),
                )
            elif len(result) == 0:
                relevant_ranks.append(0)
        binary_classification_stats.add_classification(
            (
                self.db_connection.parse_table_into_dataclass(
                    str(standardised_disease_result_path), RankedPhEvalDiseaseResult
                )
                if standardised_disease_result_path.exists()
                else []
            ),
            relevant_ranks,
        )


def assess_phenopacket_disease_prioritisation(
    phenopacket_path: Path,
    run: RunConfig,
    disease_binary_classification_stats: BinaryClassificationStats,
    disease_benchmarker: AssessDiseasePrioritisation,
) -> None:
    """
    Assess disease prioritisation for a Phenopacket by comparing PhEval standardised disease results
    against the recorded causative diseases for a proband in the Phenopacket.

    Args:
        phenopacket_path (Path): Path to the Phenopacket.
        run (RunConfig): Run configuration.
        disease_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        disease_benchmarker (AssessDiseasePrioritisation): AssessDiseasePrioritisation class instance.
    """
    standardised_disease_result_path = run.results_dir.joinpath(
        f"pheval_disease_results/{phenopacket_path.stem}-pheval_disease_result.tsv"
    )
    disease_benchmarker.assess_disease_prioritisation(
        standardised_disease_result_path,
        phenopacket_path,
        disease_binary_classification_stats,
    )


def benchmark_disease_prioritisation(
    benchmark_name: str,
    run: RunConfig,
    score_order: str,
    threshold: float,
):
    """
    Benchmark a directory based on disease prioritisation results.

    Args:
        benchmark_name (str): Name of the benchmark.
        run (RunConfig): Run configuration.
        score_order (str): The order in which scores are arranged.
        threshold (float): Threshold for assessment.

    Returns:
        BenchmarkRunResults: An object containing benchmarking results for disease prioritisation,
        including ranks and rank statistics for the benchmarked directory.
    """
    disease_binary_classification_stats = BinaryClassificationStats()
    db_connection = BenchmarkDBManager(benchmark_name)
    db_connection.initialise()
    disease_benchmarker = AssessDiseasePrioritisation(
        db_connection,
        f"{run.phenopacket_dir.parents[0].name}_disease",
        run.run_identifier,
        threshold,
        score_order,
    )
    for phenopacket_path in all_files(run.phenopacket_dir):
        assess_phenopacket_disease_prioritisation(
            phenopacket_path,
            run,
            disease_binary_classification_stats,
            disease_benchmarker,
        )
    db_connection.close()
    disease_rank_stats = RankStats()
    disease_rank_stats.add_ranks(
        benchmark_name=benchmark_name,
        table_name=f"{run.phenopacket_dir.parents[0].name}_disease",
        column_name=str(run.run_identifier),
    )
    return BenchmarkRunResults(
        rank_stats=disease_rank_stats,
        benchmark_name=run.run_identifier,
        binary_classification_stats=disease_binary_classification_stats,
        phenopacket_dir=run.phenopacket_dir,
    )
