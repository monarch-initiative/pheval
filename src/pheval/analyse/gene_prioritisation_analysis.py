from pathlib import Path

from pheval.analyse.assess_prioritisation_base import AssessPrioritisationBase
from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import RunConfig
from pheval.post_processing.post_processing import RankedPhEvalGeneResult
from pheval.utils.file_utils import all_files


class AssessGenePrioritisation(AssessPrioritisationBase):
    """Class for assessing gene prioritisation based on thresholds and scoring orders."""

    def assess_gene_prioritisation(
        self,
        standardised_gene_result_path: Path,
        phenopacket_path: Path,
        binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess gene prioritisation.
        This method assesses the prioritisation of genes based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            standardised_gene_result_path (Path): Path to the standardised gene TSV result.
            phenopacket_path (Path): Path to the Phenopacket.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        df = self.conn.execute(
            f"""SELECT * FROM {self.table_name} WHERE phenopacket = '{phenopacket_path.name}'"""
        ).fetchdf()
        for _i, row in df.iterrows():
            result = (
                self.conn.execute(
                    (
                        f"SELECT * FROM '{standardised_gene_result_path}' "
                        f"WHERE contains_entity_function(CAST(COALESCE(gene_identifier, '') AS VARCHAR), "
                        f"'{row['gene_identifier']}') "
                        f"OR contains_entity_function(CAST(COALESCE(gene_symbol, '') AS VARCHAR), "
                        f"'{row['gene_symbol']}')"
                    )
                    if standardised_gene_result_path.exists()
                    and standardised_gene_result_path.stat().st_size > 0
                    else "SELECT NULL WHERE FALSE"
                )
                .fetchdf()
                .to_dict(orient="records")
            )
            if len(result) > 0:
                gene_match = self._record_matched_entity(RankedPhEvalGeneResult(**result[0]))
                relevant_ranks.append(gene_match)
                primary_key = f"{phenopacket_path.name}-{row['gene_symbol']}"
                self.conn.execute(
                    f'UPDATE {self.table_name} SET "{self.column}" = ? WHERE identifier = ?',
                    (gene_match, primary_key),
                )
            if not result:
                relevant_ranks.append(0)
        binary_classification_stats.add_classification(
            (
                self.db_connection.parse_table_into_dataclass(
                    str(standardised_gene_result_path), RankedPhEvalGeneResult
                )
                if standardised_gene_result_path.exists()
                else []
            ),
            relevant_ranks,
        )


def assess_phenopacket_gene_prioritisation(
    phenopacket_path: Path,
    run: RunConfig,
    gene_binary_classification_stats: BinaryClassificationStats,
    gene_benchmarker: AssessGenePrioritisation,
) -> None:
    """
    Assess gene prioritisation for a Phenopacket by comparing PhEval standardised gene results
    against the recorded causative genes for a proband in the Phenopacket.

    Args:
        phenopacket_path (Path): Path to the Phenopacket.
        run (RunConfig): Run configuration.
        gene_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        gene_benchmarker (AssessGenePrioritisation): AssessGenePrioritisation class instance.
    """
    standardised_gene_result_path = run.results_dir.joinpath(
        f"pheval_gene_results/{phenopacket_path.stem}-pheval_gene_result.tsv"
    )
    gene_benchmarker.assess_gene_prioritisation(
        standardised_gene_result_path,
        phenopacket_path,
        gene_binary_classification_stats,
    )


def benchmark_gene_prioritisation(
    benchmark_name: str,
    run: RunConfig,
    score_order: str,
    threshold: float,
) -> BenchmarkRunResults:
    """
    Benchmark a directory based on gene prioritisation results.
     Args:
         benchmark_name (str): Name of the benchmark.
         run (RunConfig): Run configuration.
         score_order (str): The order in which scores are arranged.
         threshold (float): Threshold for assessment.
     Returns:
         BenchmarkRunResults: An object containing benchmarking results for gene prioritisation,
         including ranks and rank statistics for the benchmarked directory.
    """
    gene_binary_classification_stats = BinaryClassificationStats()
    db_connection = BenchmarkDBManager(benchmark_name)
    db_connection.initialise()
    gene_benchmarker = AssessGenePrioritisation(
        db_connection,
        f"{run.phenopacket_dir.parents[0].name}" f"_gene",
        run.run_identifier,
        threshold,
        score_order,
    )
    for phenopacket_path in all_files(run.phenopacket_dir):
        assess_phenopacket_gene_prioritisation(
            phenopacket_path,
            run,
            gene_binary_classification_stats,
            gene_benchmarker,
        )
    db_connection.close()
    gene_rank_stats = RankStats()
    gene_rank_stats.add_ranks(
        benchmark_name=benchmark_name,
        table_name=f"{run.phenopacket_dir.parents[0].name}_gene",
        column_name=str(run.run_identifier),
    )
    return BenchmarkRunResults(
        rank_stats=gene_rank_stats,
        benchmark_name=run.run_identifier,
        binary_classification_stats=gene_binary_classification_stats,
        phenopacket_dir=run.phenopacket_dir,
    )
