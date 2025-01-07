from pathlib import Path

from pheval.analyse.assess_prioritisation_base import AssessPrioritisationBase
from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import RunConfig
from pheval.post_processing.post_processing import RankedPhEvalVariantResult
from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import GenomicVariant


class AssessVariantPrioritisation(AssessPrioritisationBase):
    """Class for assessing variant prioritisation based on thresholds and scoring orders."""

    def assess_variant_prioritisation(
        self,
        standardised_variant_result_path: Path,
        phenopacket_path: Path,
        binary_classification_stats: BinaryClassificationStats,
    ) -> None:
        """
        Assess variant prioritisation.

        This method assesses the prioritisation of variants based on the provided criteria
        and records ranks using a PrioritisationRankRecorder.

        Args:
            standardised_variant_result_path (Path): Path to standardised variant TSV result.
            phenopacket_path (Path): Path to the phenopacket.
            binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        """
        relevant_ranks = []
        df = self.conn.execute(
            f"""SELECT * FROM {self.table_name} WHERE phenopacket = '{phenopacket_path.name}'"""
        ).fetchdf()
        for _i, row in df.iterrows():
            causative_variant = GenomicVariant(
                chrom=row["chrom"],
                pos=int(row["pos"]),
                ref=row["ref"],
                alt=row["alt"],
            )
            result = (
                self.conn.execute(
                    (
                        f"SELECT * FROM '{standardised_variant_result_path}' "
                        f"WHERE "
                        f"chromosome == '{causative_variant.chrom}' AND "
                        f"start == {causative_variant.pos} AND "
                        f"ref == '{causative_variant.ref}' AND "
                        f"alt == '{causative_variant.alt}'"
                    )
                    if standardised_variant_result_path.exists()
                    else "SELECT NULL WHERE FALSE"
                )
                .fetchdf()
                .to_dict(orient="records")
            )

            if len(result) > 0:
                variant_match = self._record_matched_entity(RankedPhEvalVariantResult(**result[0]))
                relevant_ranks.append(variant_match)
                primary_key = (
                    f"{phenopacket_path.name}-{causative_variant.chrom}-{causative_variant.pos}-"
                    f"{causative_variant.ref}-{causative_variant.alt}"
                )
                self.conn.execute(
                    f'UPDATE {self.table_name} SET "{self.column}" = ? WHERE identifier = ?',
                    (variant_match, primary_key),
                )
            elif len(result) == 0:
                relevant_ranks.append(0)
        binary_classification_stats.add_classification(
            (
                self.db_connection.parse_table_into_dataclass(
                    str(standardised_variant_result_path), RankedPhEvalVariantResult
                )
                if standardised_variant_result_path.exists()
                else []
            ),
            relevant_ranks,
        )


def assess_phenopacket_variant_prioritisation(
    phenopacket_path: Path,
    run: RunConfig,
    variant_binary_classification_stats: BinaryClassificationStats,
    variant_benchmarker: AssessVariantPrioritisation,
) -> None:
    """
    Assess variant prioritisation for a Phenopacket by comparing PhEval standardised variant results
    against the recorded causative variants for a proband in the Phenopacket.

    Args:
        phenopacket_path (Path): Path to the Phenopacket.
        run (RunConfig): Run configuration.
        variant_binary_classification_stats (BinaryClassificationStats): BinaryClassificationStats class instance.
        variant_benchmarker (AssessVariantPrioritisation): AssessVariantPrioritisation class instance.
    """
    standardised_variant_result_path = run.results_dir.joinpath(
        f"pheval_variant_results/{phenopacket_path.stem}-pheval_variant_result.tsv"
    )
    variant_benchmarker.assess_variant_prioritisation(
        standardised_variant_result_path,
        phenopacket_path,
        variant_binary_classification_stats,
    )


def benchmark_variant_prioritisation(
    benchmark_name: str,
    run: RunConfig,
    score_order: str,
    threshold: float,
):
    """
    Benchmark a directory based on variant prioritisation results.

    Args:
        benchmark_name (str): Name of the benchmark.
        run (RunConfig): Run configuration.
        score_order (str): The order in which scores are arranged.
        threshold (float): Threshold for assessment.

    Returns:
        BenchmarkRunResults: An object containing benchmarking results for variant prioritisation,
        including ranks and rank statistics for the benchmarked directory.
    """
    variant_binary_classification_stats = BinaryClassificationStats()
    db_connection = BenchmarkDBManager(benchmark_name)
    variant_benchmarker = AssessVariantPrioritisation(
        db_connection,
        f"{run.phenopacket_dir.parents[0].name}" f"_variant",
        run.run_identifier,
        threshold,
        score_order,
    )
    for phenopacket_path in all_files(run.phenopacket_dir):
        assess_phenopacket_variant_prioritisation(
            phenopacket_path,
            run,
            variant_binary_classification_stats,
            variant_benchmarker,
        )
    variant_rank_stats = RankStats()
    variant_rank_stats.add_ranks(
        benchmark_name=benchmark_name,
        table_name=f"{run.phenopacket_dir.parents[0].name}_variant",
        column_name=str(run.run_identifier),
    )
    return BenchmarkRunResults(
        benchmark_name=run.run_identifier,
        rank_stats=variant_rank_stats,
        binary_classification_stats=variant_binary_classification_stats,
        phenopacket_dir=run.phenopacket_dir,
    )
