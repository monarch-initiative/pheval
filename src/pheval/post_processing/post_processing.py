from enum import Enum
from pathlib import Path
from typing import Callable, Tuple

import polars as pl

from pheval.post_processing.phenopacket_truth_set import PhenopacketTruthSet
from pheval.post_processing.validate_result_format import ResultSchema, validate_dataframe
from pheval.utils.file_utils import all_files
from pheval.utils.logger import get_logger

logger = get_logger()

executed_results = set()


class ResultType(Enum):
    """Enumeration of the possible result types."""

    GENE = "gene"
    DISEASE = "disease"
    VARIANT = "variant"


class SortOrder(Enum):
    """Enumeration representing sorting orders."""

    ASCENDING = 1
    """Ascending sort order."""
    DESCENDING = 2
    """Descending sort order."""


def _rank_results(results: pl.DataFrame, sort_order: SortOrder) -> pl.DataFrame:
    """
    Rank results with the given sort order.
    Args:
        results (pl.DataFrame): The results to rank.
        sort_order (SortOrder): The sort order to use.
    Returns:
        pl.DataFrame: The ranked results.
    """
    sort_descending = True if sort_order == SortOrder.DESCENDING else False
    has_grouping_id = "grouping_id" in results.columns
    if has_grouping_id:
        results = (
            results.sort("score", descending=sort_descending)
            .with_columns(
                pl.struct(["score", "grouping_id"])
                .rank(method="dense", descending=sort_descending)
                .cast(pl.Int32)
                .alias("min_rank")
            )
            .with_columns(pl.col("min_rank").max().over("score").alias("rank"))
        )
    else:
        results = results.sort("score", descending=sort_descending).with_columns(
            pl.col("score").rank(method="max", descending=sort_descending).alias("rank")
        )

    return results


def _write_results_file(out_file: Path, output_df: pl.DataFrame) -> None:
    """
    Write results to compressed Parquet output.
    Args:
        out_file (Path): Output file to write to.
        output_df (pl.DataFrame): Output dataframe.
    """
    output_df.write_parquet(out_file, compression="zstd")


def _write_gene_result(ranked_results: pl.DataFrame, output_file: Path) -> None:
    """
    Write ranked PhEval gene results to a parquet file.

    Args:
        ranked_results ([PhEvalResult]): List of ranked PhEval gene results.
        output_file (Path): Path to the output file.
    """
    gene_output = ranked_results.select(
        ["rank", "score", "gene_symbol", "gene_identifier", "true_positive"]
    )
    _write_results_file(output_file, gene_output)


def _write_variant_result(ranked_results: pl.DataFrame, output_file: Path) -> None:
    """
    Write ranked PhEval variant results to a parquet file.

    Args:
        ranked_results ([PhEvalResult]): List of ranked PhEval variant results.
        output_file (Path): Path to the output file.
    """
    variant_output = ranked_results.select(
        ["rank", "score", "chrom", "start", "end", "ref", "alt", "variant_id", "true_positive"]
    )
    _write_results_file(output_file, variant_output)


def _write_disease_result(ranked_results: pl.DataFrame, output_file: Path) -> None:
    """
    Write ranked PhEval disease results to a parquet file.

    Args:
        ranked_results ([PhEvalResult]): List of ranked PhEval disease results.
        output_file (Path): Path to the output file.
    """
    disease_output = ranked_results.select(["rank", "score", "disease_identifier", "true_positive"])
    _write_results_file(output_file, disease_output)


def _get_result_type(
    result_type: ResultType, phenopacket_truth_set: PhenopacketTruthSet
) -> Tuple[Callable, Callable]:
    """
    Get the methods for extracting the entity and writing the result for a given result type.
    Args:
        result_type (ResultType): The result type.
        phenopacket_truth_set (PhenopacketTruthSet): The phenotype truth set class instance.
    Returns:
        Tuple[Callable, Callable]: The methods for extracting the entity and the write method.
    """
    match result_type:
        case ResultType.GENE:
            return phenopacket_truth_set.classified_gene, _write_gene_result
        case ResultType.VARIANT:
            return phenopacket_truth_set.classified_variant, _write_variant_result
        case ResultType.DISEASE:
            return phenopacket_truth_set.classified_disease, _write_disease_result


def create_empty_pheval_result(
    phenopacket_dir: Path, output_dir: Path, result_type: ResultType
) -> None:
    """
    Create an empty PhEval result for a given result type (gene, variant, or disease).

    Notes:
        This is necessary because some tools may not generate a result output for certain cases.
        By explicitly creating an empty result, which will contain the known entity with a rank and score of 0,
        we can track and identify false negatives  during benchmarking,
        ensuring that missing predictions are accounted for in the evaluation.

    Args:
        phenopacket_dir (Path): The directory containing the phenopackets.
        output_dir (Path): The output directory.
        result_type (ResultType): The result type.

    """
    if result_type in executed_results:
        return
    executed_results.add(result_type)
    phenopacket_truth_set = PhenopacketTruthSet(phenopacket_dir)
    classify_method, write_method = _get_result_type(result_type, phenopacket_truth_set)
    for file in all_files(phenopacket_dir):
        classified_results = classify_method(file.stem)
        write_method(
            classified_results,
            output_dir.joinpath(f"{file.stem}-{result_type.value}_result.parquet"),
        )


@validate_dataframe(ResultSchema.GENE_RESULT_SCHEMA)
def generate_gene_result(
    results: pl.DataFrame,
    sort_order: SortOrder,
    output_dir: Path,
    result_path: Path,
    phenopacket_dir: Path,
) -> None:
    """
    Generate PhEval gene results to a compressed Parquet output.
    Args:
        results (pl.DataFrame): The gene results.
        sort_order (SortOrder): The sort order to use.
        output_dir (Path): Path to the output directory
        result_path (Path): Path to the tool-specific result file.
        phenopacket_dir (Path): Path to the Phenopacket directory
    """
    output_file = output_dir.joinpath(f"pheval_gene_results/{result_path.stem}-gene_result.parquet")
    create_empty_pheval_result(
        phenopacket_dir, output_dir.joinpath("pheval_gene_results"), ResultType.GENE
    )
    ranked_results = _rank_results(results, sort_order)
    classified_results = PhenopacketTruthSet(phenopacket_dir).merge_gene_results(
        ranked_results, output_file
    )
    _write_gene_result(classified_results, output_file)


@validate_dataframe(ResultSchema.VARIANT_RESULT_SCHEMA)
def generate_variant_result(
    results: pl.DataFrame,
    sort_order: SortOrder,
    output_dir: Path,
    result_path: Path,
    phenopacket_dir: Path,
) -> None:
    """
    Generate PhEval variant results to a compressed Parquet output.
    Args:
        results (pl.DataFrame): The variant results.
        sort_order (SortOrder): The sort order to use.
        output_dir (Path): Path to the output directory
        result_path (Path): Path to the tool-specific result file.
        phenopacket_dir (Path): Path to the Phenopacket directory
    """
    output_file = output_dir.joinpath(
        f"pheval_variant_results/{result_path.stem}-variant_result.parquet"
    )
    create_empty_pheval_result(
        phenopacket_dir, output_dir.joinpath("pheval_variant_results"), ResultType.VARIANT
    )
    ranked_results = _rank_results(results, sort_order).with_columns(
        pl.concat_str(["chrom", "start", "ref", "alt"], separator="-").alias("variant_id")
    )
    classified_results = PhenopacketTruthSet(phenopacket_dir).merge_variant_results(
        ranked_results, output_file
    )
    _write_variant_result(classified_results, output_file)


@validate_dataframe(ResultSchema.DISEASE_RESULT_SCHEMA)
def generate_disease_result(
    results: pl.DataFrame,
    sort_order: SortOrder,
    output_dir: Path,
    result_path: Path,
    phenopacket_dir: Path,
) -> None:
    """
    Generate PhEval disease results to a compressed Parquet output.
    Args:
        results (pl.DataFrame): The disease results.
        sort_order (SortOrder): The sort order to use.
        output_dir (Path): Path to the output directory
        result_path (Path): Path to the tool-specific result file.
        phenopacket_dir (Path): Path to the Phenopacket directory
    """
    output_file = output_dir.joinpath(
        f"pheval_disease_results/{result_path.stem}-disease_result.parquet"
    )
    create_empty_pheval_result(
        phenopacket_dir, output_dir.joinpath("pheval_disease_results"), ResultType.DISEASE
    )
    ranked_results = _rank_results(results, sort_order)
    classified_results = PhenopacketTruthSet(phenopacket_dir).merge_disease_results(
        ranked_results, output_file
    )
    _write_disease_result(classified_results, output_file)
