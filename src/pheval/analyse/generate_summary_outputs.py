import itertools
from typing import List

from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.benchmark_generator import BenchmarkRunOutputGenerator
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.generate_plots import generate_plots


def get_new_table_name(run_identifier_1: str, run_identifier_2: str, output_prefix: str) -> str:
    """
    Get the new table name for rank comparison tables.
    Args:
        run_identifier_1: The first run identifier.
        run_identifier_2: The second run identifier.
        output_prefix: The output prefix of the table
    Returns:
        The new table name.
    """
    return f"{run_identifier_1}_vs_" f"{run_identifier_2}_" f"{output_prefix}_rank_comparison"


def create_comparison_table(
    comparison_table_name: str,
    connector: BenchmarkDBManager,
    drop_columns: List[str],
    run_identifier_1: str,
    run_identifier_2: str,
    table_name: str,
) -> None:
    """
    Create rank comparison tables.
    Args:
        comparison_table_name (str): Name of the comparison table to create.
        connector (BenchmarkDBManager): DBConnector instance.
        drop_columns (List[str]): List of columns to drop.
        run_identifier_1 (str): The first run identifier.
        run_identifier_2 (str): The second run identifier.
        table_name (str): Name of the table to extract ranks from
    """
    connector.drop_table(comparison_table_name)
    excluded_columns = tuple(drop_columns + ["identifier"]) if drop_columns else ("identifier",)
    connector.conn.execute(
        f'CREATE TABLE "{comparison_table_name}" AS SELECT * '
        f"EXCLUDE {excluded_columns} FROM {table_name}"
    )

    connector.conn.execute(
        f"""ALTER TABLE "{comparison_table_name}" ADD COLUMN rank_change VARCHAR;"""
    )
    connector.conn.execute(
        f'UPDATE "{comparison_table_name}" SET rank_change = CASE WHEN "{run_identifier_1}" = 0 '
        f'AND "{run_identifier_2}" != 0 '
        f"THEN 'GAINED' WHEN \"{run_identifier_1}\" != 0 AND \"{run_identifier_2}\" = 0 THEN 'LOST' ELSE "
        f'CAST ("{run_identifier_1}" - "{run_identifier_2}" AS VARCHAR) END;'
    )
    connector.conn.commit()


def generate_benchmark_comparison_output(
    benchmark_name: str,
    benchmarking_results: List[BenchmarkRunResults],
    run_identifiers: List[str],
    benchmark_generator: BenchmarkRunOutputGenerator,
    table_name: str,
) -> None:
    """
    Generate prioritisation outputs for benchmarking multiple runs.

    This function generates comparison outputs for benchmarking multiple runs. It compares the results
    between pairs of `BenchmarkRunResults` instances in `benchmarking_results` and generates rank
    comparison outputs using `RankComparisonGenerator` for each pair.

    Args:
        benchmark_name (str): Name of the benchmark.
        benchmarking_results (List[BenchmarkRunResults]): A list containing BenchmarkRunResults instances
            representing the benchmarking results of multiple runs.
        run_identifiers (List[str]): A list of run identifiers.
        benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
        table_name (str): The name of the table where ranks are stored.
    """
    output_prefix = benchmark_generator.prioritisation_type_string
    connector = BenchmarkDBManager(benchmark_name)
    for pair in itertools.combinations(
        [str(result.benchmark_name) for result in benchmarking_results], 2
    ):
        run_identifier_1 = pair[0]
        run_identifier_2 = pair[1]
        drop_columns = [run for run in run_identifiers if run not in pair]
        comparison_table_name = get_new_table_name(
            run_identifier_1, run_identifier_2, output_prefix
        )
        create_comparison_table(
            comparison_table_name,
            connector,
            drop_columns,
            run_identifier_1,
            run_identifier_2,
            table_name,
        )
    generate_plots(
        benchmark_name,
        benchmarking_results,
        benchmark_generator,
    )
