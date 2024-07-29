import itertools
from pathlib import Path
from typing import List
from pheval.analyse.benchmark_generator import BenchmarkRunOutputGenerator
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.generate_plots import generate_plots
from pheval.constants import RANK_COMPARISON_SUFFIX
from pheval.analyse.get_connection import DBConnector


def generate_benchmark_output(
        benchmarking_results: BenchmarkRunResults,
        plot_type: str,
        benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """
    Generate prioritisation outputs for a single benchmarking run.

    Args:
        benchmarking_results (BenchmarkRunResults): Results of a benchmarking run.
        plot_type (str): Type of plot to generate.
        benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
    """
    results_dir_name = benchmarking_results.results_dir.name
    conn = DBConnector().conn
    conn.execute(
        f"""CREATE TABLE {results_dir_name}_{benchmark_generator.prioritisation_type_string}
        {RANK_COMPARISON_SUFFIX} AS SELECT * EXCLUDE (identifier) FROM 
        {benchmarking_results.phenopacket_dir.parents[0].name}_{benchmark_generator.prioritisation_type_string}""")
    conn.close()
    generate_plots(
        [benchmarking_results],
        benchmark_generator,
        plot_type,
    )


def get_new_table_name(result_dir_1: Path, result_dir_2: Path, output_prefix: str) -> str:
    """
    Get the new table name for rank comparison tables.
    Args:
        result_dir_1: The path to the first result directory.
        result_dir_2: The path to the second result directory.
        output_prefix: The output prefix of the table
    Returns:
        The new table name.
    """
    return (f"{Path(result_dir_1).parents[0].name}_{Path(result_dir_1).name}_vs_"
            f"{Path(result_dir_2).parents[0].name}_{Path(result_dir_2).name}_"
            f"{output_prefix}{RANK_COMPARISON_SUFFIX}")


def create_comparison_table(comparison_table_name: str, connector: DBConnector, drop_columns: List[str],
                            result_dir_1: str, result_dir_2: str, table_name: str) -> None:
    connector.drop_table(comparison_table_name)
    connector.conn.execute(
        f"""CREATE TABLE "{comparison_table_name}" AS SELECT * EXCLUDE 
        ('{", ".join(drop_columns)}', identifier) FROM {table_name}""")
    connector.conn.execute(f"""ALTER TABLE "{comparison_table_name}" ADD COLUMN rank_change VARCHAR;""")
    connector.conn.execute(
        f"""UPDATE "{comparison_table_name}" SET rank_change = CASE WHEN "{result_dir_1}" = 0 
        AND "{result_dir_2}" != 0 THEN 'GAINED' WHEN "{result_dir_1}" != 0 
        AND "{result_dir_2}" = 0 THEN 'LOST' ELSE CAST ("{result_dir_1}" - "{result_dir_2}" AS VARCHAR) END;""")
    connector.conn.commit()


def generate_benchmark_comparison_output(
        benchmarking_results: List[BenchmarkRunResults],
        plot_type: str,
        benchmark_generator: BenchmarkRunOutputGenerator,
        table_name: str
) -> None:
    """
    Generate prioritisation outputs for benchmarking multiple runs.

    This function generates comparison outputs for benchmarking multiple runs. It compares the results
    between pairs of `BenchmarkRunResults` instances in `benchmarking_results` and generates rank
    comparison outputs using `RankComparisonGenerator` for each pair.

    Args:
        benchmarking_results (List[BenchmarkRunResults]): A list containing BenchmarkRunResults instances
            representing the benchmarking results of multiple runs.
        plot_type (str): The type of plot to be generated.
        benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
        table_name (str): The name of the table where ranks are stored.
    """
    output_prefix = benchmark_generator.prioritisation_type_string
    connector = DBConnector()
    run_columns = [column for column in
                   connector.conn.execute(f"PRAGMA table_info('{table_name}');").fetchdf()['name'].to_list()
                   if "/" in column]
    for pair in itertools.combinations([str(result.results_dir) for result in benchmarking_results], 2):
        result_dir_1 = pair[0]
        result_dir_2 = pair[1]
        drop_columns = [run for run in run_columns if run not in pair]
        comparison_table_name = get_new_table_name(result_dir_1, result_dir_2, output_prefix)
        create_comparison_table(comparison_table_name, connector, drop_columns, result_dir_1, result_dir_2, table_name)
    generate_plots(
        benchmarking_results,
        benchmark_generator,
        plot_type,
    )
