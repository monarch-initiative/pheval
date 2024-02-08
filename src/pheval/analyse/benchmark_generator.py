from collections import defaultdict
from dataclasses import dataclass
from typing import Callable

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.disease_prioritisation_analysis import benchmark_disease_prioritisation
from pheval.analyse.gene_prioritisation_analysis import benchmark_gene_prioritisation
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.analyse.variant_prioritisation_analysis import benchmark_variant_prioritisation
from pheval.constants import (
    DISEASE_PLOT_FILE_PREFIX,
    DISEASE_PLOT_Y_LABEL,
    GENE_PLOT_FILE_PREFIX,
    GENE_PLOT_Y_LABEL,
    VARIANT_PLOT_FILE_PREFIX,
    VARIANT_PLOT_Y_LABEL,
)


@dataclass
class BenchmarkRunOutputGenerator:
    """Base class for recording data required for generating benchmarking outputs.

    Attributes:
        prioritisation_type_file_prefix (str): Prefix for the prioritisation type output file.
        y_label (str): Label for the y-axis in benchmarking outputs.
        generate_benchmark_run_results (Callable): Callable to generate benchmark run results.
            Takes parameters: input and results directory, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file_suffix (str): Suffix for the rank comparison file.
    """

    prioritisation_type_file_prefix: str
    y_label: str
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict], BenchmarkRunResults
    ]
    stats_comparison_file_suffix: str


@dataclass
class GeneBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """
    Subclass of BenchmarkRunOutputGenerator specialised
    for producing gene prioritisation benchmarking outputs.

    This subclass inherits from BenchmarkRunOutputGenerator and specialises its attributes
    specifically for gene prioritisation benchmarking.

    Attributes:
        prioritisation_type_file_prefix (str): Prefix for the gene prioritisation type file.
            Defaults to GENE_PLOT_FILE_PREFIX.
        y_label (str): Label for the y-axis in gene prioritisation benchmarking outputs.
            Defaults to GENE_PLOT_Y_LABEL.
        generate_benchmark_run_results (Callable): Callable to generate gene prioritisation
            benchmark run results. Defaults to benchmark_gene_prioritisation.
            Takes parameters: input and results directory, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file_suffix (str): Suffix for the gene rank comparison file.
            Defaults to "-gene_summary.tsv".
    """

    prioritisation_type_file_prefix: str = GENE_PLOT_FILE_PREFIX
    y_label: str = GENE_PLOT_Y_LABEL
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict], BenchmarkRunResults
    ] = benchmark_gene_prioritisation
    stats_comparison_file_suffix: str = "-gene_summary.tsv"


@dataclass
class VariantBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """
    Subclass of BenchmarkRunOutputGenerator specialised
    for producing variant prioritisation benchmarking outputs.

    This subclass inherits from BenchmarkRunOutputGenerator and specialises its attributes
    specifically for variant prioritisation benchmarking.

    Attributes:
        prioritisation_type_file_prefix (str): Prefix for the variant prioritisation type file.
            Defaults to VARIANT_PLOT_FILE_PREFIX.
        y_label (str): Label for the y-axis in variant prioritisation benchmarking outputs.
            Defaults to VARIANT_PLOT_Y_LABEL.
        generate_benchmark_run_results (Callable): Callable to generate variant prioritisation
            benchmark run results. Defaults to benchmark_variant_prioritisation.
            Takes parameters: input and results directory, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file_suffix (str): Suffix for the variant rank comparison file.
            Defaults to "-variant_summary.tsv".

    """

    prioritisation_type_file_prefix: str = VARIANT_PLOT_FILE_PREFIX
    y_label: str = VARIANT_PLOT_Y_LABEL
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict], BenchmarkRunResults
    ] = benchmark_variant_prioritisation
    stats_comparison_file_suffix: str = "-variant_summary.tsv"


@dataclass
class DiseaseBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """
    Subclass of BenchmarkRunOutputGenerator specialised
    for producing disease prioritisation benchmarking outputs.

    This subclass inherits from BenchmarkRunOutputGenerator and specialises its attributes
    specifically for disease prioritisation benchmarking.

    Attributes:
        prioritisation_type_file_prefix (str): Prefix for the disease prioritisation type file.
            Defaults to DISEASE_PLOT_FILE_PREFIX.
        y_label (str): Label for the y-axis in disease prioritisation benchmarking outputs.
            Defaults to DISEASE_PLOT_Y_LABEL.
        generate_benchmark_run_results (Callable): Callable to generate disease prioritisation
            benchmark run results. Defaults to benchmark_disease_prioritisation.
            Takes parameters: input and results directory, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file_suffix (str): Suffix for the disease rank comparison file.
            Defaults to "-disease_summary.tsv".
    """

    prioritisation_type_file_prefix: str = DISEASE_PLOT_FILE_PREFIX
    y_label: str = DISEASE_PLOT_Y_LABEL
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict], BenchmarkRunResults
    ] = benchmark_disease_prioritisation
    stats_comparison_file_suffix: str = "-disease_summary.tsv"
