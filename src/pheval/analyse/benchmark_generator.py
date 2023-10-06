from collections import defaultdict
from dataclasses import dataclass
from typing import Callable

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.disease_prioritisation_analysis import benchmark_disease_prioritisation
from pheval.analyse.gene_prioritisation_analysis import benchmark_gene_prioritisation
from pheval.analyse.rank_stats import RankStatsWriter
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
    """Base class for recording data required for generating benchmarking outputs."""

    prioritisation_type_file_prefix: str
    y_label: str
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict, RankStatsWriter], BenchmarkRunResults
    ]
    rank_comparison_file_suffix: str


@dataclass
class GeneBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """Subclass of BenchmarkRunOutputGenerator specialised
    for producing gene prioritisation benchmarking outputs."""

    prioritisation_type_file_prefix: str = GENE_PLOT_FILE_PREFIX
    y_label: str = GENE_PLOT_Y_LABEL
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict, RankStatsWriter], BenchmarkRunResults
    ] = benchmark_gene_prioritisation
    rank_comparison_file_suffix: str = "-gene_summary.tsv"


@dataclass
class VariantBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """Subclass of BenchmarkRunOutputGenerator specialised
    for producing variant prioritisation benchmarking outputs."""

    prioritisation_type_file_prefix: str = VARIANT_PLOT_FILE_PREFIX
    y_label: str = VARIANT_PLOT_Y_LABEL
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict, RankStatsWriter], BenchmarkRunResults
    ] = benchmark_variant_prioritisation
    rank_comparison_file_suffix: str = "-variant_summary.tsv"


@dataclass
class DiseaseBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """Subclass of BenchmarkRunOutputGenerator specialised
    for producing disease prioritisation benchmarking outputs."""

    prioritisation_type_file_prefix: str = DISEASE_PLOT_FILE_PREFIX
    y_label: str = DISEASE_PLOT_Y_LABEL
    generate_benchmark_run_results: Callable[
        [TrackInputOutputDirectories, str, float, defaultdict, RankStatsWriter], BenchmarkRunResults
    ] = benchmark_disease_prioritisation
    rank_comparison_file_suffix: str = "-disease_summary.tsv"
