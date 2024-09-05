from dataclasses import dataclass
from typing import Callable

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.disease_prioritisation_analysis import benchmark_disease_prioritisation
from pheval.analyse.gene_prioritisation_analysis import benchmark_gene_prioritisation
from pheval.analyse.run_data_parser import RunConfig, SinglePlotCustomisation
from pheval.analyse.variant_prioritisation_analysis import benchmark_variant_prioritisation


@dataclass
class BenchmarkRunOutputGenerator:
    """Base class for recording data required for generating benchmarking outputs.

    Attributes:
        plot_customisation (SinglePlotCustomisation): Customisation for plot.
        prioritisation_type_string (str):  Prioritisation type string.
        y_label (str): Label for the y-axis in benchmarking outputs.
        generate_benchmark_run_results (Callable): Callable to generate benchmark run results.
            Takes parameters: input and results directory, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file (str): Suffix for the rank comparison file.
    """

    plot_customisation: SinglePlotCustomisation
    prioritisation_type_string: str
    y_label: str
    generate_benchmark_run_results: Callable[[str, RunConfig, str, float], BenchmarkRunResults]
    stats_comparison_file: str


@dataclass
class GeneBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """
    Subclass of BenchmarkRunOutputGenerator specialised
    for producing gene prioritisation benchmarking outputs.

    This subclass inherits from BenchmarkRunOutputGenerator and specialises its attributes
    specifically for gene prioritisation benchmarking.

    Attributes:
        plot_customisation (SinglePlotCustomisation): Customisation for plot.
        prioritisation_type_string (str): Prioritisation type string.
            Defaults to GENE_PRIORITISATION_TYPE_STR.
        y_label (str): Label for the y-axis in gene prioritisation benchmarking outputs.
            Defaults to GENE_PLOT_Y_LABEL.
        generate_benchmark_run_results (Callable): Callable to generate gene prioritisation
            benchmark run results. Defaults to benchmark_gene_prioritisation.
            Takes parameters: run configuration, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file (str): Suffix for the gene rank comparison file.
            Defaults to "-gene_summary".
    """

    plot_customisation: SinglePlotCustomisation = None
    prioritisation_type_string: str = "gene"
    y_label: str = "Disease-causing genes (%)"
    generate_benchmark_run_results: Callable[[str, RunConfig, str, float], BenchmarkRunResults] = (
        benchmark_gene_prioritisation
    )
    stats_comparison_file: str = "gene_summary"


@dataclass
class VariantBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """
    Subclass of BenchmarkRunOutputGenerator specialised
    for producing variant prioritisation benchmarking outputs.

    This subclass inherits from BenchmarkRunOutputGenerator and specialises its attributes
    specifically for variant prioritisation benchmarking.

    Attributes:
        plot_customisation (SinglePlotCustomisation): Customisation for plot.
        prioritisation_type_string (str): Prioritisation type string.
            Defaults to VARIANT_PRIORITISATION_TYPE_STR.
        y_label (str): Label for the y-axis in variant prioritisation benchmarking outputs.
            Defaults to VARIANT_PLOT_Y_LABEL.
        generate_benchmark_run_results (Callable): Callable to generate variant prioritisation
            benchmark run results. Defaults to benchmark_variant_prioritisation.
            Takes parameters: run configuration, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file (str): Suffix for the variant rank comparison file.
            Defaults to "-variant_summary".

    """

    plot_customisation: SinglePlotCustomisation = None
    prioritisation_type_string: str = "variant"
    y_label: str = "Disease-causing variants (%)"
    generate_benchmark_run_results: Callable[[str, RunConfig, str, float], BenchmarkRunResults] = (
        benchmark_variant_prioritisation
    )
    stats_comparison_file: str = "variant_summary"


@dataclass
class DiseaseBenchmarkRunOutputGenerator(BenchmarkRunOutputGenerator):
    """
    Subclass of BenchmarkRunOutputGenerator specialised
    for producing disease prioritisation benchmarking outputs.

    This subclass inherits from BenchmarkRunOutputGenerator and specialises its attributes
    specifically for disease prioritisation benchmarking.

    Attributes:
        plot_customisation (SinglePlotCustomisation): Customisation for plot.
        prioritisation_type_string (str): Prioritisation type string.
            Defaults to DISEASE_PRIORITISATION_TYPE_STR.
        y_label (str): Label for the y-axis in disease prioritisation benchmarking outputs.
            Defaults to DISEASE_PLOT_Y_LABEL.
        generate_benchmark_run_results (Callable): Callable to generate disease prioritisation
            benchmark run results. Defaults to benchmark_disease_prioritisation.
            Takes parameters: run configuration, score order, threshold, rank comparison,
            and returns BenchmarkRunResults.
        stats_comparison_file (str): Suffix for the disease rank comparison file.
            Defaults to "-disease_summary".
    """

    plot_customisation: SinglePlotCustomisation = None
    prioritisation_type_string: str = "disease"
    y_label: str = "Known diseases (%)"
    generate_benchmark_run_results: Callable[[str, RunConfig, str, float], BenchmarkRunResults] = (
        benchmark_disease_prioritisation
    )
    stats_comparison_file: str = "disease_summary"
