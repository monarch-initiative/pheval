from dataclasses import dataclass
from typing import Callable

from pheval.analyse.benchmarking_data import AnalysisResults, TrackRunPrioritisation
from pheval.constants import (
    DISEASE_PLOT_FILE_PREFIX,
    DISEASE_PLOT_Y_LABEL,
    GENE_PLOT_FILE_PREFIX,
    GENE_PLOT_Y_LABEL,
    VARIANT_PLOT_FILE_PREFIX,
    VARIANT_PLOT_Y_LABEL,
)


@dataclass
class BenchmarkPrioritisationOutputGenerator:
    """Base class for recording data required for generating benchmarking outputs."""

    prioritisation_type_file_prefix: str
    y_label: str
    return_function: Callable[[TrackRunPrioritisation], AnalysisResults]


@dataclass
class GeneBenchmarkPrioritisationOutputGenerator(BenchmarkPrioritisationOutputGenerator):
    """Subclass of BenchmarkPrioritisationOutputGenerator specialised
    for producing gene prioritisation benchmarking outputs."""

    prioritisation_type_file_prefix: str = GENE_PLOT_FILE_PREFIX
    y_label: str = GENE_PLOT_Y_LABEL
    return_function: Callable[
        [TrackRunPrioritisation], AnalysisResults
    ] = TrackRunPrioritisation.return_gene


@dataclass
class VariantBenchmarkPrioritisationOutputGenerator(BenchmarkPrioritisationOutputGenerator):
    """Subclass of BenchmarkPrioritisationOutputGenerator specialised
    for producing variant prioritisation benchmarking outputs."""

    prioritisation_type_file_prefix: str = VARIANT_PLOT_FILE_PREFIX
    y_label: str = VARIANT_PLOT_Y_LABEL
    return_function: Callable[
        [TrackRunPrioritisation], AnalysisResults
    ] = TrackRunPrioritisation.return_variant


@dataclass
class DiseaseBenchmarkPrioritisationOutputGenerator(BenchmarkPrioritisationOutputGenerator):
    """Subclass of BenchmarkPrioritisationOutputGenerator specialised
    for producing disease prioritisation benchmarking outputs."""

    prioritisation_type_file_prefix: str = DISEASE_PLOT_FILE_PREFIX
    y_label: str = DISEASE_PLOT_Y_LABEL
    return_function: Callable[
        [TrackRunPrioritisation], AnalysisResults
    ] = TrackRunPrioritisation.return_disease
