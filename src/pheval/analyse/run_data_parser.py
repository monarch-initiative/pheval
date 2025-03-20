from pathlib import Path
from typing import List, Optional

import yaml
from pydantic import BaseModel, field_validator

from pheval.utils.logger import get_logger


class RunConfig(BaseModel):
    """
    Store configurations for a run.

    Attributes:
        run_identifier (str): The run identifier.
        phenopacket_dir (str): The path to the phenopacket directory used for generating the results.
        results_dir (str): The path to the result directory.
        gene_analysis (bool): Whether to benchmark gene analysis results.
        variant_analysis (bool): Whether to benchmark variant analysis results.
        disease_analysis (bool): Whether to benchmark disease analysis results.
        threshold (Optional[float]): The threshold to consider for benchmarking.
        score_order (Optional[str]): The order of scores to consider for benchmarking, either ascending or descending.
    """

    run_identifier: str
    phenopacket_dir: Path
    results_dir: Path
    gene_analysis: bool
    variant_analysis: bool
    disease_analysis: bool
    threshold: Optional[float]
    score_order: Optional[str]

    @field_validator("threshold", mode="before")
    @classmethod
    def set_threshold(cls, threshold):
        return threshold or None

    @field_validator("score_order", mode="before")
    @classmethod
    def set_score_order(cls, score_order):
        return score_order or "descending"


class SinglePlotCustomisation(BaseModel):
    """
    Store customisations for plots.

    Attributes:
        plot_type (str): The plot type.
        rank_plot_title (str): The title for the rank summary plot.
        roc_curve_title (str): The title for the roc curve plot.
        precision_recall_title (str): The title for the precision-recall plot.
    """

    plot_type: Optional[str] = "bar_cumulative"
    rank_plot_title: Optional[str]
    roc_curve_title: Optional[str]
    precision_recall_title: Optional[str]

    @field_validator("plot_type", mode="before")
    @classmethod
    def set_plot_type(cls, plot_type):
        return plot_type or "bar_cumulative"


class PlotCustomisation(BaseModel):
    """
    Store customisations for all plots.
    Attributes:
        gene_plots (SinglePlotCustomisation): Customisation for all gene benchmarking plots.
        disease_plots (SinglePlotCustomisation): Customisation for all disease benchmarking plots.
        variant_plots (SinglePlotCustomisation): Customisation for all variant benchmarking plots.
    """

    gene_plots: SinglePlotCustomisation
    disease_plots: SinglePlotCustomisation
    variant_plots: SinglePlotCustomisation


class Config(BaseModel):
    """
    Store configurations for a runs.
    Attributes:
        runs (List[RunConfig]): The list of run configurations.
    """

    benchmark_name: str
    runs: List[RunConfig]
    plot_customisation: PlotCustomisation


def parse_run_config(run_config: Path) -> Config:
    """
    Parse a run configuration yaml file.
    Args:
        run_config (Path): The path to the run data yaml configuration.
    Returns:
        Config: The parsed run configurations.
    """
    logger = get_logger()
    logger.info(f"Loading benchmark configuration from {run_config}")
    with open(run_config, "r") as f:
        config_data = yaml.safe_load(f)
    f.close()
    config = Config(**config_data)
    return config
