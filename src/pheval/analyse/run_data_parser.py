from pathlib import Path
from typing import List, Optional

import yaml
from pydantic import BaseModel


class RunConfig(BaseModel):
    """
    Store configurations for a run.

    Attributes:
        run_identifier (str): The run identifier.
        phenopacket_dir (str): The path to the phenopacket directory used for generating the results.
        results_dir (str): The path to the results directory.
        gene_analysis (bool): Whether or not to benchmark gene analysis results.
        variant_analysis (bool): Whether or not to benchmark variant analysis results.
        disease_analysis (bool): Whether or not to benchmark disease analysis results.
    """

    run_identifier: str
    phenopacket_dir: Path
    results_dir: Path
    gene_analysis: bool
    variant_analysis: bool
    disease_analysis: bool


class SinglePlotCustomisation(BaseModel):
    """
    Store customisations for plots.

    Attributes:
        plot_type (str): The plot type.
        rank_plot_title (str): The title for the rank summary plot.
        roc_curve_title (str): The title for the roc curve plot.
        precision_recall_title (str): The title for the precision-recall plot.
    """

    plot_type: Optional[str]
    rank_plot_title: Optional[str]
    roc_curve_title: Optional[str]
    precision_recall_title: Optional[str]


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


def parse_run_config(run_data_path: Path) -> Config:
    """
    Parse a run configuration yaml file.
    Args:
        run_data_path (Path): The path to the run data yaml configuration.
    Returns:
        Config: The parsed run configurations.
    """
    with open(run_data_path, "r") as f:
        config_data = yaml.safe_load(f)
    f.close()
    config = Config(**config_data)
    return config
