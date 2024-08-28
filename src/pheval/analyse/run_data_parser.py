from pathlib import Path
from typing import List, Optional

import yaml
from pydantic import BaseModel, root_validator


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

    @root_validator(pre=True)
    def handle_blank_fields(cls, values: dict) -> dict:  # noqa: N805
        """
        Root validator to handle fields that may be explicitly set to None.

        This method checks if 'threshold' and 'score_order' are None and assigns default values if so.

        Args:
            values (dict): The input values provided to the model.

        Returns:
            dict: The updated values with defaults applied where necessary.
        """
        if values.get("threshold") is None:
            values["threshold"] = 0
            print("setting default threshold")
        if values.get("score_order") is None:
            values["score_order"] = "descending"
        return values


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

    @root_validator(pre=True)
    def handle_blank_fields(cls, values: dict) -> dict:  # noqa: N805
        """
        Root validator to handle fields that may be explicitly set to None.

        This method checks if 'plot_type' is None and assigns default value if so.

        Args:
            values (dict): The input values provided to the model.

        Returns:
            dict: The updated values with defaults applied where necessary.
        """
        if values.get("plot_type") is None:
            values["plot_type"] = "bar_cumulative"
        return values


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
    with open(run_config, "r") as f:
        config_data = yaml.safe_load(f)
    f.close()
    config = Config(**config_data)
    return config
