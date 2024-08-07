from dataclasses import dataclass
from pathlib import Path
from typing import List

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


class Config(BaseModel):
    """
    Store configurations for a runs.
    Attributes:
        runs (List[RunConfig]): The list of run configurations.
    """

    runs: List[RunConfig]


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


@dataclass
class TrackInputOutputDirectories:
    """
    Track the input phenopacket test data for a corresponding pheval output directory.

    Attributes:
        phenopacket_dir (Path): The directory containing input phenopackets.
        results_dir (Path): The directory containing output results from pheval.
    """

    phenopacket_dir: Path
    results_dir: Path
