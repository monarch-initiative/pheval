from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
from serde import serde
from serde.yaml import from_yaml


@serde
@dataclass
class InputDirConfig:
    """
    Class for defining the fields within the input directory config.

    Args:
        tool (str): Name of the tool implementation (e.g. exomiser/phen2gene)
        tool_version (str): Version of the tool implementation
        phenotype_only (bool): Whether the tool is run with HPO terms only (True) or with variant data (False)
        tool_specific_configuration_options (Any): Tool specific configurations


    """

    tool: str
    tool_version: str
    phenotype_only: bool
    tool_specific_configuration_options: Any


def parse_input_dir_config(input_dir: Path) -> InputDirConfig:
    """Reads the config file."""
    with open(Path(input_dir).joinpath("config.yaml"), "r") as config_file:
        config = yaml.safe_load(config_file)
    config_file.close()
    return from_yaml(InputDirConfig, yaml.dump(config))
