from dataclasses import dataclass
from pathlib import Path

import yaml
from serde import serde
from serde.yaml import from_yaml

from pheval.utils.constants import INPUT_DIR_CONFIG


@serde
@dataclass
class InputDirConfig:
    """
    Class for defining the fields within the input directory config.

    Args:
        tool (str): Name of the tool implementation (e.g. exomiser/phen2gene)
        software_directory (Path or None): Path to the software directory
        software_config (Path or None): Path to the software configuration file
        input_directory (Path): Path to the tool input data
        phenotype_only (bool): Whether the tool is run with HPO terms only (True) or with variant data (False)

    """

    tool: str
    software_directory: Path or None
    software_config: Path or None
    input_directory: Path
    phenotype_only: bool


def parse_input_dir_config(input_dir: Path) -> InputDirConfig:
    """Reads the config file."""
    with open(input_dir.joinpath(INPUT_DIR_CONFIG), "r") as config_file:
        config = yaml.safe_load(config_file)
    config_file.close()
    return from_yaml(InputDirConfig, yaml.dump(config))
