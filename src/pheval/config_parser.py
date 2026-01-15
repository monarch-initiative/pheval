from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
from serde import serde
from serde.yaml import from_yaml

from pheval.utils.logger import get_logger

logger = get_logger()

_BOOL_FIELDS = ("variant_analysis", "gene_analysis", "disease_analysis")

@serde
@dataclass
class InputDirConfig:
    """
    Class for defining the fields within the input directory config.

    Args:
        tool (str): Name of the tool implementation (e.g. exomiser/phen2gene)
        tool_version (str): Version of the tool implementation
        variant_analysis (bool): Whether to extract prioritised variants from results.
        gene_analysis (bool): Whether to extract prioritised genes from results.
        disease_analysis (bool): Whether to extract prioritised diseases from results.
        tool_specific_configuration_options (Any): Tool specific configurations

    """

    tool: str
    tool_version: str
    variant_analysis: bool
    gene_analysis: bool
    disease_analysis: bool
    tool_specific_configuration_options: Any

def _enforce_bool_fields(cfg: dict) -> None:
    """
    Validates that specified boolean fields in the configuration dictionary are present and have
    boolean values.

    Args:
        cfg (dict): The configuration dictionary to validate.
    Raises:
        ValueError: If a required boolean field is missing from the configuration
            dictionary or if the value of a boolean field is not a boolean.
    """
    for key in _BOOL_FIELDS:
        val = cfg[key]
        if not isinstance(val, bool):
            raise ValueError(
                f"Invalid value for '{key}': {val!r}. "
                f"Expected a boolean true/false (do not quote)."
            )

def parse_input_dir_config(input_dir: Path) -> InputDirConfig:
    """Reads the config file."""
    logger.info(f"Parsing config.yaml located in {input_dir}.")
    with open(Path(input_dir).joinpath("config.yaml")) as config_file:
        config = yaml.safe_load(config_file)
    config_file.close()
    _enforce_bool_fields(config)
    return from_yaml(InputDirConfig, yaml.dump(config))
