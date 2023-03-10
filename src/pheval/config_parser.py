from dataclasses import dataclass
from pathlib import Path

import yaml
from serde import serde
from serde.yaml import from_yaml


@serde
@dataclass
class InputDirConfig:
    tool: str
    software_directory: Path or None
    software_config: Path or None
    input_directory: Path
    phenotype_only: bool


def parse_input_dir_config(config_path: Path) -> InputDirConfig:
    """Reads the config file."""
    with open(config_path, "r") as config_file:
        config = yaml.safe_load(config_file)
    config_file.close()
    return from_yaml(InputDirConfig, yaml.dump(config))
