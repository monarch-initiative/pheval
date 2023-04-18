from dataclasses import dataclass
from pathlib import Path
from typing import Any

from serde import serde


@serde
@dataclass
class BasicOutputRunMetaData:
    """Class for defining variables for the run metadata.
    Args:
        tool (str): Name of the tool implementation
        tool_version (str): Version of the tool implementation
        config (Path): Path to the config file located in the input directory
        run_timestamp (int): Time taken for run to complete
        corpus (Path): Path to corpus used in pheval run
        tool_specific_configuration_options (Any): Special field that can be overwritten by tool implementations to
                                                   contain any extra tool specific configurations used in the run
    """

    tool: str
    tool_version: str
    config: Path
    run_timestamp: int
    corpus: Path
    tool_specific_configuration_options: Any = None
