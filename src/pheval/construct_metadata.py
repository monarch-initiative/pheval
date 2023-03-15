import timeit
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from serde import serde

from pheval.config_parser import InputDirConfig
from pheval.utils.constants import RUN_METADATA_FILE


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


@dataclass
class MetaDataRecord:
    run_metadata: BasicOutputRunMetaData
    run_metadata_filepath: Path


def create_run_metadata(
    input_dir_config: InputDirConfig,
    tool_version: str,
    input_dir: Path,
    start_time: int,
    corpus: Path,
    corpus_variant_path: Path,
):
    """Create the basic metadata record for a single PhEval run."""
    return MetaDataRecord(
        run_metadata=BasicOutputRunMetaData(
            tool=input_dir_config.tool,
            tool_version=tool_version,
            config=f"{Path(input_dir).parent.name}-{Path(input_dir).name}",
            run_timestamp=int(timeit.default_timer() - start_time),
            corpus=f"{Path(corpus).parent.name}-{Path(corpus).name}",
        ),
        run_metadata_filepath=Path(corpus_variant_path).joinpath(RUN_METADATA_FILE),
    )
