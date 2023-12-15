from dataclasses import dataclass
from pathlib import Path
from typing import List

import pandas as pd


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


def parse_run_data_text_file(run_data_path: Path) -> List[TrackInputOutputDirectories]:
    """
    Parse run data .txt file returning a list of input phenopacket and corresponding output directories.

    Args:
        run_data_path (Path): The path to the run data .txt file.

    Returns:
        List[TrackInputOutputDirectories]: A list of TrackInputOutputDirectories objects, containing
        input test data directories and their corresponding output directories.

    Notes:
        The run data .txt file should be formatted with tab-separated values. Each row should contain
        two columns: the first column representing the input test data phenopacket directory, and
        the second column representing the corresponding run output directory.
    """
    run_data = pd.read_csv(run_data_path, delimiter="\t", header=None)
    run_data_list = []
    for _index, row in run_data.iterrows():
        run_data_list.append(
            TrackInputOutputDirectories(phenopacket_dir=Path(row[0]), results_dir=Path(row[1]))
        )
    return run_data_list
