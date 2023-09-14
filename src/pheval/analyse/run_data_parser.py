from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class TrackInputOutputDirectories:
    """Track the input testdata for a corresponding pheval output directory"""

    phenopacket_dir: Path
    results_dir: Path


def _parse_run_data_text_file(run_data_path: Path) -> [TrackInputOutputDirectories]:
    """Parse run data .txt file returning a list of input testdata and corresponding output directories."""
    run_data = pd.read_csv(run_data_path, delimiter="\t", header=None)
    run_data_list = []
    for _index, row in run_data.iterrows():
        run_data_list.append(
            TrackInputOutputDirectories(phenopacket_dir=Path(row[0]), results_dir=Path(row[1]))
        )
    return run_data_list
