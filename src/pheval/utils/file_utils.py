import difflib
import itertools
from os import path
from pathlib import Path
from typing import List

import pandas as pd
import yaml
from serde import to_dict

from pheval.run_metadata import BasicOutputRunMetaData


def files_with_suffix(directory: Path, suffix: str):
    """Obtains all files ending in a specified suffix from a given directory."""
    files = [path for path in directory.iterdir() if path.suffix == suffix]
    files.sort()
    return files


def all_files(directory: Path) -> list[Path]:
    """Obtains all files from a given directory."""
    files = [path for path in directory.iterdir()]
    files.sort()
    return files


def is_gzipped(path: Path) -> bool:
    """Confirms whether a file is gzipped."""
    return path.name.endswith(".gz")


def obtain_closest_file_name(file_to_be_queried: Path, file_paths: list[Path]) -> Path:
    """Obtains the closest file name when given a template file name and a list of full path of files to be queried."""
    closest_file_match = Path(
        str(
            difflib.get_close_matches(
                str(file_to_be_queried.name),
                [str(file_path.name) for file_path in file_paths],
            )[0]
        )
    )
    return [
        file_path for file_path in file_paths if Path(closest_file_match) == Path(file_path.name)
    ][0]


def ensure_file_exists(*files: str):
    """Ensures the existence of files passed as parameter
    Raises:
        FileNotFoundError: If any file passed as a parameter doesn't exist a FileNotFound Exception will be raised
    """
    for file in files:
        if not path.isfile(file):
            raise FileNotFoundError(f"File {file} not found")


def ensure_columns_exists(cols: list, dataframes: List[pd.DataFrame], err_message: str = ""):
    """Ensures the columns exist in dataframes passed as argument (e.g)

    "
    ensure_columns_exists(
        cols=['column_a', 'column_b, 'column_c'],
        err_message="Custom error message if any column doesn't exist in any dataframe passed as argument",
        dataframes=[data_frame1, data_frame2],
    )
    "

    """
    flat_cols = list(itertools.chain(cols))
    if not dataframes or not flat_cols:
        return
    if err_message:
        err_msg = f"""columns: {", ".join(flat_cols[:-1])} and {flat_cols[-1]} {err_message}"""
    else:
        err_msg = f"""columns: {", ".join(flat_cols[:-1])} and {flat_cols[-1]} \
- must be present in both left and right files"""
    for dataframe in dataframes:
        if not all(x in dataframe.columns for x in flat_cols):
            raise ValueError(err_msg)


def write_metadata(output_dir: Path, meta_data: BasicOutputRunMetaData) -> None:
    """Write the metadata for a run."""
    with open(Path(output_dir).joinpath("results.yml"), "w") as metadata_file:
        yaml.dump(to_dict(meta_data), metadata_file, sort_keys=False, default_style="")
    metadata_file.close()

