import itertools
import re
import unicodedata
from os import path
from pathlib import Path
from typing import List

import pandas as pd
import yaml
from serde import to_dict

from pheval.run_metadata import BasicOutputRunMetaData


def files_with_suffix(directory: Path, suffix: str) -> list[Path]:
    """
    Obtains all files ending in a specified suffix from a given directory.

    Args:
        directory (Path): The directory path.
        suffix (str): The specified suffix to filter files.

    Returns:
        list[Path]: A list of Path objects representing files with the specified suffix.
    """
    files = [file_path for file_path in directory.iterdir() if file_path.suffix == suffix]
    files.sort()
    return files


def all_files(directory: Path) -> list[Path]:
    """
    Obtains all files from a given directory.

    Args:
        directory (Path): The directory path.

    Returns:
        list[Path]: A list of Path objects representing all files in the directory.
    """
    files = [file_path for file_path in directory.iterdir()]
    files.sort()
    return files


def is_gzipped(file_path: Path) -> bool:
    """
    Confirms whether a file is gzipped.

    Args:
        file_path (Path): The path to the file.

    Returns:
        bool: True if the file is gzipped, False otherwise.
    """
    return file_path.name.endswith(".gz")


def normalise_file_name(file_path: Path) -> str:
    """
    Normalises the file name by removing diacritical marks (accents) from Unicode characters.

    Args:
        file_path (Path): The path to the file.

    Returns:
        str: The normalised file name without diacritical marks.
    """
    normalised_file_name = unicodedata.normalize("NFD", str(file_path))
    return re.sub("[\u0300-\u036f]", "", normalised_file_name)


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
    """
    Write the metadata for a run to a YAML file.

    Args:
        output_dir (Path): The directory where the metadata file will be saved.
        meta_data (BasicOutputRunMetaData): The metadata to be written.
    """
    with open(Path(output_dir).joinpath("results.yml"), "w") as metadata_file:
        yaml.dump(to_dict(meta_data), metadata_file, sort_keys=False, default_style="")
    metadata_file.close()
