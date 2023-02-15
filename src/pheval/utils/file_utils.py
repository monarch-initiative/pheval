import itertools
import re
import unicodedata
from os import path
from pathlib import Path
from typing import List

import pandas as pd


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


def normalise_file_name(file_path: Path) -> str:
    normalised_file_name = unicodedata.normalize("NFD", str(file_path))
    return re.sub("[\u0300-\u036f]", "", normalised_file_name)


def obtain_closest_file_name(file_to_be_queried: Path, file_paths: list[Path]) -> Path:
    """Obtains the closest file name when given a template file name and a list of full path of files to be queried."""
    # closest_file_match = difflib.get_close_matches(
    #     Path(file_to_be_queried).stem,
    #     [Path(file_path).stem for file_path in file_paths],
    #     cutoff=0.4,
    # )[0]
    # return [file_path for file_path in file_paths if closest_file_match == file_path.stem][0]
    return [
        file_path
        for file_path in file_paths
        if normalise_file_name(file_to_be_queried.stem).startswith(
            normalise_file_name(file_path.stem)
        )
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
