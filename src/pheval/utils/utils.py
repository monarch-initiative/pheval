"""Contains all pheval utility methods"""

import json
import random
from datetime import datetime
from pathlib import Path
from typing import List

import pandas as pd
import requests

from pheval.utils.logger import get_logger

logger = get_logger()


RESOURCES_DIR = Path(__file__).parent.parent / "resources"
MONDO_URL = "https://data.monarchinitiative.org/mappings/latest/mondo.sssom.tsv"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

METADATA_PATH = RESOURCES_DIR / "metadata.json"


def rand(df: pd.DataFrame, min_num: int, max_num: int, scramble_factor: float) -> float:
    """
    Numeric scrambling
    Args:
        df (pd.DataFrame): dataframe records
        min_num (int): min value from this records
        max_num (int): max value from this records
        scramble_factor (float): scramble factor scalar
    Returns:
        float: randomized number
    """
    try:
        return df + (random.uniform(min_num, max_num) * scramble_factor)
    except TypeError as err:
        logger.error(df, exc_info=err)
        return df


def semsim_scramble(
    input: Path,
    output: Path,
    columns_to_be_scrambled: List[str],
    scramble_factor: float = 0.5,
) -> pd.DataFrame:
    """
    Scrambles semantic similarity profile with a magnitude between 0 and 1 (scramble_factor:
    0 means no scrambling and 1 means complete randomisation).
    It then randomises the above scores with a degree of the scramble_factor
    and returns a scrambles pandas dataframe.
        Args:
              input (Path):
              scramble_factor (float) scalar scramble factor
              columns_to_be_scrambled (List[str]):
              columns that will be scrambled in semsim file (e.g. jaccard_similarity).
              output (Path)
        Returns:
            pd.Dataframe: scrambled dataframe
    """
    semsim = pd.read_csv(input, sep="\t")
    dataframe = semsim_scramble_df(semsim, columns_to_be_scrambled, scramble_factor)
    dataframe.to_csv(output, sep="\t", index=False)


def semsim_scramble_df(
    dataframe: pd.DataFrame,
    columns_to_be_scrambled: List[str],
    scramble_factor: float,
) -> pd.DataFrame:
    """scramble_semsim_df
    Args:
        dataframe (pd.DataFrame): dataframe that contains semsim profile
        scramble_factor (float) scalar scramble factor
        columns_to_be_scrambled (List[str]):
    Returns:
        pd.Dataframe: scrambled dataframe
    """
    for col in columns_to_be_scrambled:
        min_num = dataframe[col].min()
        max_num = dataframe[col].max()
        dataframe[col] = dataframe[col].apply(rand, args=(min_num, max_num, scramble_factor))
    return dataframe


def _update_metadata(file_name: str) -> None:
    """
    Update metadata.json with the current UTC timestamp for a file.
    Args:
        file_name (str): The file name.
    """
    timestamp = datetime.utcnow().isoformat() + "Z"
    metadata = {}
    if METADATA_PATH.exists():
        metadata = json.loads(METADATA_PATH.read_text())
    metadata[file_name] = timestamp
    METADATA_PATH.write_text(json.dumps(metadata, indent=2))
    logger.info(f"Updated metadata for {file_name}")


def _download_file(url: str, target_path: Path) -> None:
    """
    Helper to download a file from a URL to a target path.
    Args:
        url (str): url to download.
        target_path (Path): Path to download to.
    """
    try:
        logger.info(f"Downloading: {url}")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        RESOURCES_DIR.mkdir(parents=True, exist_ok=True)
        target_path.write_bytes(response.content)
        logger.info(f"Saved to: {target_path}")
        _update_metadata(target_path.name)
    except requests.RequestException as e:
        logger.error(f"Failed to download {url}: {e}")
        raise


def download_mondo_mapping() -> None:
    """Download latest MONDO SSSOM mapping file."""
    _download_file(MONDO_URL, RESOURCES_DIR / "mondo.sssom.tsv")


def download_hgnc_data() -> None:
    """Download latest HGNC complete set file."""
    _download_file(HGNC_URL, RESOURCES_DIR / "hgnc_complete_set.txt")


def get_resource_timestamp(file_name: str) -> str | None:
    """
    Return the ISO timestamp when the resource file was last updated.
    Args:
        file_name (str): The file name.
    """
    if METADATA_PATH.exists():
        with open(METADATA_PATH, "r") as f:
            return json.load(f).get(file_name)
    return None
