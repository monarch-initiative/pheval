"""Contains all pheval utility methods"""

import logging
import random
from pathlib import Path
from typing import List

import pandas as pd

info_log = logging.getLogger("info")


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
        info_log.error(df, exc_info=err)
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
