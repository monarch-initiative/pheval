"""Contains all pheval utility methods"""
import logging
import random
from pathlib import Path
from typing import List

import pandas as pd

info_log = logging.getLogger("info")


def rand(self: pd.DataFrame, min_num: int, max_num: int, scramble_factor: float) -> float:
    """
    Numeric scrambling
    Args:
        self (pd.DataFrame): dataframe records
        min_num (int): min value from this records
        max_num (int): max value from this records
        scramble_factor (float): scramble factor scalar
    Returns:
        float: randomized number
    """
    try:
        return self + (random.uniform(min_num, max_num) * scramble_factor)
    except TypeError as err:
        info_log.error(self, exc_info=err)
        return self


def semsim_randomisation(
    semsim_left: Path,
    semsim_right: Path,
    output: Path,
    scramble_factor: float = 0.5,
    times: int = 1,
    columns_to_be_scrambled: List[str] = ["jaccard_similarity"],
) -> pd.DataFrame:
    """
    Scrambles semantic similarity profile with a magnitude between 0 and 1 (scramble_factor:
    0 means no scrambling and 1 means complete randomisation).
    It then randomises the above scores with a degree of the scramble_factor
    and returns a scrambles pandas dataframe.
        Args:
              semsim_left (Path):
              semsim_right (Path):
              scramble_factor (float) scalar scramble factor
              columns_to_be_scrambled (List[str], optional): [description].
              Defaults to ['jaccard_similarity'].
              output (Path)
        Returns:
            pd.Dataframe: scrambled dataframe
    """
    semsim_left = pd.read_csv(semsim_left, sep="\t")
    semsim_right = pd.read_csv(semsim_right, sep="\t")
    semsim_left = scramble_semsim_df(
        semsim_left, scramble_factor, times, columns_to_be_scrambled, "1"
    )
    semsim_right = scramble_semsim_df(
        semsim_right, scramble_factor, times, columns_to_be_scrambled, "2"
    )
    dataframe = semsim_left.merge(semsim_right, on=["subject_id", "object_id"])
    dataframe = dataframe.filter(regex="subject_id|object_id|^SEM.", axis=1)
    dataframe.to_csv(output, sep="\t", index=False)


def scramble_semsim_df(
    dataframe: pd.DataFrame,
    scramble_factor: float,
    times: int,
    columns_to_be_scrambled: List[str],
    prefix: str,
) -> pd.DataFrame:
    """scramble_semsim_df

    Args:
        dataframe (pd.DataFrame): [description]
        scramble_factor (float): [description]
        times (int): [description]
        columns_to_be_scrambled (List[str]): [description]
        prefix (str): [description]

    Returns:
        pd.DataFrame: [description]
    """
    dict_of_cols = {}
    for col in columns_to_be_scrambled:
        min_num = dataframe[col].min()
        max_num = dataframe[col].max()
        for time in range(1, times + 1, 1):
            dict_of_cols[f"SEM{prefix}.{time}"] = dataframe[col].apply(
                rand, args=(min_num, max_num, scramble_factor)
            )
        dataframe = pd.concat([dataframe, pd.DataFrame(dict_of_cols)], axis=1)
    return dataframe
