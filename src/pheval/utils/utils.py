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
        return self + (random.uniform(min_num, max_num) * scramble_factor)  # noqa
    except TypeError as err:
        info_log.error(self, exc_info=err)
        return self


def semsim_randomisation(
    input: Path,
    output: Path,
    scramble_factor: float = 0.5,
    columns_to_be_scrambled: List[str] = tuple(["jaccard_similarity"]),
) -> pd.DataFrame:
    """
    Scrambles semantic similarity profile with a magnitude between 0 and 1 (scramble_factor:
    0 means no scrambling and 1 means complete randomisation).
    It then randomises the above scores with a degree of the scramble_factor
    and returns a scrambles pandas dataframe.
        Args:
              input (Path):
              scramble_factor (float) scalar scramble factor
              columns_to_be_scrambled (List[str], optional): [description].
              Defaults to ['jaccard_similarity'].
              output (Path)
        Returns:
            pd.Dataframe: scrambled dataframe
    """
    semsim = pd.read_csv(input, sep="\t")
    dataframe = scramble_semsim_df(semsim, scramble_factor, columns_to_be_scrambled)
    dataframe.to_csv(output, sep="\t", index=False)


def scramble_semsim_df(
    dataframe: pd.DataFrame,
    scramble_factor: float,
    columns_to_be_scrambled: List[str],
) -> pd.DataFrame:
    """scramble_semsim_df
    Args:
        dataframe (pd.DataFrame): [description]
        scramble_factor (float): [description]
        columns_to_be_scrambled (List[str]): [description]
    Returns:
        pd.DataFrame: [description]
    """
    for col in columns_to_be_scrambled:
        min_num = dataframe[col].min()
        max_num = dataframe[col].max()
        dataframe[col] = dataframe[col].apply(rand, args=(min_num, max_num, scramble_factor))
    return dataframe


def semsimconvert(input: Path, output: Path, subject_prefix: str, object_prefix: str):
    """convert semsim profile to an exomiser specie mapping schema
    Example of a mapping schema:

    MAPPING_ID;HP_ID;HP_TERM;MP_ID;MP_TERM;SIMJ;IC;SCORE;LCS_ID;LCS_TERM

    To generate the example above, the parameter subject_prefix needs to be "HP" and object_prefix needs to be "MP"

    Args:
        input (Path): Path to the semsim profile
        output (Path): Path where sql file will be written
        subject_prefix: (str): mapping subject prefix (e.g. HP)
        object_prefx (str): mapping object prefix (e.g. MP)
    """
    input_data = pd.read_csv(input, sep="\t")
    input_data = input_data.rename(
        columns={
            "subject_id": f"{subject_prefix}_ID",
            "subject_label": f"{subject_prefix}_TERM",
            "object_id": f"{object_prefix}_ID",
            "object_label": f"{object_prefix}_TERM",
            "jaccard_similarity": "SIMJ",
            "ancestor_information_content": "IC",
            "phenodigm_score": "SCORE",
            "ancestor_id": "LCS_ID",
            "ancestor_label": "LCS_TERM",
        }
    )[
        [
            f"{subject_prefix}_ID",
            f"{subject_prefix}_TERM",
            f"{object_prefix}_ID",
            f"{object_prefix}_TERM",
            "SIMJ",
            "IC",
            "SCORE",
            "LCS_ID",
            "LCS_TERM",
        ]
    ]
    input_data.insert(0, "MAPPING_ID", None)
    input_data["SCORE"].replace("None", "NULL", inplace=True)
    input_data.to_csv(output, index=False, sep="\t")
    semsim2h2(input_data, f"{output}", subject_prefix, object_prefix)


def semsim2h2(input: Path, output: Path, subject_prefix: str, object_prefix: str):
    """This function converts the exomiser mapping table to sql format.
    Args:
        input (Path): Path to the semsim profile
        output (Path): Path where sql file will be written
        subject_prefix: (str): mapping subject prefix (e.g. HP)
        object_prefx (str): mapping object prefix (e.g. MP)
    """
    with open(output, "w") as file:
        insert = f"""
TRUNCATE TABLE EXOMISER.{subject_prefix}_{object_prefix}_MAPPINGS;
INSERT INTO EXOMISER.{subject_prefix}_{object_prefix}_MAPPINGS
(MAPPING_ID, {subject_prefix}_ID, {subject_prefix}_TERM, {object_prefix}_ID, {object_prefix}_TERM, SIMJ, IC, SCORE, LCS_ID, LCS_TERM)
VALUES
"""  # noqa
        file.write(insert)
        for idx, data in input.iterrows():
            values = f"""
({idx}, '{data[f'{subject_prefix}_ID']}', '{data[f'{subject_prefix}_TERM']}', '{data[f'{object_prefix}_ID']}', '{data[f'{object_prefix}_TERM']}', {data['SIMJ']}, {data['IC']}, {data['SCORE']}, '{data['LCS_ID']}', '{data['LCS_TERM']}')"""  # noqa
            if idx > (len(input) < 1):
                file.write(",")
            file.write(values)
        file.write(";")
