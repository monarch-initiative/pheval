"""
Contains all pheval utility methods
"""
import itertools
from os import path
from typing import List

import pandas as pd
from click import Path


def ensure_file_exists(*files: str):
    """ensure_file_exists"""
    for file in files:
        if not path.isfile(file):
            raise FileNotFoundError(f"File {file} not found")


def ensure_columns_exists(**kwargs):
    """ensure_columns_exists"""
    flat_cols = list(itertools.chain(*kwargs.get("cols")))
    dataframes = kwargs.get("dataframes")
    if not dataframes or not flat_cols:
        return
    err_msg = f"""columns: {", ".join(flat_cols[:-1])} and {flat_cols[-1]} \
- must be present in both left and right files"""
    for dataframe in kwargs.get("dataframes", [pd.DataFrame()]):
        if not all(x in dataframe.columns for x in flat_cols):
            raise ValueError(err_msg)


def intersect_scores(
    left: Path,
    right: Path,
    join_columns: List[str],
    comparison_columns: List[str],
    sort_columns: List[str],
):
    """intersect_scores"""
    left_df = pd.read_csv(left, sep="\t")
    right_df = pd.read_csv(right, sep="\t")
    ensure_columns_exists(
        cols=[join_columns, comparison_columns, sort_columns], dataframes=[left_df, right_df]
    )
    merged = left_df.merge(right_df, on=join_columns, suffixes=("_left", "_right"))
    return merged


def calc_semsim(
    left: Path,
    right: Path,
    join_columns: List[str],
    comparison_columns: List[str],
    sort_columns: List[str],
):
    """calc_semsim"""
    inter_scores = intersect_scores(left, right, join_columns, comparison_columns, sort_columns)

    for col in comparison_columns:
        inter_scores[f"{col}_diff"] = inter_scores[f"{col}_left"] - inter_scores[f"{col}_right"]
    sort_columns_suffix = list(
        itertools.chain(*[[f"{j}_{i}" for j in sort_columns] for i in ["left", "right"]])
    )
    ascending = list(map(lambda sort_col: sort_col[0] != "-", sort_columns_suffix))

    inter_scores.sort_values(by=sort_columns_suffix, ascending=ascending, inplace=True, axis=0)
    flat_cols = list(
        itertools.chain(
            *[
                [f"{j}_{i}" for j in itertools.chain(*[comparison_columns, sort_columns])]
                for i in ["left", "right"]
            ],
            inter_scores.columns[inter_scores.columns.str.endswith("_diff")],
        )
    )
    flat_sort = sorted(
        flat_cols,
        key=lambda x: x.split("_")[0],
    )
    return inter_scores.loc[:, list(itertools.chain(*[join_columns, flat_sort]))]
