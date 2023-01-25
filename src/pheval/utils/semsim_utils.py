"""
Contains all pheval utility methods
"""


import click
import numpy
import pandas as pd
import plotly.express as px

import pheval.utils.file_utils as file_utils


def filter_non_0_score(data: pd.DataFrame, col: str) -> pd.DataFrame:
    """Gets rid of data with a score of 0 in the selected column

    Args:
        data (pd.DataFrame): Dirty dataframe
        col (str): Column to be filtered

    Returns:
        pd.DataFrame: Filtered dataframe
    """
    return data[data[col] != 0]


def parse_semsim(df: pd.DataFrame, cols: list) -> pd.DataFrame:
    """Parses semantic similarity profiles converting the score column as a numeric value and dropping the null ones

    Args:
        df (pd.DataFrame): semantic similarity profile dataframe
        cols (list): list of columns that will be selected on semsim data

    Returns:
        pd.Dataframe: parsed semantic similarity dataframe
    """
    df[cols[-1]] = pd.to_numeric(df[cols[-1]], errors="coerce")
    df.replace("None", numpy.nan).dropna(subset=cols[-1], inplace=True)
    return df


def diff_semsim(data1: pd.DataFrame, data2: pd.DataFrame, score_column: str) -> pd.DataFrame:
    """Calculates score difference between two semantic similarity profiles

    Args:
        data1 (pd.DataFrame): first semantic similarity dataframe
        data2 (pd.DataFrame): second semantic similarity dataframe
        score_column (str): [description]

    Returns:
        pd.DataFrame: A dataframe with terms and its scores differences
    """
    keys = ["subject_id", "object_id"]
    df = pd.merge(data1, data2, on=keys, how="outer")

    df["diff"] = df[f"{score_column}_x"] - df[f"{score_column}_y"]
    return df[["subject_id", "object_id", "diff"]]


def semsim_heatmap_plot(file1: click.Path, file2: click.Path, score_column: str):
    """Plots semantic similarity profiles heatmap

    Args:
        file1 (click.Path): File path of the first semantic similarity profile
        file2 (click.Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
    """
    if file1 == file2:
        errmsg = "Semantic similarity profiles are equal. Make sure you have selected different files to analyze"
        raise Exception(errmsg)
    file_utils.ensure_file_exists(file1, file2)
    cols = ["subject_id", "object_id", score_column]
    data1 = pd.read_csv(file1, sep="\t")
    data2 = pd.read_csv(file2, sep="\t")
    file_utils.ensure_columns_exists(
        cols=[cols],
        message="must exist in semsim dataframes",
        dataframes=[data1, data2],
    )
    data1 = parse_semsim(data1, cols)
    data2 = parse_semsim(data2, cols)
    diff_df = diff_semsim(data1, data2, score_column)
    clean_df = filter_non_0_score(diff_df, "diff")
    df = clean_df.pivot(index="subject_id", columns="object_id", values="diff")
    fig = px.imshow(df, text_auto=True)
    fig.show()
