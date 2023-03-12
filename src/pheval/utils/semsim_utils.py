"""
Contains all pheval utility methods
"""
from pathlib import Path

import numpy
import pandas as pd
import plotly.express as px

import pheval.utils.file_utils as file_utils


def filter_non_0_score(data: pd.DataFrame, col: str) -> pd.DataFrame:
    """Removes rows that have value equal to 0 based on the given column passed by col parameter

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


def diff_semsim(
    semsim_left: pd.DataFrame, semsim_right: pd.DataFrame, score_column: str, absolute_diff: bool
) -> pd.DataFrame:
    """Calculates score difference between two semantic similarity profiles

    Args:
        semsim_left (pd.DataFrame): first semantic similarity dataframe
        semsim_right (pd.DataFrame): second semantic similarity dataframe
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        absolute_diff (bool, optional): Whether the difference is absolute (True) or percentage (False).
        Defaults to True.

    Returns:
        pd.DataFrame: A dataframe with terms and its scores differences
    """
    df = pd.merge(semsim_left, semsim_right, on=["subject_id", "object_id"], how="outer")
    if absolute_diff:
        df["diff"] = df[f"{score_column}_x"] - df[f"{score_column}_y"]
        return df[["subject_id", "object_id", "diff"]]
    df["diff"] = df.apply(
        lambda row: get_percentage_diff(row[f"{score_column}_x"], row[f"{score_column}_y"]), axis=1
    )
    return df[["subject_id", "object_id", f"{score_column}_x", f"{score_column}_y", "diff"]]


def percentage_diff(semsim_left: Path, semsim_right: Path, score_column: str, output: Path):
    """Compares two semantic similarity profiles

    Args:
        semsim_left (Path): File path of the first semantic similarity profile
        semsim_right (Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        output (Path): Output path for the difference tsv file
    """
    clean_df = semsim_analysis(semsim_left, semsim_right, score_column, absolute_diff=False)
    clean_df.sort_values(by="diff", ascending=False).to_csv(output, sep="\t", index=False)


def semsim_heatmap_plot(semsim_left: Path, semsim_right: Path, score_column: str):
    """Plots semantic similarity profiles heatmap

    Args:
        semsim_left (Path): File path of the first semantic similarity profile
        semsim_right (Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
    """
    clean_df = semsim_analysis(semsim_left, semsim_right, score_column)
    df = clean_df.pivot(index="subject_id", columns="object_id", values="diff")
    fig = px.imshow(df, text_auto=True)
    fig.show()


def semsim_analysis(
    semsim_left: Path, semsim_right: Path, score_column: str, absolute_diff=True
) -> pd.DataFrame:
    """semsim_analysis

    Args:
        semsim_left (Path): File path of the first semantic similarity profile
        semsim_right (Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        absolute_diff (bool, optional): Whether the difference is absolute (True) or percentage (False).
        Defaults to True.

    Returns:
        [pd.DataFrame]: DataFrame with the differences between two semantic similarity profiles
    """
    validate_semsim_file_comparison(semsim_left, semsim_right)
    cols = ["subject_id", "object_id", score_column]
    semsim_left = pd.read_csv(semsim_left, sep="\t")
    semsim_right = pd.read_csv(semsim_right, sep="\t")
    file_utils.ensure_columns_exists(
        cols=cols,
        err_message="must exist in semsim dataframes",
        dataframes=[semsim_left, semsim_right],
    )
    semsim_left = parse_semsim(semsim_left, cols)
    semsim_right = parse_semsim(semsim_right, cols)
    diff_df = diff_semsim(semsim_left, semsim_right, score_column, absolute_diff)
    return filter_non_0_score(diff_df, "diff")


def validate_semsim_file_comparison(semsim_left: Path, semsim_right: Path):
    """Checks if files exist and whether they're different
    Args:
        semsim_left (Path): File path of the first semantic similarity profile
        semsim_right (Path): File path of the second semantic similarity profile
    Raises:
        Exception: FileNotFoundException
    """
    if semsim_left == semsim_right:
        errmsg = "Semantic similarity profiles are equal. Make sure you have selected different files to analyze"
        raise Exception(errmsg)
    file_utils.ensure_file_exists(semsim_left, semsim_right)


def get_percentage_diff(current_number: float, previous_number: float) -> float:
    """Gets the percentage difference between two numbers

    Args:
        current_number (float): second number in comparison
        previous_number (float): first number in comparison

    Returns:
        float: percentage difference between two numbers
    """
    try:
        if current_number == previous_number:
            return "{:.2%}".format(0)
        if current_number > previous_number:
            number = (1 - ((current_number / previous_number))) * 100
        else:
            number = (100 - ((previous_number / current_number) * 100)) * -1
        return "{:.2%}".format(number / 100)
    except ZeroDivisionError:
        return None
