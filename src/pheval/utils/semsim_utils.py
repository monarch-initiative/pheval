"""
Contains all pheval utility methods
"""
import logging
import subprocess
from itertools import combinations
from pathlib import Path
from typing import List

import pandas as pd
import plotly.express as px
import polars as pl
import seaborn as sns
from matplotlib import pyplot as plt

from tqdm import tqdm

import pheval.utils.file_utils as file_utils

info_log = logging.getLogger("info")


def semsim_comparison(input: List[Path], score_column: str, analysis: str, output: Path):
    """Makes a paired semantic similarity profiles comparison based on a chosen score column

    Args:
        input (List[Path]): semsim profiles path's
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        analysis (str): There are three types of analysis:
        heatmap - Generates a heatmap plot that shows the differences between the semantic similarity profiles using the
        score column for this purpose. Defaults to "heatmap".
        percentage_diff - Calculates the score column percentage difference between the semantic similarity profiles.
        distribution - Plot showing the semsim score's distributions

    """
    for s in set(combinations(input, 2)):
        semsim_left = s[0]
        semsim_right = s[1]
        if analysis == "heatmap":
            semsim_heatmap_plot(semsim_left, semsim_right, score_column)
        if analysis == "percentage_diff":
            percentage_diff(semsim_left, semsim_right, score_column, output)
    if analysis == "distribution":
        distribution(input, score_column, output)


def filter_non_0_score(data: pl.DataFrame, col: str) -> pd.DataFrame:
    """Removes rows that have value equal to 0 based on the given column passed by col parameter

    Args:
        data (pl.DataFrame): Dirty dataframe
        col (str): Column to be filtered

    Returns:
        pl.DataFrame: Filtered dataframe
    """
    return data.filter(pl.col(col) != 0)


def parse_semsim(df: pl.DataFrame, cols: list) -> pd.DataFrame:
    """Parses semantic similarity profiles converting the score column as a numeric value and dropping the null ones

    Args:
        df (pl.DataFrame): semantic similarity profile dataframe
        cols (list): list of columns that will be selected on semsim data

    Returns:
        pd.Dataframe: parsed semantic similarity dataframe
    """
    df.with_columns(pl.col(cols[-1]).cast(pl.Float64))
    df[cols[-1]].set(df[cols[-1]].is_null(), None)

    df.drop_nulls(cols[-1])
    return df


def diff_semsim(
    semsim_left: pl.DataFrame, semsim_right: pl.DataFrame, score_column: str, absolute_diff: bool
) -> pl.DataFrame:
    """Calculates score difference between two semantic similarity profiles

    Args:
        semsim_left (pl.DataFrame): first semantic similarity dataframe
        semsim_right (pl.DataFrame): second semantic similarity dataframe
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        absolute_diff (bool, optional): Whether the difference is absolute (True) or percentage (False).
        Defaults to True.

    Returns:
        pl.DataFrame: A dataframe with terms and its scores differences
    """
    df = semsim_left.join(semsim_right, on=["subject_id", "object_id"], how="outer")
    df.drop_nulls(score_column)
    df.drop_nulls(f"{score_column}_right")
    if absolute_diff:
        df = df.with_columns((pl.col(score_column) - pl.col(f"{score_column}_right")).alias("diff"))
        return df[["subject_id", "object_id", f"{score_column}", f"{score_column}_right", "diff"]]
    df = df.with_columns(
        # horizontal sum with a custom apply
        pl.struct([score_column, f"{score_column}_right"])
        .apply(lambda x: get_percentage_diff(x[score_column], x[f"{score_column}_right"]))
        .alias("diff")
    )
    return df[["subject_id", "object_id", f"{score_column}", f"{score_column}_right", "diff"]]


def distribution(input: List[Path], score_column: str, output: Path):
    df_list = []
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.figsize"] = [20, 3.50 * len(input)]
    _, axes = plt.subplots(len(input), 1)
    for idx, i in enumerate(input):
        print(f"Reading {Path(i).stem}")
        df = pl.read_csv(i, separator="\t")
        df = df[["subject_id", "object_id", f"{score_column}"]]
        df = df.with_columns(semsim=pl.lit(Path(i).stem))
        df_list.append(df)
        axes[idx].ticklabel_format(style="plain", axis="both")
        axes[idx].set_xlabel(score_column)
        sns.histplot(df[score_column], bins=20, ax=axes[idx]).set_title(Path(i).stem)
    plt.setp(axes, ylim=axes[0].get_ylim())
    print("Concatenating data")
    df_concat = pl.concat(df_list)
    print(f"Saving plot in {output}/bars.png")
    plt.savefig(f"{output}/bars.png")
    plt.clf()
    sns.histplot(
        df_concat,
        x=score_column,
        bins=10,
        multiple="dodge",
        fill=True,
        kde=True,
        alpha=0.5,
        hue="semsim",
    )
    print(f"Saving plot in {output}/dist.png")
    plt.savefig(f"{output}/dist.png")


def percentage_diff(semsim_left: Path, semsim_right: Path, score_column: str, output: Path):
    """Compares two semantic similarity profiles

    Args:
        semsim_left (Path): File path of the first semantic similarity profile
        semsim_right (Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        output (Path): Output path for the difference tsv file
    """
    fname_left = Path(semsim_left).stem
    fname_right = Path(semsim_right).stem
    fname = f"{output}/{fname_left}-{fname_right}.diff.tsv"
    Path(fname).unlink(missing_ok=True)
    for idx, clean_df in enumerate(
        semsim_analysis(semsim_left, semsim_right, score_column, absolute_diff=False)
    ):
        (
            clean_df.drop_nulls("diff")
            .sort("diff", descending=True)
            .rename(
                {
                    score_column: f"{fname_left}_{score_column}",
                    f"{score_column}_right": f"{fname_right}_{score_column}",
                }
            )
            .write_csv(fname, has_header=idx == 0, separator="\t")
        )


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
    fig.update_layout(
        title=f"{Path(semsim_left).stem} - {Path(semsim_right).stem}", xaxis_nticks=36
    )
    fig.show()


def semsim_analysis(semsim_left: Path, semsim_right: Path, score_column: str, absolute_diff=True):
    """semsim_analysis

    Args:
        semsim_left (Path): File path of the first semantic similarity profile
        semsim_right (Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        absolute_diff (bool, optional): Whether the difference is absolute (True) or percentage (False).
        Defaults to True.

    Yields:
        pd.DataFrame: DataFrame with the differences between two semantic similarity profiles
    """
    validate_semsim_file_comparison(semsim_left, semsim_right)
    cols = ["subject_id", "object_id", score_column]
    batch_size = 100000
    count = int(subprocess.check_output(["wc", "-l", semsim_left]).split()[0])
    reader_left = pl.read_csv_batched(semsim_left, separator="\t", batch_size=batch_size)
    reader_right = pl.read_csv_batched(semsim_right, separator="\t", batch_size=batch_size)
    # file_utils.ensure_columns_exists(
    #     cols=cols,
    #     err_message="must exist in semsim dataframes",
    #     dataframes=[reader_left, reader_right],
    # )
    batches_left = reader_left.next_batches(5)
    batches_right = reader_right.next_batches(5)
    with tqdm(total=count - 1) as bar:
        while batches_left or batches_right:
            for input_data in zip(batches_left, batches_right):
                semsim_left = parse_semsim(input_data[0], cols)
                semsim_right = parse_semsim(input_data[1], cols)
                diff_df = diff_semsim(semsim_left, semsim_right, score_column, absolute_diff)
                bar.update(input_data[0].shape[0])
                if not absolute_diff:
                    yield diff_df
                else:
                    yield filter_non_0_score(diff_df, "diff")
            batches_left = reader_left.next_batches(5)
            batches_right = reader_right.next_batches(5)


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
        if not current_number or not previous_number:
            return None
        if current_number == previous_number:
            return "{:.2%}".format(0)
        if current_number > previous_number:
            number = (1 - ((current_number / previous_number))) * 100
        else:
            number = (100 - ((previous_number / current_number) * 100)) * -1
        return "{:.2%}".format(number / 100)
    except ZeroDivisionError:
        return None
