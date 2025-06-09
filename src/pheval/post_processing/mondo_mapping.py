from pathlib import Path

import polars as pl


def parse_mondo_mapping_table() -> pl.DataFrame:
    """
    Parse the Mondo SSSOM table.
    Returns:
        pl.DataFrame: Mondo SSSOM table.
    """
    return pl.read_csv(
        Path(__file__).parent.parent / "resources" / "mondo.sssom.tsv",
        separator="\t",
        comment_prefix="#",
    )


def map_disease_id(disease_identifier: str, mondo_mapping_table: pl.DataFrame) -> str:
    """
    Map a disease identifier to MONDO ID using the Mondo SSSOM mapping.
    Args:
        disease_identifier (str): The disease identifier to map to MONDO.
        mondo_mapping_table (pl.DataFrame): The Mondo SSSOM table.
    Returns:
        str: The MONDO ID.
    """
    mapped_identifier = mondo_mapping_table.filter(pl.col("object_id") == disease_identifier)
    if mapped_identifier.height > 0:
        return mapped_identifier["subject_id"].item()
    return disease_identifier
