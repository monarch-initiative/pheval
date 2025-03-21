from enum import Enum
from functools import wraps
from typing import Callable

import polars as pl


class ResultSchema(Enum):
    """
    Enum for different result schema formats.
    Attributes:
        GENE_RESULT_SCHEMA (pl.Schema): Schema for gene-based results.
        VARIANT_RESULT_SCHEMA (pl.Schema): Schema for variant-based results.
        DISEASE_RESULT_SCHEMA (pl.Schema): Schema for disease-based results.
    """

    GENE_RESULT_SCHEMA = pl.Schema(
        {
            "gene_symbol": pl.String,
            "gene_identifier": pl.String,
            "score": pl.Float64,
            "grouping_id": pl.Utf8,
        }
    )
    VARIANT_RESULT_SCHEMA = pl.Schema(
        {
            "chrom": pl.String,
            "start": pl.Int64,
            "end": pl.Int64,
            "ref": pl.String,
            "alt": pl.String,
            "score": pl.Float64,
            "grouping_id": pl.Utf8,
        }
    )
    DISEASE_RESULT_SCHEMA = pl.Schema(
        {
            "disease_identifier": pl.String,
            "score": pl.Float64,
            "grouping_id": pl.Utf8,
        }
    )

    def validate(self, results: pl.DataFrame) -> bool:
        """
        Validate that a DataFrame follows the expected schema.
        Args:
            results (pl.DataFrame): The DataFrame to validate.
        Raises:
            ValueError: If a required column is missing or the grouping_id column contains a null value.
            TypeError: If a column exists but has an incorrect data type.
        Returns:
            bool: True if the DataFrame is valid according to the schema.
        """
        expected_schema = self.value

        if "grouping_id" in results.columns and results["grouping_id"].null_count() > 0:
            raise ValueError("'grouping_id' column should not contain null values if provided.")

        for col_name, expected_type in expected_schema.items():
            if col_name not in results.schema:
                if col_name == "grouping_id":
                    continue
                raise ValueError(f"Missing required column: {col_name}")

            if results.schema[col_name] != expected_type:
                raise TypeError(
                    f"Column '{col_name}' has type {results.schema[col_name]}, expected {expected_type}"
                )

        return True


def validate_dataframe(schema: ResultSchema) -> Callable:
    """
    Decorator to validate DataFrame input based on a ResultSchema.
    Args:
        schema (ResultSchema): The expected schema from the `ResultSchema` enum.
    Returns:
        Callable: A wrapped function that validates the DataFrame before execution.
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(results: pl.DataFrame, *args, **kwargs):
            schema.validate(results)
            return func(results, *args, **kwargs)

        return wrapper

    return decorator
