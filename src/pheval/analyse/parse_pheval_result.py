from pathlib import Path

import pandas as pd

from pheval.post_processing.post_processing import PhEvalResult


def read_standardised_result(standardised_result_path: Path) -> dict:
    """Read the standardised result output and return a list of dictionaries."""
    try:
        return pd.read_csv(standardised_result_path, delimiter="\t").to_dict("records")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def parse_pheval_result(data_class_type: PhEvalResult, pheval_result: [dict]) -> [PhEvalResult]:
    """Parse PhEval result into specified dataclass type."""
    return [data_class_type(**row) for row in pheval_result]
