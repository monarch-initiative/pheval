import logging
from pathlib import Path
from typing import List

import pandas as pd

from pheval.post_processing.post_processing import PhEvalResult

info_log = logging.getLogger("info")


def read_standardised_result(standardised_result_path: Path) -> List[dict]:
    """
    Read the standardised result output and return a list of dictionaries.

    Args:
        standardised_result_path (Path): The path to the file containing the standardised result output.

    Returns:
        List[dict]: A list of dictionaries representing the content of the standardised result file.
    """
    if standardised_result_path.is_file():
        return pd.read_csv(standardised_result_path, delimiter="\t").to_dict("records")
    else:
        info_log.info(f"Could not find {standardised_result_path}")
        return pd.DataFrame().to_dict("records")


def parse_pheval_result(
    data_class_type: PhEvalResult, pheval_result: List[dict]
) -> List[PhEvalResult]:
    """
    Parse PhEval result into specified dataclass type.

    Args:
        data_class_type (PhEvalResult): The data class type to parse the result into.
        pheval_result (List[dict]): A list of dictionaries representing the PhEval result.

    Returns:
        List[PhEvalResult]: A list of instances of the specified data class type,
        each instance representing a row in the PhEval result.
    """
    return [data_class_type(**row) for row in pheval_result]
