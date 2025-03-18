from enum import Enum
from typing import List, NamedTuple


class BenchmarkOutputType(NamedTuple):
    """
    Represents the structure of benchmark output types.

    Attributes:
        prioritisation_type_string (str): The type of prioritisation being performed.
        y_label (str): The label for the y-axis in performance evaluation plots.
        columns (List[str]): The list of column names relevant to the benchmark output.
        result_directory (str): The directory where benchmark results are stored.
    """

    prioritisation_type_string: str
    y_label: str
    columns: List[str]
    result_directory: str


class BenchmarkOutputTypeEnum(Enum):
    """
    Enumeration of benchmark output types, representing different entities.

    Attributes:
        GENE (BenchmarkOutputType): Benchmark output type for gene prioritisation.
        VARIANT (BenchmarkOutputType): Benchmark output type for variant prioritisation.
        DISEASE (BenchmarkOutputType): Benchmark output type for disease prioritisation.
    """

    GENE = BenchmarkOutputType(
        "gene",
        "Disease-causing genes (%)",
        ["gene_identifier", "gene_symbol"],
        "pheval_gene_results",
    )
    VARIANT = BenchmarkOutputType(
        "variant", "Disease-causing variants (%)", ["variant_id"], "pheval_variant_results"
    )
    DISEASE = BenchmarkOutputType(
        "disease", "Known diseases (%)", ["disease_identifier"], "pheval_disease_results"
    )
