from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from pheval.analyse.prioritisation_result_types import (
    DiseasePrioritisationResult,
    GenePrioritisationResult,
    VariantPrioritisationResult,
)


@dataclass
class PrioritisationRankRecorder:
    """
    Record ranks for different types of prioritisation results.

    Attributes:
        index (int): The index representing the run.
        directory (Path): The result directory path.
        prioritisation_result (Union[GenePrioritisationResult, VariantPrioritisationResult,
            DiseasePrioritisationResult]): The prioritisation result object.
        run_comparison (defaultdict): The comparison dictionary to record ranks.
    """

    index: int
    directory: Path
    prioritisation_result: Union[
        GenePrioritisationResult, VariantPrioritisationResult, DiseasePrioritisationResult
    ]
    run_comparison: defaultdict

    def _record_gene_rank(self) -> None:
        """
        Record gene prioritisation rank.

        This method updates the 'Gene' key in the run comparison dictionary with the gene
        information extracted from the correct prioritisation result.
        """
        self.run_comparison[self.index]["Gene"] = self.prioritisation_result.gene

    def _record_variant_rank(self) -> None:
        """
        Record variant prioritisation rank.

        This method updates the 'Variant' key in the run comparison dictionary with the variant
        information extracted from the correct prioritisation result.
        """
        variant = self.prioritisation_result.variant
        self.run_comparison[self.index]["Variant"] = "-".join(
            [variant.chrom, str(variant.pos), variant.ref, variant.alt]
        )

    def _record_disease_rank(self) -> None:
        """
        Record disease prioritisation rank.

        This method updates the 'Disease' key in the run comparison dictionary with the disease
        information extracted from the correct prioritisation result.
        """
        self.run_comparison[self.index][
            "Disease"
        ] = self.prioritisation_result.disease.disease_identifier

    def record_rank(self) -> None:
        """
        Record the prioritisation ranks for different runs.

        It assigns the prioritisation rank and associated details such as phenopacket name
        and prioritisation result type ('Gene', 'Variant', or 'Disease') to the run comparison
        dictionary for each respective run, allowing comparison and analysis of the ranks of correct results
        across different runs.
        """
        self.run_comparison[self.index][
            "Phenopacket"
        ] = self.prioritisation_result.phenopacket_path.name
        if type(self.prioritisation_result) is GenePrioritisationResult:
            self._record_gene_rank()
        elif type(self.prioritisation_result) is VariantPrioritisationResult:
            self._record_variant_rank()
        elif type(self.prioritisation_result) is DiseasePrioritisationResult:
            self._record_disease_rank()
        self.run_comparison[self.index][self.directory] = self.prioritisation_result.rank
