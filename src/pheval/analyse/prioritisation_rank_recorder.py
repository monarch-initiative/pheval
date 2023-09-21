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
    """Compare the ranks of different runs."""

    index: int
    directory: Path
    prioritisation_result: Union[
        GenePrioritisationResult, VariantPrioritisationResult, DiseasePrioritisationResult
    ]
    run_comparison: defaultdict

    def _record_gene_rank(self) -> None:
        """Record gene prioritisation rank."""
        self.run_comparison[self.index]["Gene"] = self.prioritisation_result.gene

    def _record_variant_rank(self) -> None:
        """Record variant prioritisation rank."""
        variant = self.prioritisation_result.variant
        self.run_comparison[self.index]["Variant"] = "_".join(
            [variant.chrom, str(variant.pos), variant.ref, variant.alt]
        )

    def _record_disease_rank(self) -> None:
        """Record disease prioritisation rank."""
        self.run_comparison[self.index][
            "Disease"
        ] = self.prioritisation_result.disease.disease_identifier

    def record_rank(self) -> None:
        """Records the rank for different runs."""
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
