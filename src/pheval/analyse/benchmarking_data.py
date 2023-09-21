from dataclasses import dataclass
from pathlib import Path

from pheval.analyse.rank_stats import RankStats


@dataclass
class AnalysisResults:
    """Analysis results for a run."""

    results_dir: Path
    ranks: dict
    rank_stats: RankStats


@dataclass
class TrackRunPrioritisation:
    """Track prioritisation for a run."""

    gene_prioritisation: AnalysisResults = None
    variant_prioritisation: AnalysisResults = None
    disease_prioritisation: AnalysisResults = None

    def return_gene(self) -> AnalysisResults:
        """Return gene prioritisation analysis results."""
        return self.gene_prioritisation

    def return_variant(self) -> AnalysisResults:
        """Return variant prioritisation analysis results."""
        return self.variant_prioritisation

    def return_disease(self) -> AnalysisResults:
        """Return disease prioritisation analysis results."""
        return self.disease_prioritisation
