from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import click
import pandas as pd

from pheval.analyse.generate_plots import AnalysisResults, TrackRunPrioritisation
from pheval.analyse.generate_summary_outputs import (
    RankStatsWriter,
    generate_benchmark_comparison_disease_output,
    generate_benchmark_comparison_gene_output,
    generate_benchmark_comparison_variant_output,
    generate_benchmark_disease_output,
    generate_benchmark_gene_output,
    generate_benchmark_variant_output,
)
from pheval.analyse.parse_pheval_result import parse_pheval_result, read_standardised_result
from pheval.analyse.rank_stats import RankStats
from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)
from pheval.prepare.custom_exceptions import InputError
from pheval.utils.file_utils import all_files, files_with_suffix, obtain_closest_file_name
from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    PhenopacketUtil,
    ProbandCausativeGene,
    ProbandDisease,
    phenopacket_reader,
)


@dataclass
class GenePrioritisationResult:
    """Store rank data for causative genes."""

    phenopacket_path: Path
    gene: str
    rank: int = 0


@dataclass
class VariantPrioritisationResult:
    """Store rank data for causative variants."""

    phenopacket_path: Path
    variant: GenomicVariant
    rank: int = 0


@dataclass
class DiseasePrioritisationResult:
    """Store rank data for known diseases."""

    phenopacket_path: Path
    disease: ProbandDisease
    rank: int = 0


@dataclass
class PrioritisationRankRecorder:
    """Compare the ranks of different runs."""

    index: int
    directory: Path
    prioritisation_result: VariantPrioritisationResult or GenePrioritisationResult
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


@dataclass
class TrackInputOutputDirectories:
    """Track the input testdata for a corresponding pheval output directory"""

    phenopacket_dir: Path
    results_dir: Path


def _parse_run_data_text_file(run_data_path: Path) -> [TrackInputOutputDirectories]:
    """Parse run data .txt file returning a list of input testdata and corresponding output directories."""
    run_data = pd.read_csv(run_data_path, delimiter="\t", header=None)
    run_data_list = []
    for _index, row in run_data.iterrows():
        run_data_list.append(
            TrackInputOutputDirectories(phenopacket_dir=Path(row[0]), results_dir=Path(row[1]))
        )
    return run_data_list


class AssessGenePrioritisation:
    """Assess gene prioritisation."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_gene_results: [RankedPhEvalGeneResult],
        threshold: float,
        score_order: str,
        proband_causative_genes: [ProbandCausativeGene],
    ):
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_gene_results = standardised_gene_results
        self.threshold = threshold
        self.score_order = score_order
        self.proband_causative_genes = proband_causative_genes

    def _record_gene_prioritisation_match(
        self,
        gene: ProbandCausativeGene,
        result_entry: RankedPhEvalGeneResult,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if found within results."""
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol, rank)

    def _assess_gene_with_threshold_ascending_order(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _assess_gene_with_threshold(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _record_matched_gene(
        self, gene: ProbandCausativeGene, rank_stats: RankStats, standardised_gene_result: pd.Series
    ) -> GenePrioritisationResult:
        """Return the gene rank result - dealing with the specification of a threshold."""
        if float(self.threshold) == 0.0:
            return self._record_gene_prioritisation_match(
                gene, standardised_gene_result, rank_stats
            )
        else:
            return (
                self._assess_gene_with_threshold(standardised_gene_result, gene, rank_stats)
                if self.score_order != "ascending"
                else self._assess_gene_with_threshold_ascending_order(
                    standardised_gene_result, gene, rank_stats
                )
            )

    def assess_gene_prioritisation(self, rank_stats: RankStats, rank_records: defaultdict) -> None:
        """Assess gene prioritisation."""
        for gene in self.proband_causative_genes:
            rank_stats.total += 1
            gene_match = GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol)
            for standardised_gene_result in self.standardised_gene_results:
                if (
                    gene.gene_identifier == standardised_gene_result.gene_identifier
                    or gene.gene_symbol == standardised_gene_result.gene_symbol
                ):
                    gene_match = self._record_matched_gene(
                        gene, rank_stats, standardised_gene_result
                    )
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                GenePrioritisationResult(self.phenopacket_path, gene.gene_symbol)
                if gene_match is None
                else gene_match,
                rank_records,
            ).record_rank()


class AssessVariantPrioritisation:
    """Assess variant prioritisation."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_variant_results: [RankedPhEvalVariantResult],
        threshold: float,
        score_order: str,
        proband_causative_variants: [GenomicVariant],
    ):
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_variant_results = standardised_variant_results
        self.threshold = threshold
        self.score_order = score_order
        self.proband_causative_variants = proband_causative_variants

    def _record_variant_prioritisation_match(
        self,
        result_entry: RankedPhEvalVariantResult,
        rank_stats: RankStats,
    ) -> VariantPrioritisationResult:
        """Record the variant prioritisation rank if found within results."""
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return VariantPrioritisationResult(
            self.phenopacket_path,
            GenomicVariant(
                chrom=result_entry.chromosome,
                pos=result_entry.start,
                ref=result_entry.ref,
                alt=result_entry.alt,
            ),
            rank,
        )

    def _assess_variant_with_threshold_ascending_order(
        self, result_entry: RankedPhEvalVariantResult, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """Record the variant prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _assess_variant_with_threshold(
        self, result_entry: pd.Series, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """Record the variant prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _record_matched_variant(
        self, rank_stats: RankStats, standardised_variant_result: pd.Series
    ) -> VariantPrioritisationResult:
        """Return the variant rank result - dealing with the specification of a threshold."""
        if float(self.threshold) == 0.0:
            return self._record_variant_prioritisation_match(
                standardised_variant_result, rank_stats
            )
        else:
            return (
                self._assess_variant_with_threshold(standardised_variant_result, rank_stats)
                if self.score_order != "ascending"
                else self._assess_variant_with_threshold_ascending_order(
                    standardised_variant_result, rank_stats
                )
            )

    def assess_variant_prioritisation(
        self, rank_stats: RankStats, rank_records: defaultdict
    ) -> None:
        """Assess variant prioritisation."""
        for variant in self.proband_causative_variants:
            rank_stats.total += 1
            variant_match = VariantPrioritisationResult(self.phenopacket_path, variant)
            for result in self.standardised_variant_results:
                result_variant = GenomicVariant(
                    chrom=result.chromosome,
                    pos=result.start,
                    ref=result.ref,
                    alt=result.alt,
                )
                if variant == result_variant:
                    variant_match = self._record_matched_variant(rank_stats, result)
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                VariantPrioritisationResult(self.phenopacket_path, variant)
                if variant_match is None
                else variant_match,
                rank_records,
            ).record_rank()


class AssessDiseasePrioritisation:
    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_disease_results: [RankedPhEvalDiseaseResult],
        threshold: float,
        score_order: str,
        proband_diseases: [ProbandDisease],
    ):
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_disease_results = standardised_disease_results
        self.threshold = threshold
        self.score_order = score_order
        self.proband_diseases = proband_diseases

    def _record_disease_prioritisation_match(
        self,
        disease: ProbandDisease,
        result_entry: RankedPhEvalDiseaseResult,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """Record the disease prioritisation rank if found within results."""
        rank = result_entry.rank
        rank_stats.add_rank(rank)
        return DiseasePrioritisationResult(self.phenopacket_path, disease, rank)

    def _assess_disease_with_threshold_ascending_order(
        self,
        result_entry: RankedPhEvalDiseaseResult,
        disease: ProbandDisease,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """Record the disease prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.score):
            return self._record_disease_prioritisation_match(disease, result_entry, rank_stats)

    def _assess_disease_with_threshold(
        self,
        result_entry: RankedPhEvalDiseaseResult,
        disease: ProbandDisease,
        rank_stats: RankStats,
    ) -> DiseasePrioritisationResult:
        """Record the disease prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.score):
            return self._record_disease_prioritisation_match(disease, result_entry, rank_stats)

    def _record_matched_disease(
        self,
        disease: ProbandDisease,
        rank_stats: RankStats,
        standardised_disease_result: RankedPhEvalDiseaseResult,
    ) -> DiseasePrioritisationResult:
        """Return the gene rank result - dealing with the specification of a threshold."""
        if float(self.threshold) == 0.0:
            return self._record_disease_prioritisation_match(
                disease, standardised_disease_result, rank_stats
            )
        else:
            return (
                self._assess_disease_with_threshold(
                    standardised_disease_result, disease, rank_stats
                )
                if self.score_order != "ascending"
                else self._assess_disease_with_threshold_ascending_order(
                    standardised_disease_result, disease, rank_stats
                )
            )

    def assess_disease_prioritisation(
        self, rank_stats: RankStats, rank_records: defaultdict
    ) -> None:
        """Assess disease prioritisation."""
        for disease in self.proband_diseases:
            rank_stats.total += 1
            disease_match = DiseasePrioritisationResult(self.phenopacket_path, disease)
            for standardised_disease_result in self.standardised_disease_results:
                if (
                    disease.disease_identifier == standardised_disease_result.disease_identifier
                    or disease.disease_name == standardised_disease_result.disease_name
                ):
                    disease_match = self._record_matched_disease(
                        disease, rank_stats, standardised_disease_result
                    )
                    break
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                DiseasePrioritisationResult(self.phenopacket_path, disease)
                if disease_match is None
                else disease_match,
                rank_records,
            ).record_rank()


def _obtain_causative_genes(phenopacket_path: Path) -> [ProbandCausativeGene]:
    """Obtain causative genes from a phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_genes()


def _obtain_causative_variants(phenopacket_path: Path) -> [GenomicVariant]:
    """Obtain causative variants from a phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_variants()


def _obtain_causative_diseases(phenopacket_path: Path) -> [ProbandDisease]:
    """Obtain known diseases from a phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnoses()


def _assess_phenopacket_gene_prioritisation(
    standardised_gene_result: Path,
    score_order: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    gene_rank_stats: RankStats,
    gene_rank_comparison: defaultdict,
) -> None:
    """Assess gene prioritisation for a phenopacket."""
    phenopacket_path = obtain_closest_file_name(
        standardised_gene_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    pheval_gene_result = read_standardised_result(standardised_gene_result)
    proband_causative_genes = _obtain_causative_genes(phenopacket_path)
    AssessGenePrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_gene_results/"),
        parse_pheval_result(RankedPhEvalGeneResult, pheval_gene_result),
        threshold,
        score_order,
        proband_causative_genes,
    ).assess_gene_prioritisation(gene_rank_stats, gene_rank_comparison)


def _assess_phenopacket_variant_prioritisation(
    standardised_variant_result: Path,
    score_order: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    variant_rank_stats: RankStats,
    variant_rank_comparison: defaultdict,
) -> None:
    """Assess variant prioritisation for a phenopacket"""
    phenopacket_path = obtain_closest_file_name(
        standardised_variant_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    proband_causative_variants = _obtain_causative_variants(phenopacket_path)
    pheval_variant_result = read_standardised_result(standardised_variant_result)
    AssessVariantPrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_variant_results/"),
        parse_pheval_result(RankedPhEvalVariantResult, pheval_variant_result),
        threshold,
        score_order,
        proband_causative_variants,
    ).assess_variant_prioritisation(variant_rank_stats, variant_rank_comparison)


def _assess_phenopacket_disease_prioritisation(
    standardised_disease_result: Path,
    score_order: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    disease_rank_stats: RankStats,
    disease_rank_comparison: defaultdict,
) -> None:
    """Assess gene prioritisation for a phenopacket."""
    phenopacket_path = obtain_closest_file_name(
        standardised_disease_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    pheval_disease_result = read_standardised_result(standardised_disease_result)
    proband_diseases = _obtain_causative_diseases(phenopacket_path)
    AssessDiseasePrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_disease_results/"),
        parse_pheval_result(RankedPhEvalDiseaseResult, pheval_disease_result),
        threshold,
        score_order,
        proband_diseases,
    ).assess_disease_prioritisation(disease_rank_stats, disease_rank_comparison)


def _assess_prioritisation_for_results_directory(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    gene_rank_comparison: defaultdict,
    variant_rank_comparison: defaultdict,
    disease_rank_comparison: defaultdict,
    gene_stats_writer: RankStatsWriter,
    variants_stats_writer: RankStatsWriter,
    disease_stats_writer: RankStatsWriter,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
) -> TrackRunPrioritisation:
    """Assess prioritisation for a single results directory."""
    gene_rank_stats, variant_rank_stats, disease_rank_stats = RankStats(), RankStats(), RankStats()
    if gene_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_gene_results/"), ".tsv"
        ):
            _assess_phenopacket_gene_prioritisation(
                standardised_result,
                score_order,
                results_directory_and_input,
                threshold,
                gene_rank_stats,
                gene_rank_comparison,
            )
        gene_stats_writer.write_row(results_directory_and_input.results_dir, gene_rank_stats)
    if variant_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_variant_results/"),
            ".tsv",
        ):
            _assess_phenopacket_variant_prioritisation(
                standardised_result,
                score_order,
                results_directory_and_input,
                threshold,
                variant_rank_stats,
                variant_rank_comparison,
            )
        variants_stats_writer.write_row(results_directory_and_input.results_dir, variant_rank_stats)
    if disease_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_disease_results/"),
            ".tsv",
        ):
            _assess_phenopacket_disease_prioritisation(
                standardised_result,
                score_order,
                results_directory_and_input,
                threshold,
                disease_rank_stats,
                disease_rank_comparison,
            )
        disease_stats_writer.write_row(results_directory_and_input.results_dir, disease_rank_stats)
    return TrackRunPrioritisation(
        gene_prioritisation=AnalysisResults(
            results_dir=results_directory_and_input.results_dir,
            ranks=gene_rank_comparison,
            rank_stats=gene_rank_stats,
        ),
        variant_prioritisation=AnalysisResults(
            results_dir=results_directory_and_input.results_dir,
            ranks=variant_rank_comparison,
            rank_stats=variant_rank_stats,
        ),
        disease_prioritisation=AnalysisResults(
            results_dir=results_directory_and_input.results_dir,
            ranks=disease_rank_comparison,
            rank_stats=disease_rank_stats,
        ),
    )


def benchmark_directory(
    results_dir_and_input: TrackInputOutputDirectories,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
) -> None:
    """Benchmark prioritisation performance for a single directory."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    disease_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-disease_summary.tsv")) if disease_analysis else None
    )
    gene_rank_comparison, variant_rank_comparison, disease_rank_comparison = (
        defaultdict(dict),
        defaultdict(dict),
        defaultdict(dict),
    )
    prioritisation_data = _assess_prioritisation_for_results_directory(
        results_dir_and_input,
        score_order,
        threshold,
        gene_rank_comparison,
        variant_rank_comparison,
        disease_rank_comparison,
        gene_stats_writer,
        variants_stats_writer,
        disease_stats_writer,
        gene_analysis,
        variant_analysis,
        disease_analysis,
    )
    generate_benchmark_gene_output(prioritisation_data, plot_type) if gene_analysis else None
    generate_benchmark_variant_output(prioritisation_data, plot_type) if variant_analysis else None
    generate_benchmark_disease_output(prioritisation_data, plot_type) if disease_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None
    disease_stats_writer.close() if disease_analysis else None


def benchmark_runs(
    results_directories: [TrackInputOutputDirectories],
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
) -> None:
    """Benchmark several result directories."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    disease_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-disease_summary.tsv")) if disease_analysis else None
    )
    prioritisation_stats_for_runs = []
    for results_dir_and_input in results_directories:
        gene_rank_comparison, variant_rank_comparison, disease_rank_comparison = (
            defaultdict(dict),
            defaultdict(dict),
            defaultdict(dict),
        )
        prioritisation_stats = _assess_prioritisation_for_results_directory(
            results_dir_and_input,
            score_order,
            threshold,
            gene_rank_comparison,
            variant_rank_comparison,
            disease_rank_comparison,
            gene_stats_writer,
            variants_stats_writer,
            disease_stats_writer,
            gene_analysis,
            variant_analysis,
            disease_analysis,
        )
        prioritisation_stats_for_runs.append(prioritisation_stats)
    generate_benchmark_comparison_gene_output(
        prioritisation_stats_for_runs, plot_type
    ) if gene_analysis else None
    generate_benchmark_comparison_variant_output(
        prioritisation_stats_for_runs, plot_type
    ) if variant_analysis else None
    generate_benchmark_comparison_disease_output(
        prioritisation_stats_for_runs, plot_type
    ) if disease_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None
    disease_stats_writer.close() if disease_analysis else None


@click.command()
@click.option(
    "--directory",
    "-d",
    required=True,
    metavar="PATH",
    help="General results directory to be benchmarked, assumes contains subdirectories of pheval_gene_results/"
    "pheval_variant_results and the tool specific results directory. ",
    type=Path,
)
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    metavar="PATH",
    help="Full path to directory containing input phenopackets.",
    type=Path,
)
@click.option(
    "--output-prefix",
    "-o",
    metavar="<str>",
    required=True,
    help=" Output file prefix. ",
)
@click.option(
    "--score-order",
    "-so",
    required=True,
    help="Ordering of results for ranking.",
    type=click.Choice(["ascending", "descending"]),
    default="descending",
    show_default=True,
)
@click.option(
    "--threshold",
    "-t",
    metavar="<float>",
    default=float(0.0),
    required=False,
    help="Score threshold.",
    type=float,
)
@click.option(
    "--gene-analysis/--no-gene-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for variant prioritisation",
)
@click.option(
    "--disease-analysis/--no-disease-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for disease prioritisation",
)
@click.option(
    "--plot-type",
    "-y",
    default="bar_stacked",
    show_default=True,
    type=click.Choice(["bar_stacked", "bar_cumulative", "bar_non_cumulative"]),
    help="Bar chart type to output.",
)
def benchmark(
    directory: Path,
    phenopacket_dir: Path,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
):
    """Benchmark the gene/variant prioritisation performance for a single run."""
    if not gene_analysis and not variant_analysis and not disease_analysis:
        raise InputError("Need to specify gene analysis and/or variant and/or disease analysis.")
    benchmark_directory(
        TrackInputOutputDirectories(results_dir=directory, phenopacket_dir=phenopacket_dir),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
        disease_analysis,
        plot_type,
    )


@click.command()
@click.option(
    "--run-data",
    "-r",
    required=True,
    metavar="PATH",
    help="Path to .txt file containing testdata directory and corresponding results directory separated by tab."
    "Each run contained to a new line with the input testdata listed first and on the same line separated by a tab"
    "the results directory.",
    type=Path,
)
@click.option(
    "--output-prefix",
    "-o",
    metavar="<str>",
    required=True,
    help=" Output file prefix. ",
)
@click.option(
    "--score-order",
    "-so",
    required=True,
    help="Ordering of results for ranking.",
    type=click.Choice(["ascending", "descending"]),
    default="descending",
    show_default=True,
)
@click.option(
    "--threshold",
    "-t",
    metavar="<float>",
    default=float(0.0),
    required=False,
    help="Score threshold.",
    type=float,
)
@click.option(
    "--gene-analysis/--no-gene-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for variant prioritisation",
)
@click.option(
    "--disease-analysis/--no-disease-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for disease prioritisation",
)
@click.option(
    "--plot-type",
    "-y",
    default="bar_stacked",
    show_default=True,
    type=click.Choice(["bar_stacked", "bar_cumulative", "bar_non_cumulative"]),
    help="Bar chart type to output.",
)
def benchmark_comparison(
    run_data: Path,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
):
    """Benchmark the gene/variant prioritisation performance for two runs."""
    if not gene_analysis and not variant_analysis and not disease_analysis:
        raise InputError("Need to specify gene analysis and/or variant and/or disease analysis.")
    benchmark_runs(
        _parse_run_data_text_file(run_data),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
        disease_analysis,
        plot_type,
    )
