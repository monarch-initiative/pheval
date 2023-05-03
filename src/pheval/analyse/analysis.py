# #!/usr/bin/python
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import click
import pandas as pd

from pheval.analyse.generate_plots import (
    TrackGenePrioritisation,
    TrackPrioritisation,
    TrackVariantPrioritisation,
)
from pheval.analyse.generate_summary_outputs import (
    RankStatsWriter,
    generate_benchmark_comparison_gene_output,
    generate_benchmark_comparison_variant_output,
    generate_benchmark_gene_output,
    generate_benchmark_variant_output,
)
from pheval.analyse.rank_stats import RankStats
from pheval.post_processing.post_processing import (
    PhEvalGeneResult,
    PhEvalVariantResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)
from pheval.prepare.custom_exceptions import InputError
from pheval.utils.file_utils import all_files, files_with_suffix, obtain_closest_file_name
from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    PhenopacketUtil,
    ProbandCausativeGene,
    phenopacket_reader,
)


def _read_standardised_result(standardised_result_path: Path) -> dict:
    """Read the standardised result output and return a dictionary."""
    return pd.read_csv(standardised_result_path, delimiter="\t")


def parse_pheval_gene_result(pheval_gene_result: pd.DataFrame) -> [RankedPhEvalGeneResult]:
    """Parse PhEval gene result into RankedPhEvalGeneResult dataclass."""
    ranked_gene_results = []
    for _index, result in pheval_gene_result.iterrows():
        ranked_gene_results.append(
            RankedPhEvalGeneResult(
                pheval_gene_result=PhEvalGeneResult(
                    gene_symbol=result["gene_symbol"],
                    gene_identifier=result["gene_identifier"],
                    score=result["score"],
                ),
                rank=result["rank"],
            )
        )
    return ranked_gene_results


def parse_pheval_variant_result(pheval_variant_result: pd.DataFrame) -> [RankedPhEvalVariantResult]:
    """Parse PhEval variant result into RankedPhEvalVariantResult dataclass."""
    ranked_variant_results = []
    for _index, result in pheval_variant_result.iterrows():
        ranked_variant_results.append(
            RankedPhEvalVariantResult(
                pheval_variant_result=PhEvalVariantResult(
                    chromosome=result["chromosome"],
                    start=result["start"],
                    end=result["end"],
                    ref=result["ref"],
                    alt=result["alt"],
                    score=result["score"],
                ),
                rank=result["rank"],
            )
        )
    return ranked_variant_results


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

    def record_rank(self) -> None:
        """Records the rank for different runs."""
        self.run_comparison[self.index][
            "Phenopacket"
        ] = self.prioritisation_result.phenopacket_path.name
        self._record_gene_rank() if type(
            self.prioritisation_result
        ) is GenePrioritisationResult else self._record_variant_rank()
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
        if float(self.threshold) > float(result_entry.pheval_gene_result.score):
            return self._record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def _assess_gene_with_threshold(
        self,
        result_entry: RankedPhEvalGeneResult,
        gene: ProbandCausativeGene,
        rank_stats: RankStats,
    ) -> GenePrioritisationResult:
        """Record the gene prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.pheval_gene_result.score):
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
                    gene.gene_identifier
                    == standardised_gene_result.pheval_gene_result.gene_identifier
                    or gene.gene_symbol
                    == standardised_gene_result.pheval_gene_result.gene_identifier
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
                chrom=result_entry.pheval_variant_result.chromosome,
                pos=result_entry.pheval_variant_result.start,
                ref=result_entry.pheval_variant_result.ref,
                alt=result_entry.pheval_variant_result.alt,
            ),
            rank,
        )

    def _assess_variant_with_threshold_ascending_order(
        self, result_entry: RankedPhEvalVariantResult, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """Record the variant prioritisation rank if it meets the ascending order threshold."""
        if float(self.threshold) > float(result_entry.pheval_variant_result.score):
            return self._record_variant_prioritisation_match(result_entry, rank_stats)

    def _assess_variant_with_threshold(
        self, result_entry: pd.Series, rank_stats: RankStats
    ) -> VariantPrioritisationResult:
        """Record the variant prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry.pheval_variant_result.score):
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
                    chrom=result.pheval_variant_result.chromosome,
                    pos=result.pheval_variant_result.start,
                    ref=result.pheval_variant_result.ref,
                    alt=result.pheval_variant_result.alt,
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
    pheval_gene_result = _read_standardised_result(standardised_gene_result)
    proband_causative_genes = _obtain_causative_genes(phenopacket_path)
    AssessGenePrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_gene_results/"),
        parse_pheval_gene_result(pheval_gene_result),
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
    pheval_variant_result = _read_standardised_result(standardised_variant_result)
    AssessVariantPrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_variant_results/"),
        parse_pheval_variant_result(pheval_variant_result),
        threshold,
        score_order,
        proband_causative_variants,
    ).assess_variant_prioritisation(variant_rank_stats, variant_rank_comparison)


def _assess_prioritisation_for_results_directory(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    gene_rank_comparison: defaultdict,
    variant_rank_comparison: defaultdict,
    gene_stats_writer: RankStatsWriter,
    variants_stats_writer: RankStatsWriter,
    gene_analysis: bool,
    variant_analysis: bool,
) -> TrackPrioritisation:
    """Assess prioritisation for a single results directory."""
    gene_rank_stats, variant_rank_stats = RankStats(), RankStats()
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
    gene_stats_writer.write_row(
        results_directory_and_input.results_dir, gene_rank_stats
    ) if gene_analysis else None
    variants_stats_writer.write_row(
        results_directory_and_input.results_dir, variant_rank_stats
    ) if variant_analysis else None
    return TrackPrioritisation(
        gene_prioritisation=TrackGenePrioritisation(
            results_dir=results_directory_and_input.results_dir,
            ranks=gene_rank_comparison,
            rank_stats=gene_rank_stats,
        ),
        variant_prioritisation=TrackVariantPrioritisation(
            results_dir=results_directory_and_input.results_dir,
            ranks=variant_rank_comparison,
            rank_stats=variant_rank_stats,
        ),
    )


def benchmark_directory(
    results_dir_and_input: TrackInputOutputDirectories,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    plot_type: str,
) -> None:
    """Benchmark prioritisation performance for a single directory."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    gene_rank_comparison, variant_rank_comparison = defaultdict(dict), defaultdict(dict)
    prioritisation_data = _assess_prioritisation_for_results_directory(
        results_dir_and_input,
        score_order,
        threshold,
        gene_rank_comparison,
        variant_rank_comparison,
        gene_stats_writer,
        variants_stats_writer,
        gene_analysis,
        variant_analysis,
    )
    generate_benchmark_gene_output(prioritisation_data, plot_type) if gene_analysis else None
    generate_benchmark_variant_output(prioritisation_data, plot_type) if variant_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None


def benchmark_runs(
    results_directories: [TrackInputOutputDirectories],
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    plot_type: str,
) -> None:
    """Benchmark several result directories."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    prioritisation_stats_for_runs = []
    for results_dir_and_input in results_directories:
        gene_rank_comparison, variant_rank_comparison = defaultdict(dict), defaultdict(dict)
        prioritisation_stats = _assess_prioritisation_for_results_directory(
            results_dir_and_input,
            score_order,
            threshold,
            gene_rank_comparison,
            variant_rank_comparison,
            gene_stats_writer,
            variants_stats_writer,
            gene_analysis,
            variant_analysis,
        )
        prioritisation_stats_for_runs.append(prioritisation_stats)
    generate_benchmark_comparison_gene_output(
        prioritisation_stats_for_runs, plot_type
    ) if gene_analysis else None
    generate_benchmark_comparison_variant_output(
        prioritisation_stats_for_runs, plot_type
    ) if variant_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None


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
    "--plot-type",
    "-p",
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
    plot_type: str,
):
    """Benchmark the gene/variant prioritisation performance for a single run."""
    if not gene_analysis and not variant_analysis:
        raise InputError("Need to specify gene analysis and/or variant analysis.")
    benchmark_directory(
        TrackInputOutputDirectories(results_dir=directory, phenopacket_dir=phenopacket_dir),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
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
    "--plot-type",
    "-p",
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
    plot_type: str,
):
    """Benchmark the gene/variant prioritisation performance for two runs."""
    if not gene_analysis and not variant_analysis:
        raise InputError("Need to specify gene analysis and/or variant analysis.")
    benchmark_runs(
        _parse_run_data_text_file(run_data),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
        plot_type,
    )
