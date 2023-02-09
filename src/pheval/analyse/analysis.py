# #!/usr/bin/python
import csv
import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean

import click
import pandas as pd

from pheval.utils.file_utils import all_files, files_with_suffix, obtain_closest_file_name
from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    PhenopacketUtil,
    ProbandCausativeGene,
    phenopacket_reader,
)


def read_standardised_result(standardised_result_path: Path) -> dict:
    """Read the standardised result output and return a dictionary."""
    return pd.read_csv(standardised_result_path, delimiter="\t")


@dataclass
class GenePrioritisationResultData:
    """Store rank data for causative genes."""

    phenopacket: Path
    gene: str
    rank: int = 0


@dataclass
class VariantPrioritisationResultData:
    """Store rank data for causative variants."""

    phenopacket: Path
    variant: GenomicVariant
    rank: int = 0


@dataclass
class PrioritisationRankRecorder:
    """Compare the ranks of different runs."""

    index: int
    directory: Path
    prioritisation_run_comparison: VariantPrioritisationResultData or GenePrioritisationResultData
    run_comparison: defaultdict

    def record_rank(self) -> None:
        """Records the rank for different runs."""
        self.run_comparison[self.index][
            "Phenopacket"
        ] = self.prioritisation_run_comparison.phenopacket.name
        if type(self.prioritisation_run_comparison) is GenePrioritisationResultData:
            self.run_comparison[self.index]["Gene"] = self.prioritisation_run_comparison.gene
        if type(self.prioritisation_run_comparison) is VariantPrioritisationResultData:
            variant_info = self.prioritisation_run_comparison.variant
            self.run_comparison[self.index]["Variant"] = "_".join(
                [variant_info.chrom, str(variant_info.pos), variant_info.ref, variant_info.alt]
            )
        self.run_comparison[self.index][self.directory] = self.prioritisation_run_comparison.rank


class RankComparisonGenerator:
    """Write the run comparison of rank assignment for prioritisation."""

    def __init__(self, run_comparison: defaultdict):
        self.run_comparison = run_comparison

    def generate_dataframe(self) -> pd.DataFrame:
        """Generate pandas dataframe."""
        return pd.DataFrame.from_dict(self.run_comparison, orient="index")

    def calculate_rank_difference(self) -> pd.DataFrame:
        """Calculate the rank decrease for runs - taking the first directory as a baseline."""
        comparison_df = self.generate_dataframe()
        comparison_df["rank_decrease"] = comparison_df.iloc[:, 3] - comparison_df.iloc[:, 2]
        return comparison_df

    def generate_gene_output(self, prefix: str) -> None:
        """Generate the output for gene prioritisation ranks."""
        self.generate_dataframe().to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")

    def generate_variant_output(self, prefix: str) -> None:
        """Generate the output for variant prioritisation ranks."""
        self.generate_dataframe().to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")

    def generate_gene_comparison_output(self, prefix: str) -> None:
        """Generate the output for gene prioritisation rank comparison."""
        self.calculate_rank_difference().to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")

    def generate_variant_comparison_output(self, prefix: str) -> None:
        """Generate the output for variant prioritisation rank comparison."""
        self.calculate_rank_difference().to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")


@dataclass
class RankStats:
    """Class for keeping track of the rank stats."""

    top: int = 0
    top3: int = 0
    top5: int = 0
    found: int = 0
    total: int = 0
    reciprocal_ranks: list = field(default_factory=list)

    def add_rank(self, rank: int) -> None:
        """Add rank for phenopacket."""
        self.reciprocal_ranks.append(1 / rank)
        self.found += 1
        if rank == 1:
            self.top += 1
        if rank != "" and rank <= 3:
            self.top3 += 1
        if rank != "" and rank <= 5:
            self.top5 += 1

    def percentage_rank(self, value: int) -> float:
        """Return a percentage rank."""
        return 100 * value / self.found

    def percentage_top(self) -> float:
        """Return percentage of top matches."""
        return self.percentage_rank(self.top)

    def percentage_top3(self) -> float:
        """Return percentage of matches in the top3."""
        return self.percentage_rank(self.top3)

    def percentage_top5(self) -> float:
        """Return percentage of matches in the top5."""
        return self.percentage_rank(self.top5)

    def percentage_found(self) -> float:
        """Return percentage of matches found."""
        return 100 * self.found / self.total

    def mean_reciprocal_rank(self) -> float:
        """Return the mean reciprocal rank."""
        return mean(self.reciprocal_ranks)


class RankStatsWriter:
    """Write the rank stats for each run."""

    def __init__(self, file: Path):
        self.file = open(file, "w")
        self.writer = csv.writer(self.file, delimiter="\t")
        self.writer.writerow(
            [
                "results_directory_path",
                "top",
                "top3",
                "top5",
                "found",
                "total",
                "mean_reciprocal_rank",
                "percentage_top",
                "percentage_top3",
                "percentage_top5",
                "percentage_found",
            ]
        )

    def write_row(self, directory: Path, rank_stats: RankStats) -> None:
        """Write summary rank stats row for run."""
        try:
            self.writer.writerow(
                [
                    directory,
                    rank_stats.top,
                    rank_stats.top3,
                    rank_stats.top5,
                    rank_stats.found,
                    rank_stats.total,
                    rank_stats.mean_reciprocal_rank(),
                    rank_stats.percentage_top(),
                    rank_stats.percentage_top3(),
                    rank_stats.percentage_top5(),
                    rank_stats.percentage_found(),
                ]
            )
        except IOError:
            print("Error writing ", self.file)

    def close(self):
        """Close file."""
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


@dataclass
class TrackInputOutputDirectories:
    """Track the input testdata for a corresponding standardised output directory"""

    phenopacket_dir: Path
    results_dir: Path


class AssessGenePrioritisation:
    """Assess gene prioritisation."""

    def __init__(
        self,
        phenopacket_path: Path,
        results_dir: Path,
        standardised_gene_results: [dict],
        threshold: float,
        ranking_method: str,
        proband_causative_genes: [ProbandCausativeGene],
    ):
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_gene_results = standardised_gene_results
        self.threshold = threshold
        self.ranking_method = ranking_method
        self.proband_causative_genes = proband_causative_genes

    def record_gene_prioritisation_match(
        self, gene: ProbandCausativeGene, result_entry: dict, rank_stats: RankStats
    ) -> GenePrioritisationResultData:
        """Record the gene prioritisation rank if found within results."""
        rank = result_entry["rank"]
        rank_stats.add_rank(rank)
        gene_match = GenePrioritisationResultData(self.phenopacket_path, gene.gene_symbol, rank)
        return gene_match

    def assess_gene_with_pvalue_threshold(
        self, result_entry: dict, gene: ProbandCausativeGene, rank_stats: RankStats
    ) -> GenePrioritisationResultData:
        """Record the gene prioritisation rank if it meets the pvalue threshold."""
        if float(self.threshold) > float(result_entry["score"]):
            return self.record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def assess_gene_with_threshold(
        self, result_entry: dict, gene: ProbandCausativeGene, rank_stats: RankStats
    ) -> GenePrioritisationResultData:
        """Record the gene prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry["score"]):
            return self.record_gene_prioritisation_match(gene, result_entry, rank_stats)

    def record_matched_gene(self, gene, rank_stats, standardised_gene_result):
        """Return the gene rank result - dealing with the specification of a threshold."""
        if float(self.threshold) == 0.0:
            return self.record_gene_prioritisation_match(gene, standardised_gene_result, rank_stats)
        else:
            return (
                self.assess_gene_with_threshold(standardised_gene_result, gene, rank_stats)
                if self.ranking_method != "pValue"
                else self.assess_gene_with_pvalue_threshold(
                    standardised_gene_result, gene, rank_stats
                )
            )

    def assess_gene_prioritisation(self, rank_stats: RankStats, rank_records: defaultdict):
        """Assess gene prioritisation."""
        for gene in self.proband_causative_genes:
            rank_stats.total += 1
            gene_match = GenePrioritisationResultData(self.phenopacket_path, gene.gene_symbol)
            for _index, standardised_gene_result in self.standardised_gene_results.iterrows():
                if (
                    gene.gene_identifier == standardised_gene_result["gene_identifier"]
                    or gene.gene_symbol == standardised_gene_result["gene_symbol"]
                ):
                    gene_match = self.record_matched_gene(
                        gene, rank_stats, standardised_gene_result
                    )
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                GenePrioritisationResultData(self.phenopacket_path, gene.gene_symbol)
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
        standardised_variant_results: [dict],
        threshold: float,
        ranking_method: str,
        proband_causative_variants: [GenomicVariant],
    ):
        self.phenopacket_path = phenopacket_path
        self.results_dir = results_dir
        self.standardised_variant_results = standardised_variant_results
        self.threshold = threshold
        self.ranking_method = ranking_method
        self.proband_causative_variants = proband_causative_variants

    def record_variant_prioritisation_match(
        self,
        result_entry: pd.Series,
        rank_stats: RankStats,
    ) -> VariantPrioritisationResultData:
        """Record the variant prioritisation rank if found within results."""
        rank = result_entry["rank"]
        rank_stats.add_rank(rank)
        variant_match = VariantPrioritisationResultData(
            self.phenopacket_path,
            GenomicVariant(
                chrom=result_entry["chromosome"],
                pos=result_entry["start"],
                ref=result_entry["ref"],
                alt=result_entry["alt"],
            ),
            rank,
        )
        return variant_match

    def assess_variant_with_pvalue_threshold(
        self, result_entry: pd.Series, rank_stats: RankStats
    ) -> VariantPrioritisationResultData:
        """Record the variant prioritisation rank if it meets the pvalue threshold."""
        if float(self.threshold) > float(result_entry["score"]):
            return self.record_variant_prioritisation_match(result_entry, rank_stats)

    def assess_variant_with_threshold(
        self, result_entry: pd.Series, rank_stats: RankStats
    ) -> VariantPrioritisationResultData:
        """Record the variant prioritisation rank if it meets the score threshold."""
        if float(self.threshold) < float(result_entry["score"]):
            return self.record_variant_prioritisation_match(result_entry, rank_stats)

    def record_matched_variant(self, rank_stats, result) -> VariantPrioritisationResultData:
        """Return the variant rank result - dealing with the specification of a threshold."""
        if float(self.threshold) == 0.0:
            return self.record_variant_prioritisation_match(result, rank_stats)
        else:
            return (
                self.assess_variant_with_threshold(result, rank_stats)
                if self.ranking_method != "pValue"
                else self.assess_variant_with_pvalue_threshold(result, rank_stats)
            )

    def assess_variant_prioritisation(self, rank_stats: RankStats, rank_records: defaultdict):
        """Assess variant prioritisation."""
        for variant in self.proband_causative_variants:
            rank_stats.total += 1
            variant_match = VariantPrioritisationResultData(self.phenopacket_path, variant)
            for _index, result in self.standardised_variant_results.iterrows():
                result_variant = GenomicVariant(
                    chrom=result["chromosome"], pos=result["start"], ref=result["ref"], alt=result["alt"]
                )
                if variant == result_variant:
                    variant_match = self.record_matched_variant(rank_stats, result)
            PrioritisationRankRecorder(
                rank_stats.total,
                self.results_dir,
                VariantPrioritisationResultData(self.phenopacket_path, variant)
                if variant_match is None
                else variant_match,
                rank_records,
            ).record_rank()


def obtain_causative_genes(phenopacket_path):
    """Obtain causative genes from a phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_genes()


def obtain_causative_variants(phenopacket_path):
    """Obtain causative variants from a phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_variants()


def assess_phenopacket_gene_prioritisation(
    standardised_gene_result: Path,
    ranking_method: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    gene_rank_stats: RankStats,
    gene_rank_comparison: defaultdict,
):
    """Assess gene prioritisation for a phenopacket."""
    phenopacket_path = obtain_closest_file_name(
        standardised_gene_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    proband_causative_genes = obtain_causative_genes(phenopacket_path)
    AssessGenePrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_gene_results/"),
        read_standardised_result(standardised_gene_result),
        threshold,
        ranking_method,
        proband_causative_genes,
    ).assess_gene_prioritisation(gene_rank_stats, gene_rank_comparison)


def assess_phenopacket_variant_prioritisation(
    standardised_variant_result: Path,
    ranking_method: str,
    results_dir_and_input: TrackInputOutputDirectories,
    threshold: float,
    variant_rank_stats: RankStats,
    variant_rank_comparison: defaultdict,
):
    """Assess variant prioritisation for a phenopacket"""
    phenopacket_path = obtain_closest_file_name(
        standardised_variant_result, all_files(results_dir_and_input.phenopacket_dir)
    )
    proband_causative_variants = obtain_causative_variants(phenopacket_path)
    AssessVariantPrioritisation(
        phenopacket_path,
        results_dir_and_input.results_dir.joinpath("pheval_variant_results/"),
        read_standardised_result(standardised_variant_result),
        threshold,
        ranking_method,
        proband_causative_variants,
    ).assess_variant_prioritisation(variant_rank_stats, variant_rank_comparison)


def assess_prioritisation_for_results_directory(
    results_directory_and_input: TrackInputOutputDirectories,
    ranking_method: str,
    threshold: float,
    gene_rank_comparison: defaultdict,
    variant_rank_comparison: defaultdict,
    gene_stats_writer: RankStatsWriter,
    variants_stats_writer: RankStatsWriter,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Assess prioritisation for a single results directory."""
    gene_rank_stats, variant_rank_stats = RankStats(), RankStats()
    if gene_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_gene_results/"), ".tsv"
        ):
            assess_phenopacket_gene_prioritisation(
                standardised_result,
                ranking_method,
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
            assess_phenopacket_variant_prioritisation(
                standardised_result,
                ranking_method,
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
    return gene_rank_comparison, variant_rank_comparison


def benchmark_directory(
    results_dir_and_input: TrackInputOutputDirectories,
    ranking_method: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Benchmark prioritisation performance for a single directory."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    gene_rank_comparison, variant_rank_comparison = defaultdict(dict), defaultdict(dict)
    assess_prioritisation_for_results_directory(
        results_dir_and_input,
        ranking_method,
        threshold,
        gene_rank_comparison,
        variant_rank_comparison,
        gene_stats_writer,
        variants_stats_writer,
        gene_analysis,
        variant_analysis,
    )
    RankComparisonGenerator(gene_rank_comparison).generate_gene_output(
        f"{results_dir_and_input.results_dir.name}"
    ) if gene_analysis else None
    RankComparisonGenerator(variant_rank_comparison).generate_variant_output(
        f"{results_dir_and_input.results_dir.name}"
    ) if variant_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None


def merge_results(result1, result2):
    """Merge two nested dictionaries containing results on commonalities."""
    for key, val in result1.items():
        if type(val) == dict:
            if key in result2 and type(result2[key] == dict):
                merge_results(result1[key], result2[key])
        else:
            if key in result2:
                result1[key] = result2[key]

    for key, val in result2.items():
        if key not in result1:
            result1[key] = val
    return result1


@dataclass
class TrackGeneComparisons:
    """Track the gene ranks for each result in a result directory."""

    directory: Path
    gene_results: dict


@dataclass
class TrackVariantComparisons:
    """Track the variant ranks for each result in a result directory."""

    directory: Path
    variant_results: dict


def generate_gene_rank_comparisons(comparison_ranks: [tuple]) -> None:
    """Generate the gene rank comparison of two result directories."""
    for pair in comparison_ranks:
        merged_results = merge_results(pair[0].gene_results, pair[1].gene_results)
        RankComparisonGenerator(merged_results).generate_gene_comparison_output(
            f"{pair[0].directory.name}__v__{pair[1].directory.name}"
        )


def generate_variant_rank_comparisons(comparison_ranks: [tuple]) -> None:
    """Generate the variant rank comparison of two result directories."""
    for pair in comparison_ranks:
        merged_results = merge_results(pair[0].variant_results, pair[1].variant_results)
        RankComparisonGenerator(merged_results).generate_variant_comparison_output(
            f"{pair[0].directory.name}__v__{pair[1].directory.name}"
        )


def benchmark_runs(
    results_directories: [TrackInputOutputDirectories],
    ranking_method: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Benchmark several result directories."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    gene_ranks_for_directories = []
    variant_ranks_for_directories = []
    for results_dir_and_input in results_directories:
        gene_rank_comparison, variant_rank_comparison = defaultdict(dict), defaultdict(dict)
        gene_ranks, variant_ranks = assess_prioritisation_for_results_directory(
            results_dir_and_input,
            ranking_method,
            threshold,
            gene_rank_comparison,
            variant_rank_comparison,
            gene_stats_writer,
            variants_stats_writer,
            gene_analysis,
            variant_analysis,
        )
        gene_ranks_for_directories.append(
            TrackGeneComparisons(results_dir_and_input.results_dir, gene_ranks)
        )
        variant_ranks_for_directories.append(
            TrackVariantComparisons(results_dir_and_input.results_dir, variant_ranks)
        )
    generate_gene_rank_comparisons(
        list(itertools.combinations(gene_ranks_for_directories, 2))
    ) if gene_analysis else None
    generate_variant_rank_comparisons(
        list(itertools.combinations(variant_ranks_for_directories, 2))
    ) if variant_analysis else None


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
    "--ranking-method",
    "-r",
    required=True,
    help="Ranking method for gene prioritisation.",
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
    default=True,
    required=False,
    type=bool,
    show_default=True,
    help="Analyse gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=True,
    required=False,
    type=bool,
    show_default=True,
    help="Analyse variant prioritisation",
)
def benchmark(
    directory: Path,
    phenopacket_dir: Path,
    ranking_method: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Benchmark the gene/variant prioritisation performance for a single run."""
    benchmark_directory(
        TrackInputOutputDirectories(results_dir=directory, phenopacket_dir=phenopacket_dir),
        ranking_method,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
    )


@click.command()
@click.option(
    "--directory1",
    "-d1",
    required=True,
    metavar="PATH",
    help="Baseline results directory for benchmarking, assumes contains subdirectories of pheval_gene_results/"
    "pheval_variant_results and the tool specific results directory.",
    type=Path,
)
@click.option(
    "--directory2",
    "-d2",
    required=True,
    metavar="PATH",
    help="Comparison results directory for benchmarking, assumes contains subdirectories of pheval_gene_results/"
    "pheval_variant_results and the tool specific results directory.",
    type=Path,
)
@click.option(
    "--phenopacket-dir1",
    "-p1",
    required=True,
    metavar="PATH",
    help="Full path to directory containing phenopackets for input for baseline directory.",
    type=Path,
)
@click.option(
    "--phenopacket-dir2",
    "-p2",
    required=True,
    metavar="PATH",
    help="Full path to directory containing phenopackets for input for comparison directory.",
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
    "--ranking-method",
    "-r",
    required=True,
    help="Ranking method for gene prioritisation.",
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
    default=True,
    required=False,
    type=bool,
    show_default=True,
    help="Analyse gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=True,
    required=False,
    type=bool,
    show_default=True,
    help="Analyse variant prioritisation",
)
def benchmark_comparison(
    directory1: Path,
    directory2: Path,
    phenopacket_dir1: Path,
    phenopacket_dir2: Path,
    ranking_method: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Benchmark the gene/variant prioritisation performance for two runs."""
    benchmark_runs(
        [
            TrackInputOutputDirectories(results_dir=directory1, phenopacket_dir=phenopacket_dir1),
            TrackInputOutputDirectories(results_dir=directory2, phenopacket_dir=phenopacket_dir2),
        ],
        ranking_method,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
    )
