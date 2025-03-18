"""PhEval utils Command Line Interface"""

from pathlib import Path
from typing import List

import click

from pheval.analyse.benchmark import benchmark_runs
from pheval.analyse.generate_plots import generate_plots_from_db
from pheval.prepare.create_noisy_phenopackets import scramble_phenopackets
from pheval.prepare.create_spiked_vcf import spike_vcfs
from pheval.prepare.custom_exceptions import InputError, MutuallyExclusiveOptionError
from pheval.prepare.prepare_corpus import prepare_corpus
from pheval.prepare.update_phenopacket import update_phenopackets
from pheval.utils.exomiser import semsim_to_exomiserdb
from pheval.utils.semsim_utils import percentage_diff, semsim_heatmap_plot
from pheval.utils.utils import semsim_scramble


@click.command("semsim-scramble")
@click.option(
    "--input",
    "-i",
    required=True,
    metavar="FILE",
    help="Path to the semantic similarity profile to be scrambled.",
    type=Path,
)
@click.option(
    "--output",
    "-o",
    metavar="FILE",
    required=True,
    help="Path where the scrambled semsim file will be written.",
    type=Path,
)
@click.option(
    "--score-column",
    "-c",
    required=True,
    multiple=True,
    type=click.Choice(
        ["jaccard_similarity", "dice_similarity", "phenodigm_score"], case_sensitive=False
    ),
    help="Score column that will be scrambled",
)
@click.option(
    "--scramble-factor",
    "-s",
    metavar=float,
    default=0.5,
    show_default=True,
    type=float,
    help="""Scramble Magnitude (noise)
    that will be applied to semantic similarity score column (e.g. jaccard similarity).""",
)
def semsim_scramble_command(
    input: Path, output: Path, score_column: List[str], scramble_factor: float
):
    """Scrambles semsim profile multiplying score value by scramble factor
    Args:
        input (Path): Path file that points out to the semsim profile
        output (Path): Path file that points out to the output file
        score_column (List[str]): Score column(s) that will be scrambled
        scramble_factor (float): Scramble Magnitude
    """
    semsim_scramble(input, output, score_column, scramble_factor)


@click.command("scramble-phenopackets")
@click.option(
    "--phenopacket-path",
    "-p",
    metavar="PATH",
    help="Path to phenopacket.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_dir"],
)
@click.option(
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    help="Path to phenopackets directory.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_path"],
)
@click.option(
    "--scramble-factor",
    "-s",
    metavar=float,
    required=True,
    default=0.5,
    show_default=True,
    help="Scramble factor for randomising phenopacket phenotypic profiles.",
    type=float,
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="noisy_phenopackets",
    type=Path,
)
@click.option(
    "--local-ontology-cache",
    "-l",
    metavar="PATH",
    required=False,
    help="Path to the local ontology cache, e.g., path to the hp.obo.",
    default=None,
    type=Path,
)
def scramble_phenopackets_command(
    phenopacket_path: Path,
    phenopacket_dir: Path,
    scramble_factor: float,
    output_dir: Path,
    local_ontology_cache: Path,
):
    """Generate noisy phenopackets from existing ones."""
    if phenopacket_path is None and phenopacket_dir is None:
        raise InputError("Either a phenopacket or phenopacket directory must be specified")
    else:
        scramble_phenopackets(
            output_dir, phenopacket_path, phenopacket_dir, scramble_factor, local_ontology_cache
        )


@click.command("semsim-comparison")
@click.option(
    "--semsim-left",
    "-L",
    required=True,
    metavar="FILE",
    help="Path to the first semantic similarity profile.",
)
@click.option(
    "--semsim-right",
    "-R",
    required=True,
    metavar="FILE",
    help="Path to the second semantic similarity profile.",
)
@click.option(
    "--score-column",
    "-c",
    required=True,
    type=click.Choice(
        ["jaccard_similarity", "dice_similarity", "phenodigm_score"], case_sensitive=False
    ),
    help="Score column that will be used in comparison",
)
@click.option(
    "--analysis",
    "-a",
    required=True,
    type=click.Choice(["heatmap", "percentage_diff"], case_sensitive=False),
    help="""There are two types of analysis:
        heatmap - Generates a heatmap plot that shows the differences between the semantic similarity profiles using the
        score column for this purpose. Defaults to "heatmap".
        percentage_diff - Calculates the score column percentage difference between the semantic similarity profiles""",
)
@click.option(
    "--output",
    "-o",
    metavar="FILE",
    default="percentage_diff.semsim.tsv",
    help="Output path for the difference tsv. Defaults to percentage_diff.semsim.tsv",
)
def semsim_comparison(
    semsim_left: Path,
    semsim_right: Path,
    score_column: str,
    analysis: str,
    output: Path = "percentage_diff.semsim.tsv",
):
    """Compares two semantic similarity profiles

    Args:
        semsim-left (Path): File path of the first semantic similarity profile
        semsim-right (Path): File path of the second semantic similarity profile
        output (Path): Output path for the difference tsv. Defaults to "percentage_diff.semsim.tsv".
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        analysis (str): There are two types of analysis:
        heatmap - Generates a heatmap plot that shows the differences between the semantic similarity profiles using the
        score column for this purpose. Defaults to "heatmap".
        percentage_diff - Calculates the score column percentage difference between the semantic similarity profiles.
    """
    if analysis == "heatmap":
        return semsim_heatmap_plot(semsim_left, semsim_right, score_column)
    if analysis == "percentage_diff":
        percentage_diff(semsim_left, semsim_right, score_column, output)


@click.command("update-phenopackets")
@click.option(
    "--phenopacket-path",
    "-p",
    metavar="PATH",
    help="Path to phenopacket.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_dir"],
)
@click.option(
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    help="Path to phenopacket directory for updating.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_path"],
)
@click.option(
    "--output-dir",
    "-o",
    metavar="PATH",
    required=True,
    help="Path to write phenopacket.",
    type=Path,
)
@click.option(
    "--gene-identifier",
    "-g",
    required=False,
    default="ensembl_id",
    show_default=True,
    help="Gene identifier to add to phenopacket",
    type=click.Choice(["ensembl_id", "entrez_id", "hgnc_id"]),
)
def update_phenopackets_command(
    phenopacket_path: Path, phenopacket_dir: Path, output_dir: Path, gene_identifier: str
):
    """Update gene symbols and identifiers for phenopackets."""
    if phenopacket_path is None and phenopacket_dir is None:
        raise InputError("Either a phenopacket or phenopacket directory must be specified")
    update_phenopackets(gene_identifier, phenopacket_path, phenopacket_dir, output_dir)


@click.command("create-spiked-vcfs")
@click.option(
    "--phenopacket-path",
    "-p",
    metavar="PATH",
    help="Path to phenopacket.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_dir"],
)
@click.option(
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    help="Path to phenopacket directory for updating.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_path"],
)
@click.option(
    "--hg19-template-vcf",
    "-hg19",
    metavar="PATH",
    required=False,
    help="Template hg19 VCF file",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg19_vcf_dir"],
)
@click.option(
    "--hg38-template-vcf",
    "-hg38",
    metavar="PATH",
    required=False,
    help="Template hg38 VCF file",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg38_vcf_dir"],
)
@click.option(
    "--hg19-vcf-dir",
    "-hg19-dir",
    metavar="PATH",
    required=False,
    help="Path to directory containing hg19 VCF templates.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg19_template_vcf"],
)
@click.option(
    "--hg38-vcf-dir",
    "-hg38-dir",
    metavar="PATH",
    required=False,
    help="Path to directory containing hg38 VCF templates.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg38_template_vcf"],
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="vcf",
    type=Path,
)
def create_spiked_vcfs_command(
    phenopacket_path: Path,
    phenopacket_dir: Path,
    output_dir: Path,
    hg19_template_vcf: Path = None,
    hg38_template_vcf: Path = None,
    hg19_vcf_dir: Path = None,
    hg38_vcf_dir: Path = None,
):
    """
    Create spiked VCF from either a Phenopacket or a Phenopacket directory.

    Args:
        phenopacket_path (Path): Path to a single Phenopacket file (optional).
        phenopacket_dir (Path): Path to a directory containing Phenopacket files (optional).
        output_dir (Path): The directory to store the generated spiked VCF file(s).
        hg19_template_vcf (Path): Path to the hg19 template VCF file (optional).
        hg38_template_vcf (Path): Path to the hg38 template VCF file (optional).
        hg19_vcf_dir (Path): Path to the directory containing the hg19 VCF files (optional).
        hg38_vcf_dir (Path): Path to the directory containing the hg38 VCF files (optional).
    """
    if phenopacket_path is None and phenopacket_dir is None:
        raise InputError("Either a phenopacket or phenopacket directory must be specified")
    spike_vcfs(
        output_dir,
        phenopacket_path,
        phenopacket_dir,
        hg19_template_vcf,
        hg38_template_vcf,
        hg19_vcf_dir,
        hg38_vcf_dir,
    )


@click.command()
@click.option(
    "--run-yaml",
    "-r",
    required=True,
    metavar="PATH",
    help="Path to yaml configuration file for benchmarking.",
    type=Path,
)
def benchmark(
    run_yaml: Path,
):
    """Benchmark the gene/variant/disease prioritisation performance for runs."""
    benchmark_runs(
        run_yaml,
    )


@click.command("semsim-to-exomiserdb")
@click.option(
    "--input-file",
    "-i",
    required=True,
    metavar="FILE",
    help="Semsim input file.",
    type=Path,
)
@click.option(
    "--object-prefix",
    required=True,
    metavar="object-prefix",
    help="Object Prefix. e.g. MP",
    type=str,
)
@click.option(
    "--subject-prefix",
    required=True,
    metavar="subject-prefix",
    help="Subject Prefix. e.g. HP",
    type=str,
)
@click.option(
    "--db-path",
    "-d",
    required=True,
    metavar="db-path",
    help="""Exomiser Phenotypic Database Folder Path.
    (e.g. /exomiser_folder/2209_phenotype/2209_phenotype/).
    This is the path where the phenotypic database folder will be written out.""",
    type=Path,
)
def semsim_to_exomiserdb_command(
    input_file: Path, object_prefix: str, subject_prefix: str, db_path: Path
):
    """ingests semsim file into exomiser phenotypic database

    Args:
        input_file (Path): semsim input file. e.g phenio-plus-hp-mp.0.semsimian.tsv
        object_prefix (str): object prefix. e.g. MP
        subject_prefix (str): subject prefix e.g HP
        db_path (Path): Exomiser Phenotypic Database Folder Path. (e.g. /exomiser_folder/2209_phenotype/2209_phenotype/)
    """
    semsim_to_exomiserdb(input_file, object_prefix, subject_prefix, db_path)


@click.command()
@click.option(
    "--benchmark-db",
    "-b",
    required=True,
    metavar="PATH",
    help="Path to benchmark db output by PhEval benchmark commands.",
    type=Path,
)
@click.option(
    "--run-data",
    "-r",
    required=True,
    metavar="PATH",
    help="Path to yaml configuration file for benchmarking.",
    type=Path,
)
def generate_plots(
    benchmark_db: Path,
    run_data: Path,
):
    """Generate bar plot from benchmark db."""
    generate_plots_from_db(benchmark_db, run_data)


@click.command("prepare-corpus")
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    metavar="PATH",
    help="Path to phenopacket corpus directory..",
    type=Path,
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify whether to check for complete variant records in the phenopackets.",
)
@click.option(
    "--gene-analysis/--no-gene-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify whether to check for complete gene records in the phenopackets.",
)
@click.option(
    "--disease-analysis/--no-disease-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify whether to check for complete disease records in the phenopackets.",
)
@click.option(
    "--gene-identifier",
    "-g",
    required=False,
    help="Gene identifier to update in phenopacket",
    type=click.Choice(["ensembl_id", "entrez_id", "hgnc_id"]),
)
@click.option(
    "--hg19-template-vcf",
    "-hg19",
    metavar="PATH",
    required=False,
    help="Template hg19 VCF file",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg19_vcf_dir"],
)
@click.option(
    "--hg38-template-vcf",
    "-hg38",
    metavar="PATH",
    required=False,
    help="Template hg38 VCF file",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg38_vcf_dir"],
)
@click.option(
    "--hg19-vcf-dir",
    "-hg19-dir",
    metavar="PATH",
    required=False,
    help="Path to directory containing hg19 VCF templates.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg19_template_vcf"],
)
@click.option(
    "--hg38-vcf-dir",
    "-hg38-dir",
    metavar="PATH",
    required=False,
    help="Path to directory containing hg38 VCF templates.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["hg38_template_vcf"],
)
@click.option(
    "--output-dir",
    "-o",
    metavar="PATH",
    required=True,
    help="Path to output prepared corpus.",
    default="prepared_corpus",
    type=Path,
)
def prepare_corpus_command(
    phenopacket_dir: Path,
    variant_analysis: bool,
    gene_analysis: bool,
    disease_analysis: bool,
    gene_identifier: str,
    hg19_template_vcf: Path,
    hg38_template_vcf: Path,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
    output_dir: Path,
):
    """
    Prepare a corpus of Phenopackets for analysis, optionally checking for complete variant records and updating
    gene identifiers.

    Args:
        phenopacket_dir (Path): The path to the directory containing Phenopackets.
        variant_analysis (bool): If True, check for complete variant records in the Phenopackets.
        gene_analysis (bool): If True, check for complete gene records in the Phenopackets.
        disease_analysis (bool): If True, check for complete disease records in the Phenopackets.
        gene_identifier (str): Identifier for updating gene identifiers, if applicable.
        hg19_template_vcf (Path): Path to the hg19 template VCF file (optional).
        hg38_template_vcf (Path): Path to the hg38 template VCF file (optional).
        hg19_vcf_dir (Path): Path to the directory containing the hg19 VCF files (optional).
        hg38_vcf_dir (Path): Path to the directory containing the hg38 VCF files (optional).
        output_dir (Path): The directory to save the prepared Phenopackets and, optionally, VCF files.
    Notes:
        To spike variants into VCFs for variant-based analysis at least one of hg19_template_vcf, hg38_template_vcf,
        hg19_vcf_dir or hg38_vcf_dir is required.
    """
    prepare_corpus(
        phenopacket_dir,
        variant_analysis,
        gene_analysis,
        disease_analysis,
        gene_identifier,
        hg19_template_vcf,
        hg38_template_vcf,
        hg19_vcf_dir,
        hg38_vcf_dir,
        output_dir,
    )
