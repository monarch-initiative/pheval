from pathlib import Path
from typing import List

import polars as pl

from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    PhenopacketUtil,
    ProbandCausativeGene,
    ProbandDisease,
    phenopacket_reader,
)


def calculate_end_pos(variant_start: int, variant_ref: str) -> int:
    """Calculate the end position for a variant
    Args:
        variant_start (int): The start position of the variant
        variant_ref (str): The reference allele of the variant

    Returns:
        int: The end position of the variant
    """
    return variant_start + len(variant_ref) - 1


class PhenopacketTruthSet:
    """Class for finding the causative gene/disease/variant from a phenopacket"""

    def __init__(self, phenopacket_dir: Path):
        self.phenopacket_dir = phenopacket_dir

    def _get_phenopacket_path(self, phenopacket_name: str) -> Path:
        """
        Get the phenopacket path for a given phenopacket name.
        Args:
            phenopacket_name (str): Name of the phenopacket.
        Returns:
            Path: Path to the phenopacket path.
        """
        phenopacket_path = self.phenopacket_dir.joinpath(f"{phenopacket_name}.json")
        if not phenopacket_path.exists():
            raise FileNotFoundError(phenopacket_name + " not found in corpus!")
        return phenopacket_path

    def _get_phenopacket_util(self, phenopacket_name: str) -> PhenopacketUtil:
        """
        Get the phenopacket util for a given phenopacket name.
        Args:
            phenopacket_name (str): Name of the phenopacket.
        Returns:
            PhenopacketUtil: PhenopacketUtil object.
        """
        phenopacket_path = self._get_phenopacket_path(phenopacket_name)
        phenopacket = phenopacket_reader(phenopacket_path)
        return PhenopacketUtil(phenopacket)

    def _get_causative_genes(self, phenopacket_name: str) -> List[ProbandCausativeGene]:
        """
        Get the causative genes for a given phenopacket.
        Args:
            phenopacket_name (str): Name of the phenopacket.
        Returns:
            List[ProbandCausativeGene]: List of ProbandCausativeGene.
        """
        phenopacket_util = self._get_phenopacket_util(phenopacket_name)
        return phenopacket_util.diagnosed_genes()

    def _get_causative_variants(self, phenopacket_name: str) -> List[GenomicVariant]:
        """
        Get the causative variants for a given phenopacket.
        Args:
            phenopacket_name (str): Name of the phenopacket.
        Returns:
            List[GenomicVariant]: List of GenomicVariant.
        """
        phenopacket_util = self._get_phenopacket_util(phenopacket_name)
        return phenopacket_util.diagnosed_variants()

    def _get_causative_diseases(self, phenopacket_name: str) -> List[ProbandDisease]:
        """
        Get the diseases for a given phenopacket.
        Args:
            phenopacket_name (str): Name of the phenopacket.
        Returns:
            List[ProbandDisease]: List of ProbandDisease
        """
        phenopacket_util = self._get_phenopacket_util(phenopacket_name)
        return phenopacket_util.diagnoses()

    def classified_gene(self, result_name: str) -> pl.DataFrame:
        """
        Classify gene results for a given phenopacket.
        Args:
            result_name (str): Name of the result file.
        Returns:
            pl.DataFrame: Classified ranked gene results.
        """
        causative_genes = self._get_causative_genes(result_name)
        gene_symbols = [causative_gene.gene_symbol for causative_gene in causative_genes]
        gene_identifiers = [causative_gene.gene_identifier for causative_gene in causative_genes]
        return pl.DataFrame(
            {
                "gene_symbol": [g for g in gene_symbols],
                "gene_identifier": [g for g in gene_identifiers],
            }
        ).with_columns(
            [
                pl.lit(0).cast(pl.Float64).alias("score"),
                pl.lit(0).cast(pl.Int64).alias("rank"),
                pl.lit(True).alias("true_positive"),
            ]
        )

    @staticmethod
    def merge_gene_results(ranked_results: pl.DataFrame, output_file: Path) -> pl.DataFrame:
        """
        Merge ranked gene results with the classified genes.
        Args:
            ranked_results (pl.DataFrame): Ranked gene results.
            output_file (Path): Path to the output file.
        Returns:
            pl.DataFrame: Merged ranked gene results.
        """
        classified_results = pl.read_parquet(output_file)
        return (
            ranked_results.with_columns(
                (
                    pl.col("gene_symbol").is_in(classified_results["gene_symbol"])
                    | pl.col("gene_identifier").is_in(classified_results["gene_identifier"])
                ).alias("true_positive")
            )
            .with_columns(pl.col("rank").cast(pl.Int64))
            .select(classified_results.columns)
            .vstack(
                classified_results.filter(
                    ~pl.col("gene_symbol").is_in(ranked_results["gene_symbol"])
                )
            )
        )

    def classified_variant(self, result_name: str) -> pl.DataFrame:
        """
        Classified variant results for a given phenopacket.
        Args:
            result_name (str): Name of the result file.
        Returns:
            pl.DataFrame: Classified ranked variant results.
        """
        variants = self._get_causative_variants(result_name)
        return pl.DataFrame(
            {
                "chrom": [v.chrom for v in variants],
                "start": [v.pos for v in variants],
                "end": [calculate_end_pos(v.pos, v.ref) for v in variants],
                "ref": [v.ref for v in variants],
                "alt": [v.alt for v in variants],
            }
        ).with_columns(
            [
                pl.concat_str(["chrom", "start", "ref", "alt"], separator="-").alias("variant_id"),
                pl.lit(0.0).cast(pl.Float64).alias("score"),
                pl.lit(0).cast(pl.Int64).alias("rank"),
                pl.lit(True).alias("true_positive"),
            ]
        )

    @staticmethod
    def merge_variant_results(ranked_results: pl.DataFrame, output_file: Path) -> pl.DataFrame:
        """
        Merge ranked variant results with the classified variants.
        Args:
            ranked_results (pl.DataFrame): Ranked variant results.
            output_file (Path): Path to the output file.
        Returns:
            pl.DataFrame: Merged ranked variant results.
        """
        classified_results = pl.read_parquet(output_file)
        return (
            ranked_results.with_columns(
                [
                    pl.struct(["chrom", "start", "end", "ref", "alt"])
                    .is_in(
                        classified_results.select(
                            pl.struct(["chrom", "start", "end", "ref", "alt"])
                        ).to_series()
                    )
                    .alias("true_positive")
                ]
            )
            .with_columns(pl.col("rank").cast(pl.Int64))
            .select(classified_results.columns)
            .vstack(
                classified_results.filter(
                    ~pl.struct(["chrom", "start", "end", "ref", "alt"]).is_in(
                        ranked_results.select(
                            pl.struct(["chrom", "start", "end", "ref", "alt"])
                        ).to_series()
                    )
                )
            )
        )

    def classified_disease(self, result_name: str) -> pl.DataFrame:
        """
        Classify disease results for a given phenopacket.
        Args:
            result_name (str): Name of the result file.
        Returns:
            pl.DataFrame: Classified ranked disease results.
        """
        diseases = self._get_causative_diseases(result_name)
        disease_identifiers = list(set(disease.disease_identifier for disease in diseases))
        return pl.DataFrame(
            {
                "disease_identifier": [d for d in disease_identifiers],
            }
        ).with_columns(
            [
                pl.lit(0).cast(pl.Float64).alias("score"),
                pl.lit(0).cast(pl.Int64).alias("rank"),
                pl.lit(True).alias("true_positive"),
            ]
        )

    @staticmethod
    def merge_disease_results(ranked_results: pl.DataFrame, output_file: Path) -> pl.DataFrame:
        """
        Merge ranked disease results with the classified diseases.
        Args:
            ranked_results (pl.DataFrame): Ranked disease results.
            output_file (Path): Path to the output file.
        Returns:
            pl.DataFrame: Merged ranked disease results.
        """
        classified_results = pl.read_parquet(output_file)
        return (
            ranked_results.with_columns(
                (
                    pl.col("disease_identifier").is_in(classified_results["disease_identifier"])
                ).alias("true_positive")
            )
            .with_columns(pl.col("rank").cast(pl.Int64))
            .select(classified_results.columns)
            .vstack(
                classified_results.filter(
                    ~pl.col("disease_identifier").is_in(ranked_results["disease_identifier"])
                )
            )
        )
