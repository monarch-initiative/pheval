from pathlib import Path
from typing import List

from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.benchmark_generator import (
    BenchmarkRunOutputGenerator,
    DiseaseBenchmarkRunOutputGenerator,
    GeneBenchmarkRunOutputGenerator,
    VariantBenchmarkRunOutputGenerator,
)
from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    PhenopacketUtil,
    ProbandCausativeGene,
    ProbandDisease,
    phenopacket_reader,
)


def _obtain_causative_diseases(phenopacket_path: Path) -> List[ProbandDisease]:
    """
    Obtain known diseases from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.

    Returns:
       List[ProbandDisease]: A list of known diseases associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnoses()


def _obtain_causative_variants(phenopacket_path: Path) -> List[GenomicVariant]:
    """
    Obtain known variants from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.

    Returns:
       List[GenomicVariant]: A list of known variants associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_variants()


def _obtain_causative_genes(phenopacket_path: Path) -> List[ProbandCausativeGene]:
    """
    Obtain known genes from a Phenopacket.
    Args:
       phenopacket_path (Path): Path to the Phenopacket file.
    Returns:
       List[ProbandCausativeGene]: A list of known genes associated with the proband,
       extracted from the Phenopacket.
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_util = PhenopacketUtil(phenopacket)
    return phenopacket_util.diagnosed_genes()


class CorpusParser:
    """Class for parsing phenopacket corpus and retrieving known variants/genes/diseases."""

    def __init__(self, benchmark_name: str, phenopacket_dir: Path) -> None:
        """
        Initialise the CorpusParser class.
        Args:
            phenopacket_dir (Path): Path to the Phenopacket directory.
        """
        self.phenopacket_dir = phenopacket_dir
        self.conn = BenchmarkDBManager(benchmark_name).conn
        self.table_name = phenopacket_dir.parents[0].name

    def _create_gene_table(self) -> None:
        """
        Create the Gene benchmarking table if it doesn't already exist.
        """
        self.conn.execute(
            f"""
                    CREATE TABLE IF NOT EXISTS {self.table_name}_gene (
                        identifier VARCHAR(255) PRIMARY KEY,
                        phenopacket VARCHAR,
                        gene_symbol VARCHAR,
                        gene_identifier VARCHAR
                    )
                    """
        )

    def _create_variant_table(self) -> None:
        """
        Create the Variant benchmarking table if it doesn't already exist.
        """
        self.conn.execute(
            f"""
                    CREATE TABLE IF NOT EXISTS {self.table_name}_variant (
                        identifier VARCHAR(255) PRIMARY KEY,
                        phenopacket VARCHAR,
                        chrom VARCHAR,
                        pos INTEGER,
                        "ref" VARCHAR,
                        alt VARCHAR
                    )
                    """
        )

    def _create_disease_table(self):
        """
        Create the Disease benchmarking table if it doesn't already exist.
        """
        self.conn.execute(
            f"""
                    CREATE TABLE IF NOT EXISTS {self.table_name}_disease (
                        identifier VARCHAR(255) PRIMARY KEY,
                        phenopacket VARCHAR,
                        disease_identifier VARCHAR,
                        disease_name VARCHAR
                    )
                    """
        )

    def _create_tables(self, benchmark_generator: BenchmarkRunOutputGenerator) -> None:
        """
        Create tables based on the benchmarking analysis specified.
        Args:
            benchmark_generator (BenchmarkRunOutputGenerator): Class instance of the benchmark generator type.
        """

        if isinstance(benchmark_generator, GeneBenchmarkRunOutputGenerator):
            self._create_gene_table()
        if isinstance(benchmark_generator, VariantBenchmarkRunOutputGenerator):
            self._create_variant_table()
        if isinstance(benchmark_generator, DiseaseBenchmarkRunOutputGenerator):
            self._create_disease_table()

    def _insert_genes(self, phenopacket_path: Path, genes: List[ProbandCausativeGene]) -> None:
        """
        Insert known disease-causing genes into the Gene benchmarking table.
        Args:
            phenopacket_path(Path): Path to the Phenopacket file.
            genes(List[ProbandCausativeGene]): List of known genes associated with the proband.
        """
        for gene in genes:
            identifier = f"{phenopacket_path.name}-{gene.gene_symbol}"
            self.conn.execute(
                f"""
                INSERT OR IGNORE INTO {self.table_name}_gene (identifier, phenopacket, gene_symbol, gene_identifier)
                VALUES (?, ?, ?, ?)
                """,
                (identifier, phenopacket_path.name, gene.gene_symbol, gene.gene_identifier),
            )

    def _insert_variants(self, phenopacket_path: Path, variants: List[GenomicVariant]) -> None:
        """
        Insert known variants into the Variant benchmarking table.
        Args:
            phenopacket_path (Path): Path to the Phenopacket file.:
            variants (List[GenomicVariant]): List of known variants associated with the proband.
        """
        for variant in variants:
            identifier = (
                f"{phenopacket_path.name}-{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}"
            )
            self.conn.execute(
                f"""
                INSERT OR IGNORE INTO {self.table_name}_variant (identifier, phenopacket, chrom, pos, "ref", alt)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    identifier,
                    phenopacket_path.name,
                    variant.chrom,
                    variant.pos,
                    variant.ref,
                    variant.alt,
                ),
            )

    def _insert_diseases(self, phenopacket_path: Path, diseases: List[ProbandDisease]) -> None:
        """
        Insert known diseases into the Disease benchmarking table.
        Args:
            phenopacket_path (Path): Path to the Phenopacket file.:
            diseases (List[ProbandDisease]): List of known diseases associated with the proband.
        """
        for disease in diseases:
            identifier = f"{phenopacket_path.name}-{disease.disease_identifier}"
            self.conn.execute(
                f"INSERT OR IGNORE INTO {self.table_name}_disease "
                f"(identifier, phenopacket, disease_identifier, disease_name) VALUES (?, ?, ?, ?)",
                (
                    identifier,
                    phenopacket_path.name,
                    disease.disease_identifier,
                    disease.disease_name,
                ),
            )

    def parse_corpus(self, benchmark_generator: BenchmarkRunOutputGenerator) -> None:
        """
        Parse the phenopacket corpus and add known genes/variants/diseases to relevant benchmarking tables.
        Args:
            benchmark_generator (BenchmarkRunOutputGenerator): Class instance of the benchmark generator type.
        """
        self._create_tables(benchmark_generator)
        for phenopacket_path in all_files(self.phenopacket_dir):
            if isinstance(benchmark_generator, GeneBenchmarkRunOutputGenerator):
                genes = _obtain_causative_genes(phenopacket_path)
                self._insert_genes(phenopacket_path, genes)
            if isinstance(benchmark_generator, VariantBenchmarkRunOutputGenerator):
                variants = _obtain_causative_variants(phenopacket_path)
                self._insert_variants(phenopacket_path, variants)
            if isinstance(benchmark_generator, DiseaseBenchmarkRunOutputGenerator):
                diseases = _obtain_causative_diseases(phenopacket_path)
                self._insert_diseases(phenopacket_path, diseases)
        self.conn.close()
