import ast
import re
from typing import List, Type, Union

import duckdb
from duckdb import DuckDBPyConnection

from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)


class BenchmarkDBManager:
    """
    Class to connect to database.
    """

    def __init__(self, benchmark_name: str):
        """Initialise the BenchmarkDBManager class."""
        self.conn = self.get_connection(
            f"{benchmark_name}" if str(benchmark_name).endswith(".db") else f"{benchmark_name}.db"
        )

    def initialise(self):
        """Initialise the duckdb connection."""
        self.add_contains_function()

    @staticmethod
    def get_connection(db_name: str) -> DuckDBPyConnection:
        """
        Get a connection to the database.
        Returns:
            DuckDBPyConnection: Connection to the database.
        """
        conn = duckdb.connect(db_name)
        return conn

    def add_column_integer_default(self, table_name: str, column: str, default: int = 0) -> None:
        """
        Add a column to an existing table with an integer default value.
        Args:
            table_name (str): Name of the table.
            column (str): Name of the column to add.
            default (int): Default integer value to add.
        """
        try:
            self.conn.execute(
                f'ALTER TABLE {table_name} ADD COLUMN "{column}" INTEGER DEFAULT {default}'
            )
            self.conn.execute(f'UPDATE {table_name} SET "{column}" = ?', (default,))
            self.conn.commit()
        except duckdb.CatalogException:
            pass

    def drop_table(self, table_name: str) -> None:
        """
        Drop a table from the database.
        Args:
            table_name: Name of the table to drop from the database
        """
        self.conn.execute(f"""DROP TABLE IF EXISTS "{table_name}";""")

    @staticmethod
    def contains_entity_function(entity: str, known_causative_entity: str) -> bool:
        """
        Determines if a known causative entity is present within an entity or list of entities.
        Args:
            entity (str): The entity to be checked. It can be a single entity or a string representation of a list.
            known_causative_entity (str): The entity to search for within the `entity`.

        Returns:
            bool: `True` if `known_causative_entity` is found in `entity` (or its list representation),
                `False` otherwise.
        """
        list_pattern = re.compile(r"^\[\s*(?:[^\[\],\s]+(?:\s*,\s*[^\[\],\s]+)*)?\s*]$")
        entity = entity.replace("nan", "None").replace("NaN", "None")
        if list_pattern.match(str(entity)):
            list_representation = ast.literal_eval(entity)
            if isinstance(list_representation, list):
                return known_causative_entity in list_representation
        return known_causative_entity == entity

    def add_contains_function(self) -> None:
        """
        Adds a custom `contains_entity_function` to the DuckDB connection if it does not already exist.
        """
        result = self.conn.execute(
            "SELECT * FROM duckdb_functions() WHERE function_name = ?", ["contains_entity_function"]
        ).fetchall()
        if not result:
            self.conn.create_function("contains_entity_function", self.contains_entity_function)

    def parse_table_into_dataclass(
        self,
        table_name: str,
        dataclass: Union[
            Type[RankedPhEvalGeneResult],
            Type[RankedPhEvalVariantResult],
            Type[RankedPhEvalDiseaseResult],
        ],
    ) -> Union[
        List[RankedPhEvalGeneResult],
        List[RankedPhEvalVariantResult],
        List[RankedPhEvalDiseaseResult],
    ]:
        """
        Parses a DuckDB table into a list of dataclass instances.
        Args:
            table_name (str): The name of the DuckDB table to be parsed.
            dataclass (Union[Type[RankedPhEvalGeneResult], Type[RankedPhEvalVariantResult],
            Type[RankedPhEvalDiseaseResult]]):
                The dataclass type to which each row in the table should be mapped.

        Returns:
            List[dataclass]: A list of instances of the provided dataclass, each representing a row from the table.
        """
        result = (
            self.conn.execute(f"SELECT * FROM '{table_name}'").fetchdf().to_dict(orient="records")
        )
        return [dataclass(**row) for row in result]

    def check_table_exists(self, table_name: str) -> bool:
        """
        Check if a table exists in the connected DuckDB database.
        Args:
            table_name (str): The name of the table to check for existence.
        Returns:
            bool: Returns `True` if the table exists in the database, `False` otherwise.
        """
        result = self.conn.execute(
            f"SELECT * FROM information_schema.tables WHERE table_name = '{table_name}'"
        ).fetchall()
        if result:
            return True
        return False

    def close(self):
        """Close the connection to the database."""
        self.conn.close()
