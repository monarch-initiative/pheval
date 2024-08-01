import duckdb
from duckdb import DuckDBPyConnection


class DBConnector:
    """
    Class to connect to database.
    """

    def __init__(self):
        """Initialize the DBConnector class."""
        self.conn = self.get_connection()

    @staticmethod
    def get_connection() -> DuckDBPyConnection:
        """
        Get a connection to the database.
        Returns:
            DuckDBPyConnection: Connection to the database.
        """
        conn = duckdb.connect("analysis.db")
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

    def close(self):
        """Close the connection to the database."""
        self.conn.close()
