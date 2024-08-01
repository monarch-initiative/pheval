import duckdb
from duckdb import DuckDBPyConnection


class DBConnector:

    def __init__(self):
        self.conn = self.get_connection()

    @staticmethod
    def get_connection() -> DuckDBPyConnection:
        conn = duckdb.connect("analysis.db")
        return conn

    def add_column_integer_default(self, table_name: str, column: str, default: int = 0) -> None:
        try:
            self.conn.execute(
                f'ALTER TABLE {table_name} ADD COLUMN "{column}" INTEGER DEFAULT {default}'
            )
            self.conn.execute(f'UPDATE {table_name} SET "{column}" = ?', (default,))
            self.conn.commit()
        except duckdb.CatalogException:
            pass

    def drop_table(self, table_name: str) -> None:
        self.conn.execute(f"""DROP TABLE IF EXISTS "{table_name}";""")

    def close(self):
        self.conn.close()
