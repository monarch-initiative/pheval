import polars as pl
from duckdb import DuckDBPyConnection

from pheval.utils.logger import get_logger

logger = get_logger()


def load_table_lazy(table_name: str, conn: DuckDBPyConnection) -> pl.LazyFrame:
    logger.info(f"Loading table {table_name}")
    return pl.from_arrow(conn.execute(f"SELECT * FROM {table_name}").fetch_arrow_table()).lazy()


def write_table(conn: DuckDBPyConnection, df: pl.DataFrame, table_name: str) -> None:
    """
    Write table to DuckDB database.
    Args:
        conn (DuckDBPyConnection): DuckDB connection.
        df (pl.DataFrame): DuckDB dataframe.
        table_name (str): Table name.
    """
    logger.info(f"Storing results in {table_name}.")
    conn.execute(f"""CREATE TABLE "{table_name}" AS SELECT * FROM df""")
