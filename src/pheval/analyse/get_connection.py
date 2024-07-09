import duckdb
from duckdb import DuckDBPyConnection


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def get_connection() -> DuckDBPyConnection:
    """
    Get a connection to the benchmarking results database.
    """
    conn = duckdb.connect("analysis.db")
    return conn
