# https://blog.actorsfit.com/a?ID=00700-7b9922d4-a0c0-4007-a62b-e15a7abcdc34

# -*- coding: cp936 -*-
from __future__ import generators
import jaydebeapi
import pandas as pd
from .decorators import *
import logging as log
from tqdm import tqdm
from .randomisation import *

log.basicConfig(
    level=log.INFO,
    format="%(asctime)-8s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d",
)


class DBConnector:
    def __init__(self, jar, driver, server, database, user, password):

        self.jar = jar
        self.driver = driver
        self.server = server
        self.database = database
        self.user = user
        self.password = password
        self.dbconn = None

    def create_connection(self):
        return jaydebeapi.connect(
            self.driver,
            f"{self.server}{self.database}",
            [self.user, self.password],
            self.jar,
        )

    def __enter__(self):
        self.dbconn = self.create_connection()
        return self.dbconn

    def __exit__(self, *other):
        self.dbconn.close()


class DBConnection:
    connection = None

    def __init__(cls, connection):
        DBConnection.connection = connection

    @classmethod
    def get_connection(cls):
        """Creates return new Singleton database connection"""
        return DBConnection.connection

    @classmethod
    def get_cursor(cls):
        connection = cls.get_connection()
        return connection.cursor()

    @classmethod
    def execute_query(cls, query, is_fetch=False):
        connection = cls.get_connection()
        result = None
        try:
            cursor = connection.cursor()
        except Exception:
            connection = cls.get_connection(new=True)  # Create new connection
            cursor = connection.cursor()
        cursor.execute(query)
        if is_fetch:
            result = cursor.fetchall()
        cursor.close()
        return result

    @classmethod
    def execute_many(cls, query, data):
        connection = cls.get_connection()
        try:
            cursor = connection.cursor()
        except Exception:
            connection = cls.get_connection(new=True)  # Create new connection
            cursor = connection.cursor()
        for d in data:
            try:
                cursor.execute(query, d)
            except Exception:
                pass
        cursor.close()


connector = DBConnector(
    jar="./lib/h2.jar",
    driver="org.h2.Driver",
    server="jdbc:h2:/home/data/exomiser-data/2202_phenotype/",
    user="sa",
    password="",
    database="2202_phenotype",
)


def ResultIterator(query: str, conn: DBConnection, arraysize: int = 1000) -> tuple:
    """
    Records rows generator
    Args:
        query (str): SQL query
        conn (DBConnection): Database connection
        arraysize (int, optional): chunk size Defaults to 1000.

    Returns:
        tuple: rows as tuple

    Yields:
        Iterator[tuple]: rows as tuple
    """
    cursor = conn.get_cursor()
    cursor.execute(query)
    while True:
        results = cursor.fetchmany(arraysize)
        if not results or len(results) == 0:
            break
        df = pd.DataFrame(results)
        df.columns = [i[0] for i in cursor.description]
        for result in df.itertuples():
            yield result


def create_table(conn: DBConnection, table_name: str):
    """
    Creating table
    Args:
        conn (DBConnection): Database connection
        table_name (str): table name
    """
    drop_query = f"DROP TABLE IF EXISTS EXOMISER.{table_name}_SCRAMBLE;"
    conn.execute_query(drop_query)
    create_query = f"""CREATE TABLE EXOMISER.{table_name}_SCRAMBLE (
      MAPPING_ID INTEGER,
      HP_ID VARCHAR(10),
      HP_TERM VARCHAR(200) NULL,
      ZP_ID VARCHAR(10),
      ZP_TERM VARCHAR(200),
      SIMJ DOUBLE NULL,
      IC DOUBLE NULL,
      SCORE DOUBLE NULL,
      LCS_ID VARCHAR(40),
      LCS_TERM VARCHAR(150)
    );"""
    conn.execute_query(create_query)


def count_total(conn: DBConnection, table_name: str) -> int:
    """count number of records from table {table_name}

    Args:
        conn (DBConnection): Database connection
        table_name (str): Table name

    Returns:
        int: number of records from table
    """
    count_query = f"SELECT COUNT(1) FROM EXOMISER.{table_name}"
    cnt_res = conn.execute_query(count_query, is_fetch=True)
    return cnt_res[0][0]


@measure_time
def insert(
    conn: DBConnection,
    table_name: str,
    scramble_factor: float,
    count: int,
    chunk: int = 10000,
):
    """insert modified lines in scrambled table

    Args:
        conn (DBConnection): Database connection
        table_name (str): table name
        scramble_factor (float): scramble factor
        count (int): total of rows from original table
        chunk (int, optional): chunk size. Defaults to 10000.
    """
    i = 1
    end = 0
    for j in tqdm(range(end, count // chunk)):
        start = (i - 1) * chunk
        end = (chunk * i) - 1
        data_to_update = select(conn, table_name, chunk, start)
        mod = ssp_randomisation(data_to_update, scramble_factor)
        sql = f"""INSERT INTO EXOMISER.{table_name}_SCRAMBLE (MAPPING_ID,HP_ID,HP_TERM,ZP_ID,ZP_TERM,SIMJ,IC,SCORE,LCS_ID,LCS_TERM) VALUES(?,?,?,?,?,?,?,?,?,?);"""
        mod.set_index("MAPPING_ID", inplace=True)
        try:
            conn.execute_many(sql, mod.values.tolist())
        except Exception as err:
            print(err)
        finally:
            i += 1


def select(
    conn: DBConnection, table_name: str, chunk: int, start: int = 0
) -> pd.DataFrame:
    """
    Select rows from original phenotypic database
    Args:
        conn (DBConnection): Database connection
        table_name (str): Original Table name
        chunk (int): Number of processed lines by execution
        start (int, optional): offset number

    Returns:
        pd.DataFrame: Original rows
    """
    select_query = f"SELECT * FROM EXOMISER.{table_name} LIMIT {chunk} OFFSET {start};"
    data_to_update = pd.DataFrame(ResultIterator(select_query, conn))
    return data_to_update


@memory(percentage=0.95)
def scramble_table(table_name: str, scramble_factor: float) -> None:
    with connector as cnn:
        conn = DBConnection(cnn)
        log.info(f"counting records from {table_name}")
        count_original = count_total(conn, table_name)
        log.info(f"Original records length: {count_original}")
        log.info(
            f"Scrambling records from table: {table_name} using {scramble_factor} magnitude"
        )
        insert(conn, table_name, scramble_factor, count_original, 3000)
        count_scramble = count_total(conn, f"{table_name}_SCRAMBLE")
        log.info(f"Scrambled records length: {count_scramble}")
        log.info("Done")
