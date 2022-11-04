# https://blog.actorsfit.com/a?ID=00700-7b9922d4-a0c0-4007-a62b-e15a7abcdc34

# -*- coding: cp936 -*-
from __future__ import generators
import jaydebeapi
import pandas as pd
import pheval.randomisation as randomisation
from pheval.decorators import *
import logging as log
from tqdm import tqdm
import os
from pheval.memory import *
import subprocess
from pheval.mappings import *

info_log = log.getLogger("info")
info_debug = log.getLogger("debug")


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

    def close(cls):
        return cls.connection.close()

    @classmethod
    def get_cursor(cls):
        connection = cls.get_connection()
        return connection.cursor()

    @classmethod
    def execute_query(cls, query, array_size=1000, is_fetch=False):
        connection = cls.get_connection()
        result = None
        try:
            cursor = connection.cursor()
        except Exception:
            connection = cls.get_connection(new=True)  # Create new connection
            cursor = connection.cursor()
        cursor.execute(query)
        if is_fetch:
            if array_size > 0:
                result = cursor.fetchmany(array_size)
            else:
                result = cursor.fetchall()
        cursor.close()
        return result

    @classmethod
    def execute_many(cls, query, data):
        cursor = cls.get_cursor()
        for d in data:
            try:
                cursor.execute(query, d)
            except Exception as err:
                pass
                # info_debug.error(err)
        cursor.close()


connector = DBConnector(
    jar=f"{os.path.dirname(__file__)}/../raw/lib/h2.jar",
    driver="org.h2.Driver",
    server=f"jdbc:h2:{os.path.dirname(__file__)}/../raw/data/2209_phenotype/2209_phenotype/",
    user="sa",
    password="",
    database="2209_phenotype",
)


def ResultIterator(conn: DBConnection, query: str, arraysize: int = 1000) -> tuple:
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
    drop_query1 = f"DROP TABLE IF EXISTS EXOMISER.{table_name}_SCRAMBLE;"
    conn.execute_query(drop_query1)
    create_query = create_mapping(table_name)
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
def insert(conn: DBConnection, table_name: str, df: pd.DataFrame):
    """insert modified lines in scrambled table

    Args:
        conn (DBConnection): Database connection
        table_name (str): table name
        scramble_factor (float): scramble factor
        df (pd.DataFrame): Updated data to be inserted
    """
    """"""
    sql = insert_mapping(table_name)
    try:
        conn.execute_many(sql, df.values.tolist())
    except Exception as err:
        debug_log.debug(err)


@measure_time
def process_from_db(conn, table_name, scramble_factor, chunksize, count):
    i = 1
    end = 0
    for j in tqdm(range(end, count // chunksize)):
        start = (i - 1) * chunksize
        end = (chunksize * i) - 1
        try:
            data_to_update = select(conn, table_name, chunksize, start)
            mod = randomisation.ssp_randomisation(data_to_update, scramble_factor)
            insert(conn, table_name, mod)
        except Exception as err:
            info_debug.error(err)
        finally:
            i += 1


@measure_time
def process_from_file(conn, table_name, input, scramble_factor, chunksize):
    try:
        data = read_file(input, chunksize)
        # UNFAIR, I KNOW!
        count = int(subprocess.check_output(["wc", "-l", input]).split()[0])
        for data_to_update in tqdm(data, total=count // chunksize):
            mod = randomisation.ssp_randomisation(data_to_update, scramble_factor)
            insert(conn, table_name, mod)
    except Exception as err:
        debug_log.error(err)


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
    data_to_update = pd.DataFrame(ResultIterator(conn, select_query))
    return data_to_update


def rename_table(old_table_name: str, new_table_name: str):
    with connector as cnn:
        conn = DBConnection(cnn)
        create_query = f"""ALTER TABLE EXOMISER.{old_table_name} RENAME TO EXOMISER.{new_table_name};"""
        conn.execute_query(create_query)


def read_file(file_name, chunksize=10**6):
    for chunk in pd.read_csv(file_name, chunksize=chunksize, sep="\t"):
        yield chunk


def clean_aux_table(table_name: str):
    with connector as cnn:
        conn = DBConnection(cnn)
        drop_query1 = f"DROP TABLE IF EXISTS EXOMISER.{table_name}_ORIGINAL;"
        conn.execute_query(drop_query1)


def scramble_table(table_name: str, scramble_factor: float, input: str) -> None:
    with connector as cnn:
        try:
            conn = DBConnection(cnn)
            create_table(conn, table_name)
            info_log.info(
                f"Scrambling records from table: {table_name} using {scramble_factor} magnitude"
            )
            process_from_file(conn, table_name, input, scramble_factor, 100000)
            count_scramble = count_total(conn, f"{table_name}_SCRAMBLE")
            info_log.info(f"Scrambled records length: {count_scramble}")
            info_log.info("Done")
        except Exception as err:
            debug_log.error(err)
