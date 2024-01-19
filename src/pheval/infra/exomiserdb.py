# -*- coding: cp936 -*-
import logging as log
import os
from pathlib import Path

import jaydebeapi
import polars as pl
from tqdm import tqdm

info_log = log.getLogger("info")
info_debug = log.getLogger("debug")


class DBConnector:
    def __init__(
        self, jar: Path, driver: str, server: str, database: str, user: str, password: str
    ):
        self.jar = jar
        self.driver = driver
        self.server = server
        self.database = database
        self.user = user
        self.password = password
        self.dbconn = None

    def create_connection(self) -> jaydebeapi.Connection:
        """creates h2 database connection"""
        return jaydebeapi.connect(
            self.driver,
            f"{self.server}{self.database}",
            [self.user, self.password],
            self.jar,
        )

    def __enter__(self) -> jaydebeapi.Connection:
        self.dbconn = self.create_connection()
        return self.dbconn

    def __exit__(self, *other):
        self.dbconn.close()


class DBConnection:
    connection = None

    def __init__(self, connection):
        DBConnection.connection = connection

    @classmethod
    def get_connection(cls) -> jaydebeapi.Connection:
        """Creates return new Singleton database connection"""
        return DBConnection.connection

    def close(self):
        return self.connection.close()

    @classmethod
    def get_cursor(cls) -> jaydebeapi.Cursor:
        connection = cls.get_connection()
        return connection.cursor()


class ExomiserDB:
    def __init__(self, db_path: Path):
        try:
            self.connector = DBConnector(  # noqa
                jar=os.path.join(os.path.dirname(__file__), "../../../lib/h2-1.4.199.jar"),
                driver="org.h2.Driver",
                server=f"jdbc:h2:{db_path}",
                user="sa",
                password="",
                database="",
            )
        except Exception as e:
            print("An exception occurred", e)

    def import_from_semsim_file(self, input_file: Path, subject_prefix: str, object_prefix: str):
        """imports semsim tsv profile into exomiser phenotype database

        Args:
            input_file (Path): semsim profile
            subject_prefix (str): Subject Prefix. e.g HP
            object_prefix (str): Object Prefix. e.g MP
        """
        with self.connector as cnn:
            conn = DBConnection(cnn)
            reader = pl.read_csv_batched(input_file, separator="\t")
            batch_length = 5
            batches = reader.next_batches(batch_length)
            cursor = conn.get_cursor()
            # # TODO: Refactor this
            with open(input_file, "r") as f:
                total = sum(1 for line in f)
            pbar = tqdm(total=total - 1)
            mapping_id = 1
            while batches:
                input_data = pl.concat(batches)
                sql = _semsim2h2(input_data, object_prefix, subject_prefix, mapping_id=mapping_id)
                cursor.execute(sql)
                len_input_data = len(input_data)
                mapping_id += len_input_data
                pbar.update(len_input_data)

                batches = reader.next_batches(batch_length)


def _format_row(mapping_id, data):
    """format row in a exomiser database way

    Args:
        mapping_id (_type_): row sequencial id
        data (_type_): row data
    """
    # TODO:Improve string escaping. Replace this code with parametrised query
    return f"""({mapping_id}, '{data['subject_id']}', '{data['subject_label'].replace("'", "")}', '{data['object_id']}', '{data['object_label'].replace("'", "")}', {data['jaccard_similarity']}, {data['ancestor_information_content']}, {data['phenodigm_score']}, '{data['ancestor_id'].split(",")[0]}', '{data['ancestor_label'].replace("'", "")}')"""  # noqa


def _semsim2h2(
    input_data: pl.DataFrame, subject_prefix: str, object_prefix: str, mapping_id=1
) -> None:
    """This function is responsible for generate sql insertion query for each semsim profile row

    Args:
        input_data (pl.DataFrame): input data. (e.g. semantic similarity profile file)
        subject_prefix (str): subject prefix. (e.g HP)
        object_prefix (str): object prefix. (e.g MP)
        mapping_id (int, optional): MAPPING_ID.
    """
    sql = ""
    if mapping_id == 1:
        sql += f"TRUNCATE TABLE EXOMISER.{subject_prefix}_{object_prefix}_MAPPINGS;\n"

    object_id = (
        f"{object_prefix}_ID_HIT" if subject_prefix == object_prefix else f"{object_prefix}_ID"
    )
    object_term = (
        f"{object_prefix}_HIT_TERM" if subject_prefix == object_prefix else f"{object_prefix}_TERM"
    )
    sql += f"""INSERT INTO EXOMISER.{subject_prefix}_{object_prefix}_MAPPINGS
(MAPPING_ID, {subject_prefix}_ID, {subject_prefix}_TERM, {object_id}, {object_term}, SIMJ, IC, SCORE, LCS_ID, LCS_TERM)
VALUES"""
    rows = [
        _format_row(data=frame, mapping_id=mapping_id + jdx)
        for jdx, frame in enumerate(input_data.iter_rows(named=True))
    ]
    sql += ",\n".join(rows) + ";"
    return sql
