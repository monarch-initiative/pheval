import os
import shutil
import unittest
from pathlib import Path

from pheval.infra.exomiserdb import DBConnection, ExomiserDB, _format_row
from pheval.utils.exomiser import semsim_to_exomiserdb

PHENO_FOLDER = os.path.abspath("./testdata/phenotype/2302_phenotype")
H2_JAR = os.path.abspath("./tests/lib/h2-1.4.199.jar")


class TestExomiserDBIngestion(unittest.TestCase):
    def setUp(self) -> None:
        _set_bkp_db()

    def tearDown(self):
        _clean_db()

    def select_data(self, query: str):
        edb = ExomiserDB(f"{PHENO_FOLDER}/2302_phenotype_test", H2_JAR)
        with edb.connector as cnn:
            conn = DBConnection(cnn)
            cursor = conn.get_cursor()
            cursor.execute(query)
            res = cursor.fetchall()
            return list(res)

    def test_semsim_to_exomiserdb(self):
        "tests semsim ingestion into exomiser phenotypic db"
        input_path = os.path.abspath("./testdata/semsim/hp-mp.semsim.tsv")
        query = "SELECT * FROM EXOMISER.HP_MP_MAPPINGS"
        prev_res = self.select_data(query)
        self.assertEqual(len(prev_res), 0)
        semsim_to_exomiserdb(
            input_path=input_path,
            object_prefix="MP",
            subject_prefix="HP",
            db_path=f"{PHENO_FOLDER}/2302_phenotype_test",
            h2_jar=H2_JAR,
        )
        new_res = self.select_data(query)
        self.assertEqual(len(new_res), 9)

    def test_format_row(self):
        mapping_id = 1
        data = {
            "subject_id": "subject1",
            "subject_label": "label1",
            "object_id": "object1",
            "object_label": "label2",
            "jaccard_similarity": 0.5,
            "ancestor_information_content": 0.6,
            "phenodigm_score": 0.7,
            "ancestor_id": "ancestor1,ancestor2",
            "ancestor_label": "ancestor_label1",
        }
        expected_output = "(1, 'subject1', 'label1', 'object1', 'label2', 0.5, 0.6, 0.7, 'ancestor1', 'ancestor_label1')"  # noqa
        self.assertEqual(_format_row(mapping_id, data), expected_output)


def _set_bkp_db():
    """copy an empty database to be populated in each test"""
    shutil.copyfile(
        f"{PHENO_FOLDER}/2302_phenotype_test.mv.db.bkp", f"{PHENO_FOLDER}/2302_phenotype_test.mv.db"
    )
    shutil.copyfile(
        f"{PHENO_FOLDER}/2302_phenotype_test.trace.db.bkp",
        f"{PHENO_FOLDER}/2302_phenotype_test.trace.db",
    )


def _clean_db():
    """remove populated database in each test"""
    Path(f"{PHENO_FOLDER}/2302_phenotype_test.mv.db").unlink(missing_ok=True)
    Path(f"{PHENO_FOLDER}/2302_phenotype_test.trace.db").unlink(missing_ok=True)
