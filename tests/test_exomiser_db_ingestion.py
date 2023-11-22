"__summary__"
import os
import shutil
import unittest
from pathlib import Path

from pheval.utils.exomiser import semsim_to_exomiserdb

PHENO_FOLDER = os.path.abspath("./testdata/phenotype/2302_phenotype")


class TestExomiserDBIngestion(unittest.TestCase):
    def setUp(self) -> None:
        _set_bkp_db()

    def tearDown(self):
        _clean_db()

    def test_semsim_to_exomiserdb(self):
        "tests semsim ingestion into exomiser phenotypic db"
        input_path = os.path.abspath("./testdata/semsim/hp-mp.semsim.tsv")
        semsim_to_exomiserdb(
            input_path=input_path,
            object_prefix="MP",
            subject_prefix="HP",
            db_path=f"{PHENO_FOLDER}/2302_phenotype_test",
        )


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
