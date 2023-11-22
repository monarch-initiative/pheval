import os
import tempfile
import unittest

import pandas as pd

from pheval.utils.utils import semsim_scramble


class TestSemsimUtils(unittest.TestCase):
    "Tests Semsim utilitary functions"

    def setUp(self) -> None:
        self._temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._temp_dir.cleanup()

    def test_semsim_scramble(self):
        "semsim_scramble function test"
        semsim_file = "./testdata/semsim/hp-mp.semsim.tsv"
        output_file = os.path.join(self._temp_dir.name, "test.sql")
        semsim_scramble(
            input=semsim_file,
            output=output_file,
            columns_to_be_scrambled=["jaccard_similarity"],
            scramble_factor=1,
        )
        exomiser_columns = [
            "subject_id",
            "subject_label",
            "subject_source",
            "object_id",
            "object_label",
            "object_source",
            "ancestor_id",
            "ancestor_label",
            "ancestor_source",
            "object_information_content",
            "subject_information_content",
            "ancestor_information_content",
            "jaccard_similarity",
            "dice_similarity",
            "phenodigm_score",
        ]

        with open(output_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            input_file = pd.read_csv(semsim_file, sep="\t")
            self.assertEqual(len(lines), len(input_file) + 1)
            self.assertEqual(list(input_file.columns.values), exomiser_columns)
