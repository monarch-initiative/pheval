import os
import tempfile
import unittest

import pandas as pd

from pheval.utils.utils import semsim_convert, semsim_scramble


class TestSemsimUtils(unittest.TestCase):
    def setUp(self) -> None:
        self._temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._temp_dir.cleanup()

    def test_semsim_convert_exomiser(self):
        semsim_file = "./testdata/semsim/hp-mp.semsim.tsv"
        output_file = os.path.join(self._temp_dir.name, "test.tsv")
        semsim_convert(
            input=semsim_file,
            output=output_file,
            object_prefix="HP",
            subject_prefix="MP",
            format="exomiserdb",
        )
        with open(output_file, "r") as f:
            lines = f.readlines()
            input = pd.read_csv(semsim_file, sep="\t")
            assert len(lines) == len(input) + 1

    def test_semsim_convert_invalid_format(self):
        semsim_file = "./testdata/semsim/hp-mp.semsim.tsv"
        output_file = os.path.join(self._temp_dir.name, "test.tsv")

        with self.assertRaises(ValueError):
            semsim_convert(
                input=semsim_file,
                output=output_file,
                object_prefix="HP",
                subject_prefix="MP",
                format="invalid_format",
            )

    def test_semsim_scramble(self):
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

        with open(output_file, "r") as f:
            lines = f.readlines()
            input = pd.read_csv(semsim_file, sep="\t")
            assert len(lines) == len(input) + 1
            assert list(input.columns.values) == exomiser_columns
