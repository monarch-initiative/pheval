import os
import tempfile
import unittest

from pheval.utils.utils import semsim_randomisation, semsimconvert


class TestSemsimUtils(unittest.TestCase):
    def setUp(self) -> None:
        self._temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._temp_dir.cleanup()

    def test_semsimconvert(self):
        semsim_file = "./testdata/semsim/hp-mp.semsim.tsv"
        output_file = os.path.join(self._temp_dir.name, "test.tsv")
        semsimconvert(
            input=semsim_file, output=output_file, object_prefix="HP", subject_prefix="MP"
        )
        with open(output_file, "r") as f:
            lines = f.readlines()
            self.assertTrue("INSERT" in " ".join(lines))

    def test_semsim_randomisation(self):
        semsim_file = "./testdata/semsim/hp-mp.semsim.tsv"
        output_file = os.path.join(self._temp_dir.name, "test.sql")
        semsim_randomisation(input=semsim_file, output=output_file, scramble_factor=1)
        with open(output_file, "r") as f:
            lines = f.readlines()
            self.assertTrue("subject_id" in " ".join(lines))
