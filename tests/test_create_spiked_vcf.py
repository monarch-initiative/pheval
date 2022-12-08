import os
import unittest
from pathlib import Path

from pheval.prepare.create_spiked_vcf import (
    ProbandVariantChecker,
    VcfHeader,
    VcfParser,
    VcfPicker,
    VcfSpiker,
)
from pheval.utils.phenopacket_utils import CausativeVariant


class TestVcfPicker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_file = Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_1.vcf")
        cls.vcf_dir = Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir")

    def test_pick_file_from_dir(self):
        self.assertTrue(
            "test_1.vcf" or "test_2.vcf" == VcfPicker(None, self.vcf_dir).pick_file_from_dir()
        )

    def test_pick_file(self):
        self.assertEqual("test_1.vcf", VcfPicker(self.vcf_file, None).pick_file().name)


class TestVcfParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.hg19_file = VcfParser(
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_1.vcf"
            ),
            False,
        )
        cls.hg38_file = VcfParser(
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_2.vcf"
            ),
            False,
        )

    def test_parse_assembly(self):
        self.assertEqual(self.hg19_file.parse_assembly(), ("GRCh37", True))
        self.assertEqual(self.hg38_file.parse_assembly(), ("GRCh38", False))

    def test_parse_sample_id(self):
        self.assertEqual(self.hg19_file.parse_sample_id(), "TEMPLATE")
        self.assertEqual(self.hg38_file.parse_sample_id(), "HG001")

    def test_parse_header(self):
        self.assertTrue(type(self.hg19_file.parse_header()) == VcfHeader)


class TestProbandVariantChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.causative_variants = [
            CausativeVariant(
                Path(
                    os.path.dirname(os.path.realpath(__file__))
                    + "/input_dir/test_phenopacket_1.json"
                ),
                "TEST1",
                "GRCh37",
                "1",
                886190,
                "G",
                "A",
                "heterozygous",
            )
        ]
        cls.vcf_header_compatible = VcfHeader("TEMPLATE", "GRCh37", True)
        cls.vcf_header_incompatible = VcfHeader("HG001", "GRCh38", False)

    def test_check_variant_assembly(self):
        self.assertTrue(
            ProbandVariantChecker(
                self.causative_variants, self.vcf_header_compatible
            ).check_variant_assembly()
            == []
        )
        self.assertTrue(
            ProbandVariantChecker(
                self.causative_variants, self.vcf_header_incompatible
            ).check_variant_assembly()
            != []
        )


class TestVcfSpiker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_spiker = VcfSpiker(
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_phenopacket_1.json"
            ),
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_1.vcf"
            ),
            Path("test_output_dir/"),
            [
                CausativeVariant(
                    Path(
                        os.path.dirname(os.path.realpath(__file__))
                        + "/input_dir/test_phenopacket_1.json"
                    ),
                    "TEST1",
                    "GRCh37",
                    "1",
                    886190,
                    "G",
                    "A",
                    "heterozygous",
                )
            ],
            VcfHeader("TEMPLATE", "GRCh37", True),
        )
        cls.variant = CausativeVariant(
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_phenopacket_1.json"
            ),
            "TEST1",
            "GRCh37",
            "1",
            886190,
            "G",
            "A",
            "heterozygous",
        )

    def test_construct_variant(self):
        variant_record = self.vcf_spiker.construct_variant(self.variant)
        self.assertEqual(type(variant_record), list)
        self.assertEqual(len(variant_record), 10)
        self.assertTrue(
            "chr1" in variant_record and "886190" in variant_record and "0/1" in variant_record[9]
        )

    def test_construct_vcf_records(self):
        vcf_contents = self.vcf_spiker.construct_vcf_records(False)
        self.assertTrue(len(vcf_contents) == 68)
        self.assertTrue(
            "chr1" in vcf_contents[59],
            "886190" in vcf_contents[59] and "SPIKED_VARIANT" in vcf_contents[59],
        )

    def test_construct_header(self):
        updated_vcf_contents = self.vcf_spiker.construct_header(False)
        for line in updated_vcf_contents:
            if line.startswith("#CHROM"):
                self.assertTrue("TEST1" in line and "TEMPLATE" not in line)


#
