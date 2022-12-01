import os
import unittest
from pathlib import Path

from pheval.prepare import create_spiked_vcf


class TestVcfPicker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_file = Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_1.vcf")
        cls.vcf_dir = Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir")

    def test_vcf_picker(self):
        self.assertEqual(
            "test_1.vcf",
            create_spiked_vcf.VcfPicker(self.vcf_file, None).pick_file().name,
        )
        selected_file = create_spiked_vcf.VcfPicker(None, self.vcf_dir).pick_file()
        self.assertTrue("test_1.vcf" == selected_file.name or "test_2.vcf" == selected_file.name)


class TestVcfParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.hg19_file = create_spiked_vcf.VcfParser(
            Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_1.vcf")
        )
        cls.hg38_file = create_spiked_vcf.VcfParser(
            Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_2.vcf")
        )

    def test_parse_header(self):
        self.assertEqual(
            create_spiked_vcf.VcfHeader("TEMPLATE", "GRCh37", True),
            self.hg19_file.parse_header(),
        )
        self.assertEqual(
            create_spiked_vcf.VcfHeader("HG001", "GRCh38", True), self.hg38_file.parse_header()
        )


class TestProbandVariantChecker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.causative_variants = [
            create_spiked_vcf.CausativeVariant(
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
        cls.vcf_header_compatible = create_spiked_vcf.VcfHeader("TEMPLATE", "GRCh37", True)
        cls.vcf_header_incompatible = create_spiked_vcf.VcfHeader("HG001", "GRCh38", True)

    def test_check_variant_assembly(self):
        self.assertTrue(
            create_spiked_vcf.ProbandVariantChecker(
                self.causative_variants, self.vcf_header_compatible
            ).check_variant_assembly()
            == []
        )
        self.assertTrue(
            create_spiked_vcf.ProbandVariantChecker(
                self.causative_variants, self.vcf_header_incompatible
            ).check_variant_assembly()
            != []
        )


class TestVcfSpiker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_spiker = create_spiked_vcf.VcfSpiker(
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_phenopacket_1.json"
            ),
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_1.vcf"
            ),
            Path("test_output_dir/"),
            [
                create_spiked_vcf.CausativeVariant(
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
            create_spiked_vcf.VcfHeader("TEMPLATE", "GRCh37", True),
        )
        cls.variant = create_spiked_vcf.CausativeVariant(
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
        vcf_contents = self.vcf_spiker.construct_vcf_records()
        self.assertTrue(len(vcf_contents) == 68)
        self.assertTrue(
            "chr1" in vcf_contents[59],
            "886190" in vcf_contents[59] and "SPIKED_VARIANT" in vcf_contents[59],
        )

    def test_construct_header(self):
        updated_vcf_contents = self.vcf_spiker.construct_header()
        for line in updated_vcf_contents:
            if line.startswith("#CHROM"):
                self.assertTrue("TEST1" in line and "TEMPLATE" not in line)


class TestVcfWriter(unittest.TestCase):
    vcf_contents = create_spiked_vcf.VcfSpiker(
        Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_phenopacket_1.json"),
        Path(os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_vcf_dir/test_1.vcf"),
        Path("test_output_dir/"),
        [
            create_spiked_vcf.CausativeVariant(
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
        create_spiked_vcf.VcfHeader("TEMPLATE", "GRCh37", True),
    ).construct_header()

    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_writer = create_spiked_vcf.VcfWriter(
            Path(
                os.path.dirname(os.path.realpath(__file__)) + "/input_dir/test_phenopacket_1.json"
            ),
            cls.vcf_contents,
            Path(Path(os.path.dirname(os.path.realpath(__file__)) + "/test_copied_vcf.json")),
        )

    @classmethod
    def tearDownClass(cls) -> None:
        Path.unlink(
            Path(Path(os.path.dirname(os.path.realpath(__file__)) + "/test_copied_vcf.json"))
        )

    def test_copy_template_to_new_file(self):
        self.vcf_writer.copy_template_to_new_file()
        self.assertTrue(
            os.path.exists(
                Path(Path(os.path.dirname(os.path.realpath(__file__)) + "/test_copied_vcf.json"))
            )
        )

    def test_write_vcf(self):
        self.vcf_writer.write_vcf()
        with open(
            Path(Path(os.path.dirname(os.path.realpath(__file__)) + "/test_copied_vcf.json"))
        ) as f:
            self.assertEqual(len(f.readlines()), 68)
        f.close()
