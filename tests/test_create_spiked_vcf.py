import unittest
from pathlib import Path

from pheval.prepare.create_spiked_vcf import (
    VcfHeader,
    VcfHeaderParser,
    VcfPicker,
    VcfSpiker,
    check_variant_assembly,
    read_vcf,
)
from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    IncompatibleGenomeAssemblyError,
    ProbandCausativeVariant,
)

hg19_vcf = read_vcf(Path("./tests/input_dir/test_vcf_dir/test_1.vcf"))
hg38_vcf = read_vcf(Path("./tests/input_dir/test_vcf_dir/test_2.vcf.gz"))


class TestVcfPicker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_file = Path("./tests/input_dir/test_vcf_dir/test_1.vcf")
        cls.vcf_dir = Path("./tests/input_dir/test_vcf_dir/test_2.vcf.gz")

    def test_pick_file_from_dir(self):
        self.assertTrue(
            "test_1.vcf" or "test_2.vcf.gz" == VcfPicker(None, self.vcf_dir).pick_file_from_dir()
        )

    def test_pick_file(self):
        self.assertEqual("test_1.vcf", VcfPicker(self.vcf_file, None).pick_file().name)


class TestReadVcf(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_file = Path("./tests/input_dir/test_vcf_dir/test_1.vcf")
        cls.gzipped_vcf_file = Path("./tests/input_dir/test_vcf_dir/test_2.vcf.gz")

    def test_read_vcf(self):
        self.assertTrue(type(read_vcf(self.vcf_file)), list[str])

    def test_read_vcf_gzipped(self):
        for line in read_vcf(self.gzipped_vcf_file):
            self.assertFalse(isinstance(line, (bytes, bytearray)))


class TestVcfHeaderParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.hg19_file = VcfHeaderParser(hg19_vcf)
        cls.hg38_file = VcfHeaderParser(hg38_vcf)

    def test_parse_assembly(self):
        self.assertEqual(self.hg19_file.parse_assembly(), ("GRCh37", True))
        self.assertEqual(self.hg38_file.parse_assembly(), ("GRCh38", False))

    def test_parse_sample_id(self):
        self.assertEqual(self.hg19_file.parse_sample_id(), "TEMPLATE")
        self.assertEqual(self.hg38_file.parse_sample_id(), "HG001")

    def test_parse_vcf_header(self):
        self.assertEqual(self.hg19_file.parse_vcf_header(), VcfHeader("TEMPLATE", "GRCh37", True))
        self.assertEqual(self.hg38_file.parse_vcf_header(), VcfHeader("HG001", "GRCh38", False))


class TestCheckVariantAssembly(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.causative_variants = [
            ProbandCausativeVariant(
                "TEST1",
                "GRCh37",
                GenomicVariant("1", 886190, "G", "A"),
                "heterozygous",
            )
        ]
        cls.incorrect_causative_variant_assembly = [
            ProbandCausativeVariant(
                "TEST1",
                "hg10",
                GenomicVariant("1", 886190, "G", "A"),
                "heterozygous",
            )
        ]
        cls.multiple_assemblies = [
            ProbandCausativeVariant(
                "TEST1",
                "hg19",
                GenomicVariant("1", 886190, "G", "A"),
                "heterozygous",
            ),
            ProbandCausativeVariant(
                "TEST1",
                "hg38",
                GenomicVariant("1", 886190, "G", "A"),
                "heterozygous",
            ),
        ]
        cls.vcf_header_compatible = VcfHeader("TEMPLATE", "GRCh37", True)
        cls.vcf_header_incompatible = VcfHeader("HG001", "GRCh38", False)

    def test_check_variant_assembly(self):
        self.assertEqual(
            check_variant_assembly(
                self.causative_variants, self.vcf_header_compatible, Path("/path/to/phenopacket")
            ),
            None,
        )

    def test_check_variant_assembly_phenopacket_assemblies(self):
        with self.assertRaises(ValueError):
            check_variant_assembly(
                self.multiple_assemblies, self.vcf_header_compatible, Path("/path/to/phenopacket")
            )

    def test_check_variant_assembly_incompatible_vcf(self):
        with self.assertRaises(IncompatibleGenomeAssemblyError):
            check_variant_assembly(
                self.causative_variants, self.vcf_header_incompatible, Path("/path/to/phenopacket")
            )

    def test_check_variant_assembly_incorrect_variant_assembly(self):
        with self.assertRaises(IncompatibleGenomeAssemblyError):
            check_variant_assembly(
                self.incorrect_causative_variant_assembly,
                self.vcf_header_incompatible,
                Path("/path/to/phenopacket"),
            )


class TestVcfSpiker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.vcf_spiker = VcfSpiker(
            hg19_vcf,
            [
                ProbandCausativeVariant(
                    "TEST1",
                    "GRCh37",
                    GenomicVariant("1", 886190, "G", "A"),
                    "heterozygous",
                )
            ],
            VcfHeader("TEMPLATE", "GRCh37", True),
        )
        cls.vcf_spiker_multiple_variants = VcfSpiker(
            hg19_vcf,
            [
                ProbandCausativeVariant(
                    "TEST1",
                    "GRCh37",
                    GenomicVariant("1", 886190, "G", "A"),
                    "heterozygous",
                ),
                ProbandCausativeVariant(
                    "TEST1",
                    "GRCh37",
                    GenomicVariant("3", 61580860, "G", "A"),
                    "homozygous",
                ),
            ],
            VcfHeader("TEMPLATE", "GRCh37", True),
        )
        cls.variant = ProbandCausativeVariant(
            "TEST1",
            "GRCh37",
            GenomicVariant("1", 886190, "G", "A"),
            "heterozygous",
        )

    def test_construct_variant(self):
        self.assertEqual(
            self.vcf_spiker.construct_variant_entry(self.variant),
            [
                "chr1",
                "886190",
                ".",
                "G",
                "A",
                "100",
                "PASS",
                "SPIKED_VARIANT_HETEROZYGOUS",
                "GT:AD:DP:GQ:PL",
                "0/1:0,2:2:12:180,12,0\n",
            ],
        )

    def test_construct_vcf_records_single_variant(self):
        self.assertEqual(
            self.vcf_spiker.construct_vcf_records()[59],
            "chr1\t886190\t.\tG\tA\t100\tPASS\tSPIKED_VARIANT_HETEROZYGOUS\t"
            "GT:AD:DP:GQ:PL\t0/1:0,2:2:12:180,12,0\n",
        )

    def test_construct_vcf_records_multiple_variants(self):
        updated_records = self.vcf_spiker_multiple_variants.construct_vcf_records()
        self.assertEqual(
            updated_records[59],
            "chr1\t886190\t.\tG\tA\t100\tPASS\tSPIKED_VARIANT_HETEROZYGOUS\t"
            "GT:AD:DP:GQ:PL\t0/1:0,2:2:12:180,12,0\n",
        )
        self.assertEqual(
            updated_records[64],
            "chr3\t61580860\t.\tG\tA\t100\tPASS\tSPIKED_VARIANT_HOMOZYGOUS\t"
            "GT:AD:DP:GQ:PL\t1/1:0,2:2:12:180,12,0\n",
        )

    def test_construct_header(self):
        self.assertEqual(
            [
                line
                for line in self.vcf_spiker.construct_header(hg19_vcf)
                if line.startswith("#CHROM")
            ],
            ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST1\n"],
        )
