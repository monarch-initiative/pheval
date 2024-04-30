import unittest
from pathlib import Path

from pheval.prepare.create_spiked_vcf import (
    VcfHeader,
    VcfHeaderParser,
    VcfSpiker,
    check_variant_assembly,
    read_vcf,
)
from pheval.utils.phenopacket_utils import (
    GenomicVariant,
    IncompatibleGenomeAssemblyError,
    ProbandCausativeVariant,
)

hg19_vcf = [
    "##fileformat=VCFv4.1\n",
    '##FILTER=<ID=PASS,Description="All filters passed">\n',
    '##FILTER=<ID=LowQual,Description="Low quality">\n',
    "##FORMAT=<ID=AD,Number=.,Type=Integer,Description="
    '"Allelic depths for the ref and alt alleles in the order listed">\n',
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description="
    '"Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description="
    '"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n',
    "##contig=<ID=chrM,length=16571>\n",
    "##contig=<ID=chr1,length=249250621>\n",
    "##contig=<ID=chr2,length=243199373>\n",
    "##contig=<ID=chr3,length=198022430>\n",
    "##contig=<ID=chr4,length=191154276>\n",
    "##contig=<ID=chr5,length=180915260>\n",
    "##contig=<ID=chr6,length=171115067>\n",
    "##contig=<ID=chr7,length=159138663>\n",
    "##contig=<ID=chr8,length=146364022>\n",
    "##contig=<ID=chr9,length=141213431>\n",
    "##contig=<ID=chr10,length=135534747>\n",
    "##contig=<ID=chr11,length=135006516>\n",
    "##contig=<ID=chr12,length=133851895>\n",
    "##contig=<ID=chr13,length=115169878>\n",
    "##contig=<ID=chr14,length=107349540>\n",
    "##contig=<ID=chr15,length=102531392>\n",
    "##contig=<ID=chr16,length=90354753>\n",
    "##contig=<ID=chr17,length=81195210>\n",
    "##contig=<ID=chr18,length=78077248>\n",
    "##contig=<ID=chr19,length=59128983>\n",
    "##contig=<ID=chr20,length=63025520>\n",
    "##contig=<ID=chr21,length=48129895>\n",
    "##contig=<ID=chr22,length=51304566>\n",
    "##contig=<ID=chrX,length=155270560>\n",
    "##contig=<ID=chrY,length=59373566>\n",
    "##reference=file:///share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.5/hg19/ucsc.hg19.fasta\n",
    "##bcftools_viewVersion=1.9+htslib-1.9\n",
    "##bcftools_viewCommand=view -s NIST7035 project.NIST.hc.snps.indels.vcf; Date=Mon Mar 25 21:52:10 2019\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEMPLATE\n",
    "chr1\t886179\trs34136865\tCA\tC\t190.31\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=1;FS=0;MLEAC=4;MLEAF=1;MQ=29;MQ0=0;QD=29.42\tGT:AD:DP:GQ:PL\t1/1:.:.:12:180,12,0\n",
    "chr1\t886182\trs35678314\tTG\tT\t190.31\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=1;FS=0;MLEAC=4;MLEAF=1;MQ=29;MQ0=0;QD=34.28\tGT:AD:DP:GQ:PL\t1/1:.:.:12:180,12,0\n",
    "chr1\t886186\trs145470050\tC\tCAG\t190.31\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=5;FS=0;MLEAC=4;MLEAF=1;MQ=35.54;MQ0=0;QD=19.03\tGT:AD:DP:GQ:PL\t1/1:0,4:4:12:180,12,0\n",
    "chr1\t886788\trs10465242\tG\tA\t49.17\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=2;FS=0;MLEAC=2;MLEAF=1;MQ=60;MQ0=0;QD=24.58\tGT:AD:DP:GQ:PL\t1/1:0,2:2:6:75,6,0\n",
    "chr3\t61575424\trs4688637\tC\tA\t72.31\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=4;FS=0;MLEAC=4;MLEAF=1;MQ=49.84;MQ0=0;QD=18.08\tGT:AD:DP:GQ:PL\t1/1:0,2:2:6:49,6,0\n",
    "chr3\t61578327\trs10866049\tT\tC\t55.17\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=2;FS=0;MLEAC=2;MLEAF=1;MQ=60;MQ0=0;QD=27.58\tGT:AD:DP:GQ:PL\t1/1:0,2:2:6:81,6,0\n",
    "chr3\t61580859\trs4447758\tG\tT\t233.26\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=6;FS=0;MLEAC=4;MLEAF=1;MQ=60;MQ0=0;QD=31.48\tGT:AD:DP:GQ:PL\t1/1:0,1:1:3:45,3,0\n",
    "chr3\t61580861\trs4621322\tT\tC\t233.26\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=6;FS=0;MLEAC=4;MLEAF=1;MQ=60;MQ0=0;QD=31.41\tGT:AD:DP:GQ:PL\t1/1:0,1:1:3:45,3,0\n",
    "chr3\t61580912\trs4447759\tC\tT\t134.35\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=5;FS=0;MLEAC=4;MLEAF=1;MQ=60;MQ0=0;QD=26.87\tGT:AD:DP:GQ:PL\t1/1:0,1:1:3:38,3,0\n",
    "chr3\t61582116\trs145762950\tT\tTTTG\t144.42\t.\t"
    "AC=0;AF=1;AN=0;DB;DP=4;FS=0;MLEAC=2;MLEAF=1;MQ=60;MQ0=0;QD=12.04\tGT:AD:DP:GQ:PL\t./.:.:.:.:.\n",
    "chr3\t61593146\trs35485211\tT\tTGC\t55.13\t.\t"
    "AC=2;AF=1;AN=2;DB;DP=2;FS=0;MLEAC=2;MLEAF=1;MQ=60;MQ0=0;QD=13.78\tGT:AD:DP:GQ:PL\t1/1:0,2:2:6:90,6,0",
]

hg38_vcf = [
    "##fileformat=VCFv4.2\n",
    "##fileDate=20160824\n",
    '##CL=vcffilter -i - -o - --javascript "function record() {HG001.PS=\\".\\";}"\n',
    "##contig=<ID=1,length=248956422,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=2,length=242193529,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=3,length=198295559,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=4,length=190214555,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=5,length=181538259,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=6,length=170805979,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=7,length=159345973,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=8,length=145138636,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=9,length=138394717,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=10,length=133797422,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=11,length=135086622,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=12,length=133275309,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=13,length=114364328,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=14,length=107043718,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=15,length=101991189,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=16,length=90338345,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=17,length=83257441,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=18,length=80373285,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=19,length=58617616,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=20,length=64444167,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=21,length=46709983,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=22,length=50818468,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=X,length=156040895,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=Y,length=57227415,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "##contig=<ID=M,length=16569,assembly=human_GRCh38_no_alt_analysis_set.fasta>\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG001\n",
    "10\t19819808\t.\tC\tT\t50\tPASS\tplatforms=5;platformnames=Illumina,PacBio,CG,10X,Solid;datasets=5;"
    "datasetnames=HiSeqPE300x,CCS15kb_20kb,CGnormal,10XomiumLR,SolidSE75bp;callsets=7;"
    "callsetnames=HiSeqPE300xGATK,CCS15kb_20kbDV,CCS15kb_20kbGATK4,CGnormal,HiSeqPE300xfreebayes,10XLRGATK,"
    "SolidSE75GATKHC;datasetsmissingcall=IonExome;callable="
    "CS_HiSeqPE300xGATK_callable,CS_CCS15kb_20kbDV_callable,CS_10XLRGATK_callable,CS_CCS15kb_20kbGATK4_callable,"
    "CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable\tGT:PS:DP:ADALL:AD:GQ\t0/1:.:964:173,188:239,249:855\n",
    "10\t19830688\t.\tCT\tC\t50\tPASS\tplatforms=5;platformnames=PacBio,CG,Illumina,10X,Solid;datasets=5;"
    "datasetnames=CCS15kb_20kb,CGnormal,HiSeqPE300x,10XomiumLR,SolidSE75bp;callsets=7;"
    "callsetnames=CCS15kb_20kbDV,CCS15kb_20kbGATK4,CGnormal,HiSeqPE300xfreebayes,HiSeqPE300xGATK,10XLRGATK,"
    "SolidSE75GATKHC;datasetsmissingcall=IonExome;callable="
    "CS_CCS15kb_20kbDV_callable,CS_10XLRGATK_callable,CS_CCS15kb_20kbGATK4_callable,CS_CGnormal_callable,"
    "CS_HiSeqPE300xfreebayes_callable\tGT:PS:DP:ADALL:AD:GQ\t1/1:.:769:0,292:70,137:410\n",
    "10\t19830740\t.\tG\tGT\t50\tPASS\tplatforms=3;platformnames=Illumina,PacBio,10X;datasets=3;"
    "datasetnames=HiSeqPE300x,CCS15kb_20kb,10XomiumLR;callsets=4;"
    "callsetnames=HiSeqPE300xGATK,CCS15kb_20kbDV,10XLRGATK,HiSeqPE300xfreebayes;datasetsmissingcall="
    "CCS15kb_20kb,CGnormal,IonExome,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable;filt=CS_CGnormal_filt;"
    "difficultregion=GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5,GRCh38_SimpleRepeat_imperfecthomopolgt10_slop5"
    "\tGT:PS:DP:ADALL:AD:GQ\t1/1:.:575:30,271:7,237:128",
]


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
        cls.vcf_spiker_new_variant_chrom = VcfSpiker(
            hg19_vcf,
            [
                ProbandCausativeVariant(
                    "TEST1",
                    "GRCh37",
                    GenomicVariant("X", 123450, "G", "A"),
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
        cls.structural_variant_vcf_spiker = VcfSpiker(
            hg38_vcf,
            [
                ProbandCausativeVariant(
                    "TEST1",
                    "GRCh38",
                    GenomicVariant("5", 134858794, "N", "DEL"),
                    "heterozygous",
                    "SVTYPE=DEL;END=135099433;SVLEN=-269958227",
                )
            ],
            VcfHeader("TEMPLATE", "GRCh38", True),
        )
        cls.variant = ProbandCausativeVariant(
            "TEST1",
            "GRCh37",
            GenomicVariant("1", 886190, "G", "A"),
            "heterozygous",
        )
        cls.structural_variant = ProbandCausativeVariant(
            proband_id="test-subject-1",
            assembly="GRCh38",
            variant=GenomicVariant(chrom="5", pos=134858794, ref="N", alt="DEL"),
            genotype="heterozygous",
            info="SVTYPE=DEL;END=135099433;SVLEN=-269958227",
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
                ".",
                "GT",
                "0/1\n",
            ],
        )

    def test_construct_variant_structural_variant(self):
        self.assertEqual(
            self.structural_variant_vcf_spiker.construct_variant_entry(self.structural_variant),
            [
                "chr5",
                "134858794",
                ".",
                "N",
                "<DEL>",
                "100",
                "PASS",
                "SVTYPE=DEL;END=135099433;SVLEN=-269958227",
                "GT",
                "0/1\n",
            ],
        )

    def test_construct_vcf_records_single_variant(self):
        self.assertEqual(
            self.vcf_spiker.construct_vcf_records("template.vcf")[40],
            "chr1\t886190\t.\tG\tA\t100\tPASS\t.\t" "GT\t0/1\n",
        )

    def test_construct_vcf_records_multiple_variants(self):
        updated_records = self.vcf_spiker_multiple_variants.construct_vcf_records("template.vcf")
        self.assertEqual(
            updated_records[40],
            "chr1\t886190\t.\tG\tA\t100\tPASS\t.\t" "GT\t0/1\n",
        )
        self.assertEqual(
            updated_records[45],
            "chr3\t61580860\t.\tG\tA\t100\tPASS\t.\t" "GT\t1/1\n",
        )

    def test_construct_vcf_records_new_variant_pos(self):
        updated_records = self.vcf_spiker_new_variant_chrom.construct_vcf_records("template.vcf")
        self.assertEqual(
            updated_records[48],
            "chrX\t123450\t.\tG\tA\t100\tPASS\t.\tGT\t0/1\n",
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
