# PheVal Benchmarking

Assesses the percentage of the top, top3 and top5 genes from Exomiser runs - retrieving the data from phenopackets containing the known variant and gene variants. 

## Set up

```

git clone --branch PhEval-Benchmarking https://github.com/monarch-initiative/pheval.git

cd pheval/pheval-benchmark

poetry shell

poetry install

```

## Basic usage

Details on commands:

```

python src/pheval_benchmark/cli.py

```
or

```

poetry run pheval

```
or

```

pheval

```

To create spiked VCF files for disease phenopackets:

```

pheval spike-vcf --phenopacket-dir src/pheval_benchmark/resources/disease_phenopackets --template-vcf template_exome_vcf_hg19.vcf.gz --output-dir spiked_vcf

```
To create noisy phenopackets from original disease phenopackets:

Example:

```

pheval create-noisy-phenopackets --phenopacket-dir src/pheval_benchmark/resources/disease_phenopackets --max-real-id 3 --number-of-parent-terms 2 --number-of-random-terms 4 --output-file-suffix added_noise

```

To create batch files to run with Exomiser:

Example:

```

pheval prepare-exomiser-batch --analysis PATH/TO/ANALYSIS/YML --phenopacket-dir src/pheval_benchmark/resources/disease_phenopackets --vcf-dir PATH/TO/SPIKED/VCF/DIR --batch-prefix Exomiser-run

```

To analyse the prioritisation of two Exomiser runs:

Example:

```

pheval assess-prioritisation --directory-list /PATH/TO/LIST/OF/DIRECTORIES --phenopacket-dir src/pheval_benchmark/resources/disease_phenopackets --output-prefix Comparison

```

First line of directory list should be the directory used as a baseline.

