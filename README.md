# PhEval - Phenotypic Inference Evaluation Framework

![PyPI](https://img.shields.io/pypi/v/pheval)
![Build Status](https://img.shields.io/github/actions/workflow/status/monarch-initiative/pheval/pypi-publish.yml?branch=main)
![License](https://img.shields.io/github/license/monarch-initiative/pheval)
![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)
![Issues](https://img.shields.io/github/issues/monarch-initiative/pheval)

## Overview

The absence of standardised benchmarks and data standardisation for Variant and Gene Prioritisation Algorithms (VGPAs) presents a significant challenge in the field of genomic research. To address this, we developed PhEval, a novel framework designed to streamline the evaluation of VGPAs that incorporate phenotypic data. PhEval offers several key benefits:

- Automated Processes: Reduces manual effort by automating various evaluation tasks, thus enhancing efficiency.
- Standardisation: Ensures consistency and comparability in evaluation methodologies, leading to more reliable and standardised assessments.
- Reproducibility: Facilitates reproducibility in research by providing a standardised platform, allowing for consistent validation of algorithms.
- Comprehensive Benchmarking: Enables thorough benchmarking of algorithms, providing well-founded comparisons and deeper insights into their performance. 

PhEval is a valuable tool for researchers looking to improve the accuracy and reliability of VGPA evaluations through a structured and standardised approach.

For more information please see the full [documentation](https://monarch-initiative.github.io/pheval/).

## Download and Installation

1. Ensure you have Python 3.8 or greater installed.
2. Install with `pip`:
```bash
pip install pheval
```
3. See list of all PhEval utility commands:
```bash
pheval-utils --help
```

## Usage

The PhEval CLI offers a variety of commands categorised into two main types: **Runner Implementations** and **Utility Commands**. Below is an overview of each category, detailing how they can be utilised to perform various tasks within PhEval.

### Runner Implementations

The primary command used within PhEval is `pheval run`. This command is responsible for executing concrete VGPA runner implementations, that we sometimes term as plugins. By using pheval run, users can leverage these runner implementations to: execute the VGPA on a set of test corpora, produce tool-specific result outputs, and post-process tool-specific outputs to PhEval standardised TSV outputs.

Some concrete PhEval runner implementations include the [Exomiser runner](https://github.com/monarch-initiative/pheval.exomiser) and the [Phen2Gene runner](https://github.com/monarch-initiative/pheval.phen2gene). The full list of currently implemented runners can be found [here](https://monarch-initiative.github.io/pheval/plugins/)

Please read the [documentation](https://monarch-initiative.github.io/pheval/developing_a_pheval_plugin/) for a step-by-step for creating your own PhEval plugin. 

### Utility Commands

In addition to the main `run` command, PhEval provides a set of utility commands designed to enhance the overall functionality of the CLI. These commands can be used to set up and configure experiments, streamline data preparation, and benchmark the performance of various VGPA runner implementations. By utilising these utilities, users can optimise their experimental workflows, ensure reproducibility, and compare the efficiency and accuracy of different approaches. The utility commands offer a range of options that facilitate the customisation and fine-tuning to suit diverse research objectives.

#### Example Usage

To add noise to an existing corpus of phenopackets, this could be used to assess the robustness of VGPAs when less relevant or unreliable phenotype data is introduced:
```bash
pheval-utils scramble-phenopackets --phenopacket-dir /phenopackets --scramble-factor 0.5 --output-dir /scrambled_phenopackets_0.5
```

To update the gene symbols and identifiers to a specific namespace:
```bash
pheval-utils update-phenopackets --phenopacket-dir /phenopackets --output-dir /updated_phenopackets --gene-identifier ensembl_id
```

To prepare VCF files for a corpus of phenopackets, spiking in the known causative variants:
```bash
pheval-utils create-spiked-vcfs --phenopacket-dir /phenopackets --hg19-template-vcf /template_hg19.vcf --hg38-template-vcf /template_hg38.vcf --output-dir /vcf
```

Alternatively, you can wrap all corpus preparatory commands into a single step. Specifying `--variant-analysis`/`--gene-analysis`/`--disease-analysis` will check the phenopackets for complete records documenting the known entities. If template vcf(s) are provided this will spike VCFs with the known variant for the corpus. If a `--gene-identifier` is specified then the corpus of phenopackets is updated.
```bash
pheval-utils prepare-corpus \
    --phenopacket-dir /phenopackets \
    --variant-analysis \
    --gene-analysis \
    --gene-identifier ensembl_id \
    --hg19-template-vcf /template_hg19.vcf \
    --hg38-template-vcf /template_hg38.vcf \
    --output-dir /vcf
```

See the [documentation](https://monarch-initiative.github.io/pheval/executing_a_benchmark/) for instructions on benchmarking and evaluating the performance of various VGPAs.

