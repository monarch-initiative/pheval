# PhEval Pipeline


## TLDR


The Pipeline presented on PhEval preprint (https://www.biorxiv.org/content/10.1101/2024.06.13.598672v1) was moved to a new repository - Monarch PhEval - https://github.com/monarch-initiative/monarch_pheval.

**NOTE: The default Monarch PhEval pipeline, as proposed in the paper preprint, requires approximately 1 TB of disk space. Learn how to modify the pipeline configuration in the later sections of this documentation to customize the experiment.**

### 1. Clone [Monarch PhEval](https://github.com/monarch-initiative/monarch_pheval)
  ```bash
  git clone https://github.com/monarch-initiative/monarch_pheval.git
  ```

### 2. Installing PhEval Pipeline dependencies
   Enter in the cloned folder and enter the following commands:

```bash
poetry shell
poetry install
```

### 3. Executing Pipeline

```bash
make pheval
```

## Pipeline Description

The Pipeline is divided in three main steps

### 1. Data Preparation Phase

The data preparation phase, checks the completeness of the disease, gene and variant input data and optionally preparing simulated VCF files if required, gives the user the ability to randomise phenotypic profiles using the PhEval corpus scramble command utility, allowing for the assessment of how well VGPAs handle noise and less specific phenotypic profiles when making predict.

### 2. Runner Phase

The runner phase is structured into three stages: prepare, run, and post-process.
 - The prepare step plays a crucial role in adapting the input data to meet the specific requirements of the tool. 
 - In the run step, the VGPA is executed, applying the selected algorithm to the prepared data and generating the tool-specific outputs. Within the run stage, an essential task is the generation of input command files for the algorithm. These files serve as collections of individual commands, each tailored to run the targeted VGPA on specific samples. These commands are configured with the appropriate inputs, outputs and specific configuration settings, allowing for the automated and efficient processing of large corpora. 
 - Finally, the post-processing step takes care of harmonising the tool-specific outputs into standardised PhEval TSV format, ensuring uniformity and ease of analysis of results from all VGPAs. In this context, the tool-specific output is condensed to provide only two essential elements, the entity of interest, which can either be a variant, gene, or disease, and its corresponding score. PhEval then assumes the responsibility of subsequent standardisation processes. This involves the reranking of the results in a uniform manner, ensuring that fair and comprehensive comparisons can be made between tools.

### 3. Analysis Phase

In the analysis phase, PhEval generates comprehensive statistical reports based on
standardised outputs from the runner phase.

## Customising PhEval Pipeline Experiments 

The phEval pipeline is orchestrated using a Makefile strategy. Therefore, to describe a new experiment in the pipeline, the user needs to generate a Makefile workflow based on a configuration file.

In the resources folder are the following files responsible for Makefile generation:

📦resources  
┣ 📜Makefile.j2  
┣ 📜custom.Makefile  
┣ 📜generatemakefile.sh  
┗ 📜pheval-config.yaml  

Let's start describing `pheval-config.yaml` structure.

### PhEval Configuration File

#### Directories Section

The `data` and `tmp` properties are mandatory and must be specified in this section.

- `data` property refers to the folder location where the necessary phenotypic data for the pipeline will be downloaded and extracted.
- `tmp` property points to the folder where all temporary intermediate files will be generated.

```yaml
directories:
  data: data
  tmp: data/tmp
```

#### Corpora Section


The `corpora` section specifies which corpus will be used in the experiment. In this example is defined [LIRICAL](https://pubmed.ncbi.nlm.nih.gov/32755546/) corpus, A small comparison corpus created for benchmarking the [LIRICAL](https://pubmed.ncbi.nlm.nih.gov/32755546/) system which contains 385 case reports.

The user needs to specify corpus id and it must be equals to the corpora folder structure, e.g.

📦corpora  
 ┃ ┣ 📂lirical  
 ┃ ┣ ┣ 📂small_version  
 ┃ ┣ ┣ ┣ 📂phenopackets  
 ┃ ┣ ┣ ┣ ┣ 📜PATIENT1.json  
 ┃ ┣ ┣ ┣ ┣ 📜PATIENT2.json  
 ┃ ┣ ┣ ┣ 📂vcf  
 ┃ ┣ ┣ ┣ ┣ 📜PATIENT1.vcf.gz  
 ┃ ┣ ┣ ┣ ┣ 📜PATIENT2.vcf.gz  
 ┃ ┣ ┣ ┣ 📜corpus.yml  
 ┃ ┣ ┣ ┣ 📜template_exome_hg19.vcf.gz  

```yaml
corpora:
  - id: lirical
    variant: small_version
```

#### Configs Section


The `configs` section holds all custom configurations for the different VGPAs.
It must declare:
- tool: VGPA tool name.
- id: it's an arbiratry unique identifier that will be used in the `runs` section
- version: VGPA tool version

```yaml
configs:
  - tool: phen2gene
    id: phen2gene-1.2.3
    version: 1.2.3
```

`configs` section can also deal with special VGPA data preparation steps, for example,  Semantic Similarity ingestions into Exomiser phenotypic database e.g.

```yaml
  configs:
  - tool: exomiser
    id: exomiser-semsim-ingest-13.3.0
    version: 13.3.0
    phenotype: 2309
    preprocessing:
      - phenio-monarch-hp-hp.0.4.semsimian.sql
```    
`phenotype` property describes the Exomiser phenotype database version and the `preprocessing` section will execute SQL scripts into that phenotypic database.


#### Runs Section

The "runs" section will integrate all previously described sections and pass them to pheval VGPA for concrete execution.

- `tool` property specifies which runner will be called
- `corpus` and `corpusvariant` must match properties declared on the [corpora section](#corpora-section).
- `version` should correspond to the tool version
- `configuration` must match the id described on the [configuration section](#configs-section).

```yaml
runs:
  - tool: exomiser
    corpus: lirical
    corpusvariant: small_version
    version: 13.3.0
    configuration: exomiser-semsim-ingest-13.3.0
```