# PhEval Pipeline Exomiser Runner


## Step by Step to PhEval Run Pipeline (with ExomiserRunner)

### 1. Download [Exomiser Software](https://github.com/exomiser/Exomiser/releases)
```bash
wget https://github.com/exomiser/Exomiser/releases/download/13.2.0/exomiser-cli-13.2.0-distribution.zip
```
### 2. Download Phenotype Data
```bash
wget https://data.monarchinitiative.org/exomiser/latest/2302_hg19.zip
wget https://data.monarchinitiative.org/exomiser/latest/2302_hg38.zip
wget https://data.monarchinitiative.org/exomiser/latest/2302_phenotype.zip
```

### 3. Unzip data

  ```bash
  # unzip the distribution and data files - this will create a directory called 'exomiser-cli-13.1.0' in the current working directory
  unzip exomiser-cli-13.2.0-distribution.zip
  unzip 2302_*.zip -d exomiser-cli-13.2.0/data
  ```


### 4. Clone PhEval repo and follow steps described in Pipeline Documentation:

```bash
git clone https://github.com/monarch-initiative/pheval.git
cd pheval
poetry shell
poetry install
```

### 5. Set PhEval Config YAML File

```yaml

directories:
  tmp: data/tmp
  exomiser: /path_where_exomiser_was_extracted
  phenotype: /path_where_phenotype_was_extracted
  workspace: /tmp/pheval

corpora:
  - id: small_test
    scrambled:
      - factor: 0.5
      - factor: 0.7
    custom_variants:
      - id: no_phenotype

configs:
  - tool: exomiser
    version: 13.2.0
    configuration: default
    exomiser_db: semsim1


runs:
  - tool: exomiser
    configuration: default
    corpus: small_test
    corpusvariant: scrambled-0.5
    version: 13.2.0


```

### 6. Generate Makefile based on configuration

```bash
bash ./resources/generatemakefile.sh
```


### 7. Exomiser Runner requires the following configuration

The `config.yaml` file should be formatted like the example below and must be placed in `exomiser: /path_where_exomiser_was_extracted` declared in `pheval-config.yaml` file.

```yaml
tool: exomiser
tool_version: 13.2.0
phenotype_only: False # NOTE phenotype-only preset analysis should only be run with Exomiser versions >= 13.2.0
tool_specific_configuration_options:
  environment: local
  exomiser_software_directory: .
  analysis_configuration_file: preset-exome-analysis.yml
  max_jobs: 0
  application_properties:
    remm_version:
    cadd_version:
    hg19_data_version: 2302
    hg19_local_frequency_path:
    hg38_data_version: 2302
    phenotype_data_version: 2302
    cache_type:
    cache_caffeine_spec:
  post_process:
    score_name: combinedScore
    sort_order: DESCENDING
```

### 7. PhEval Run

```bash
make pheval run
```