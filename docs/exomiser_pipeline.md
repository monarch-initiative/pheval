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
  unzip 2302_hg19.zip -d exomiser-cli-13.2.0/data
  unzip 2302_hg38.zip -d exomiser-cli-13.2.0/data
  ```


### 4. Clone PhEval repo and follow steps described in Pipeline Documentation:

```bash
git clone https://github.com/monarch-initiative/pheval.git
cd pheval
poetry shell
poetry install
pip install pheval.exomiser
```

### 5. Set PhEval Config YAML File

```yaml

directories:
  tmp: data/tmp
  exomiser: /path_where_exomiser_was_extracted
  phenotype: /path_where_phenotype_was_extracted
  workspace: /pheval's_path # path where pheval was cloned

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
variant_analysis: True
gene_analysis: True
disease_analysis: True
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


### 8. Preset Exome Analysis File

Exomiser requires a *preset-exome-analysis.yml* file saved at **/path_where_exomiser_was_extracted/preset-exome-analysis.yml**

This is an example of preset-exome-analysis.yml file

```yaml
## Exomiser Analysis Template.
# These are all the possible options for running exomiser. Use this as a template for
# your own set-up.
---
analysisMode: PASS_ONLY
inheritanceModes: {
  AUTOSOMAL_DOMINANT: 0.1,
  AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
  AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
  X_DOMINANT: 0.1,
  X_RECESSIVE_HOM_ALT: 0.1,
  X_RECESSIVE_COMP_HET: 2.0,
  MITOCHONDRIAL: 0.2
}
frequencySources: [
    THOUSAND_GENOMES,
    TOPMED,
    UK10K,

    ESP_AFRICAN_AMERICAN, ESP_EUROPEAN_AMERICAN, ESP_ALL,

    EXAC_AFRICAN_INC_AFRICAN_AMERICAN, EXAC_AMERICAN,
    EXAC_SOUTH_ASIAN, EXAC_EAST_ASIAN,
    EXAC_FINNISH, EXAC_NON_FINNISH_EUROPEAN,
    EXAC_OTHER,

    GNOMAD_E_AFR,
    GNOMAD_E_AMR,
  #        GNOMAD_E_ASJ,
    GNOMAD_E_EAS,
    GNOMAD_E_FIN,
    GNOMAD_E_NFE,
    GNOMAD_E_OTH,
    GNOMAD_E_SAS,

    GNOMAD_G_AFR,
    GNOMAD_G_AMR,
  #        GNOMAD_G_ASJ,
    GNOMAD_G_EAS,
    GNOMAD_G_FIN,
    GNOMAD_G_NFE,
    GNOMAD_G_OTH,
    GNOMAD_G_SAS
]
# Possible pathogenicitySources: (POLYPHEN, MUTATION_TASTER, SIFT), (REVEL, MVP), CADD, REMM
# REMM is trained on non-coding regulatory regions
# *WARNING* if you enable CADD or REMM ensure that you have downloaded and installed the CADD/REMM tabix files
# and updated their location in the application.properties. Exomiser will not run without this.
pathogenicitySources: [ REVEL, MVP ]
#this is the standard exomiser order.
steps: [
    failedVariantFilter: { },
    variantEffectFilter: {
      remove: [
          FIVE_PRIME_UTR_EXON_VARIANT,
          FIVE_PRIME_UTR_INTRON_VARIANT,
          THREE_PRIME_UTR_EXON_VARIANT,
          THREE_PRIME_UTR_INTRON_VARIANT,
          NON_CODING_TRANSCRIPT_EXON_VARIANT,
          NON_CODING_TRANSCRIPT_INTRON_VARIANT,
          CODING_TRANSCRIPT_INTRON_VARIANT,
          UPSTREAM_GENE_VARIANT,
          DOWNSTREAM_GENE_VARIANT,
          INTERGENIC_VARIANT,
          REGULATORY_REGION_VARIANT
      ]
    },
    frequencyFilter: { maxFrequency: 2.0 },
    pathogenicityFilter: { keepNonPathogenic: true },
    inheritanceFilter: { },
    omimPrioritiser: { },
    hiPhivePrioritiser: { }
]

```

### 9. PhEval Run

```bash
make pheval run
```
