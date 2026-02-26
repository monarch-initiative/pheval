# Data preparation utilities

This page documents **data preparation utilities** provided by PhEval.
These commands are used to prepare, normalise, and organise input data *before* running phenotype-driven tools via plugins.

This page only covers commands related to **data preparation**.
Variant spiking and other specialised workflows are documented elsewhere.

---

## Purpose

Data preparation utilities help to:

- Construct phenopacket corpora for evaluation
- Normalise gene identifiers
- Ensure consistent input structure across cohorts
- Reduce technical variability unrelated to tool performance

These steps are particularly important when benchmarking across tools, versions, or knowledge resources.

---

## Preparing a phenopacket corpus

The `prepare-corpus` command is used to prepare a directory of phenopackets for downstream analysis.

Typical use cases include:

- Validating that phenopackets contain the required records
- Preparing separate corpora for gene-, disease-, or variant-based analyses
- Optionally generating associated VCFs for variant-based workflows

### Basic example

Prepare a corpus of phenopackets for gene-based analysis:

```bash
pheval-utils prepare-corpus \
  --phenopacket-dir phenopackets/ \
  --gene-analysis \
  --output-dir prepared_corpus/
```

Prepare a corpus of phenopackets for gene-based analysis and update all gene identifiers to Ensembl IDs:

```bash
pheval-utils prepare-corpus \
  --phenopacket-dir phenopackets/ \
  --gene-analysis \
  --gene-identifier ensembl_id \
  --output-dir prepared_corpus/
```

### Variant-based analysis example

Prepare a corpus for variant-based analysis using an hg38 VCF template:

```bash
pheval-utils prepare-corpus \
  --phenopacket-dir phenopackets/ \
  --variant-analysis \
  --output-dir prepared_corpus/
```

Prepare a corpus for variant-based analysis and spike variants into an hg38 VCF template:

```bash
pheval-utils prepare-corpus \
  --phenopacket-dir phenopackets/ \
  --variant-analysis \
  --hg38-template-vcf hg38_template.vcf \
  --output-dir prepared_corpus/
```


> !!! note "Notes: "
    - At least one of `--variant-analysis`, `--gene-analysis`, or `--disease-analysis` should be specified.
    - For variant-based analysis, a VCF template or directory is required.
    - The prepared output directory is used as input to runners provided by plugins.

---

## Updating phenopackets and identifiers

The `update-phenopackets` command is used to update gene symbols and identifiers in existing phenopackets.

This is useful when:

- Phenopackets contain outdated gene identifiers
- A consistent identifier scheme is required across a cohort
- Benchmarking is performed across different database or ontology versions

### Example: update a directory of phenopackets

Update phenopackets to include Ensembl gene identifiers:

```bash
pheval-utils update-phenopackets \
  --phenopacket-dir phenopackets/ \
  --gene-identifier ensembl_id \
  --output-dir updated_phenopackets/
```

### Example: update a single phenopacket

```bash
pheval-utils update-phenopackets \
  --phenopacket-path case_001.json \
  --gene-identifier hgnc_id \
  --output-dir updated_case/
```

---

## How data preparation fits into a workflow

A typical workflow using data preparation utilities looks like:

1. Collect or generate phenopackets
2. Prepare and normalise phenopackets using data preparation utilities
3. Run tools via plugin-provided runners using `pheval run`
4. Benchmark and analyse the resulting outputs

Not all workflows require all preparation steps, but these utilities help ensure reproducibility and consistency.

