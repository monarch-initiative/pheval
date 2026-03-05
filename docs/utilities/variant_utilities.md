# Variant utilities

This page documents utilities used to work with **variant-level data** in PhEval workflows.
These commands are primarily used to construct or manipulate VCF inputs for **variant-based evaluation experiments**.

They do not execute tools directly and do not perform benchmarking.

---

## Purpose

Variant utilities are used to:

- Generate VCFs containing known ("spiked") variants
- Support controlled variant-based evaluation experiments

These utilities are typically used in conjunction with phenopacket data and plugin-provided runners.

---

## Creating spiked VCFs

The `create-spiked-vcfs` command is used to generate VCF files containing known causal variants derived from phenopackets.

This is particularly useful when:

- Evaluating variant-based prioritisation methods
- Simulating realistic diagnostic scenarios
- Benchmarking tools that require both phenotypes and variants

The command supports both single phenopackets and directories of phenopackets.

---

### Example: create spiked VCFs from a phenopacket directory (hg38)

```bash
pheval-utils create-spiked-vcfs \
  --phenopacket-dir phenopackets/ \
  --hg38-template-vcf hg38_template.vcf \
  --output-dir spiked_vcfs/
```

---

### Example: create a spiked VCF from a single phenopacket (hg19)

```bash
pheval-utils create-spiked-vcfs \
  --phenopacket-path case_001.json \
  --hg19-template-vcf hg19_template.vcf \
  --output-dir spiked_vcf/
```

---

### Example: use a directory of VCF templates

Instead of a single template file, a directory of VCF templates can be provided:

```bash
pheval-utils create-spiked-vcfs \
  --phenopacket-dir phenopackets/ \
  --hg38-vcf-dir hg38_vcf_templates/ \
  --output-dir spiked_vcfs/
```

---

> !!! note "Notes and constraints:"
    - Exactly one of `--phenopacket-path` or `--phenopacket-dir` must be provided.
    - For each genome build, either a template VCF file or a directory of template VCFs must be supplied.
    - The generated VCFs are written to the specified output directory.
    - Spiked VCFs are typically consumed by runners that support variant-based analysis.

---

## How variant utilities fit into a workflow

A typical variant-based evaluation workflow might look like:

1. Prepare and normalise phenopackets
2. Generate spiked VCFs using variant utilities
3. Run tools via plugin-provided runners using `pheval run`
4. Benchmark and analyse variant-level results