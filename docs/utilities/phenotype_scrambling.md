# Phenotype scrambling utilities

This page documents utilities used to **introduce noise or perturbations into phenotype data**.
These commands are typically used to assess the robustness and sensitivity of phenotype-driven prioritisation methods.

They operate on existing phenotype data and do not execute tools or perform benchmarking directly.

---

## Purpose

Phenotype scrambling utilities are used to:

- Simulate noisy or incomplete phenotypic observations
- Evaluate how sensitive methods are to phenotype quality
- Test robustness under controlled perturbations

These experiments are useful when comparing tools, parameterisations, or ontology versions.

---

## Scrambling phenopackets

The `scramble-phenopackets` command generates perturbed versions of existing phenopackets.

The scrambled phenopackets can then be used as inputs to runners for execution and benchmarking.

### Example: scramble a phenopacket corpus

Generate scrambled phenopackets from an existing corpus:

```bash
pheval-utils scramble-phenopackets \
  --phenopacket-dir phenopackets/ \
  --output-dir scrambled_phenopackets/ \
  --scramble-factor 0.7 \
  --local-ontology-cache ./hp.obo
```

> !!! note "Notes: "
    - The original phenopackets are not modified.
    - Scrambled outputs are written to a separate directory.
    - The resulting phenopackets can be used directly with plugin-provided runners.

---

## How phenotype scrambling fits into a workflow

A typical robustness experiment using phenotype scrambling might look like:

1. Prepare a clean phenopacket corpus
2. Generate scrambled phenopackets
3. Run tools via plugin-provided runners using `pheval run`
4. Benchmark and compare performance against the original results

Scrambling utilities are optional and primarily used in experimental or methodological evaluations.

