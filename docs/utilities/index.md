# Utilities

This section documents the **utility commands** provided with PhEval.
These utilities support data preparation, manipulation, and experimental workflows that sit *around* tool execution and benchmarking.

They are not required for every use case, but are commonly used when preparing cohorts, running robustness experiments, or standardising inputs.

## What the utilities are for

PhEval utilities are designed to help with tasks such as:

- Preparing phenopacket corpora for evaluation
- Updating identifiers and mappings in existing data
- Generating synthetic or perturbed inputs for robustness testing
- Supporting benchmarking and downstream analysis

They are provided via the `pheval-utils` command-line interface, which is installed automatically when installing PhEval.

## Scope and boundaries

Utilities are intentionally separated from:

- **tool execution**, which is handled by runners via plugins
- **benchmarking logic**, which is documented in the Benchmarking section

This separation keeps workflows modular and reproducible.

## Categories of utilities

The utilities fall broadly into the following categories.

### [Data preparation](data_preparation.md)

Commands used to prepare or normalise input data before execution, including:

- Preparing corpora of phenopackets
- Updating gene symbols and identifiers
- Ensuring consistent formats for downstream tools


### [Phenotype scrambling and noise experiments](phenotype_scrambling.md)

Commands used to introduce noise or perturbations into phenotype data.
These are commonly used to assess robustness and sensitivity of phenotype-driven methods.


### [Variant-related utilities](variant_utilities.md)

Commands that operate on variant-level data, such as creating spiked VCFs for controlled evaluation experiments.


### [Resource and mapping updates](resource_updates.md)

Commands used to download or update shared resources, such as ontology mappings and identifier tables.


## How utilities fit into a workflow

A typical workflow using utilities might look like:

1. Prepare or update phenopacket data using utilities
2. Execute tools via runners using `pheval run`
3. Benchmark and analyse results

Not all workflows require utilities; they are optional building blocks.

