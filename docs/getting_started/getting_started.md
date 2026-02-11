# Getting started with PhEval

This section helps new users orient themselves: what PhEval is, what problems it is designed to solve, and where to go next.

## Video walkthrough

If you prefer a guided walkthrough, start here:

- [**Introduction to PhEval and running a simple benchmark**](https://www.youtube.com/watch?v=nIPzVN99UWc)

## What is PhEval?

PhEval (Phenotypic Inference Evaluation Framework) supports the **systematic evaluation of phenotype-driven prioritisation methods**.

Phenotype-driven tools are sensitive to changes in:

- Ontologies and their versions
- Phenotype, gene, and disease mappings
- Semantic similarity methods
- Underlying knowledge resources and cohorts

PhEval addresses this by providing a framework for **reproducible and controlled evaluation**, with standardised outputs that can be compared across tools, versions, and experimental conditions.

## High-level workflow

At a high level, PhEval workflows involve:

1. Installing PhEval and any required plugins
2. Producing PhEval-standardised results using tool-specific runners
3. Benchmarking and analysing those results

Each of these steps is documented in more detail elsewhere.

---

## Suggested reading order

- New to PhEval: start with [Installation](installation.md), then read [Plugins and runners](../using_pheval/plugins_and_runners.md/)
- Benchmarking existing outputs: go to [Benchmarking](../benchmarking/index.md)
- Running robustness or simulation experiments: go to [Utilities](../utilities/index.md)
- Extending PhEval with a new tool: go to [Developing a PhEval plugin](../resources_for_contributors/developing_a_pheval_plugin.md)