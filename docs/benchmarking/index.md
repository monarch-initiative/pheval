# Benchmarking and analysis

This section describes how PhEval is used to **benchmark and compare phenotype-driven prioritisation methods** once tool execution has completed.

Benchmarking in PhEval is designed to support **controlled, reproducible evaluation** across:

- Tools and tool versions
- Cohorts and simulation strategies
- Ontology and knowledge-base updates

This section focuses on *analysis*, not execution. Tools are executed via runners provided by plugins.

---

## What benchmarking means in PhEval

In PhEval, benchmarking refers to the process of:

- Consuming **PhEval standardised results** produced by runners
- Computing rank-based and binary classification metrics
- Comparing performance across multiple runs
- Generating plots and summary statistics for interpretation

Benchmarking operates over one or more completed runs and assumes that tool execution has already taken place.

---

## High-level benchmarking workflow

A typical benchmarking workflow consists of:

1. **Execute one or more runners**  
   Runners produce PhEval-standardised outputs for gene, variant, and/or disease prioritisation.

2. **Configure benchmarking parameters**  
   A YAML configuration file specifies which runs to include and how benchmarking should be performed.

3. **Run benchmarking and analysis**  
   PhEval utilities compute metrics, comparisons, and plots across the specified runs.

Each of these steps is described in more detail in the following pages.

---

## What benchmarking produces

Benchmarking generates:

- Ranking-based statistics
- Binary classification statistics
- Comparative summaries between runs
- Plots for visual comparison
- A singular DuckDB database containing all computed metrics and comparisons

These outputs support both exploratory analysis and formal evaluation.
