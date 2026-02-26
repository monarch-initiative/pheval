# Executing a Benchmark

This page describes how to execute a benchmark, configure benchmarking parameters, and interpret the resulting outputs.

It assumes that one or more PhEval runs have already been completed using plugin-provided runners.

---

## After runner execution

After executing a run, an output directory structure similar to the following is produced:

```tree
.
├── pheval_disease_results
│   ├── patient_1-disease_result.parquet
├── pheval_gene_results
│   ├── patient_1-gene_result.parquet
├── pheval_variant_results
│   ├── patient_1-variant_result.parquet
├── raw_results
│   ├── patient_1.json
├── results.yml
└── tool_input_commands
    └── tool_input_commands.txt
```

Which result directories are present depends on the configuration used during runner execution.

The contents of the `pheval_*_results` directories are consumed during benchmarking.

---

## Benchmarking configuration file

Benchmarking is configured using a YAML file supplied to the CLI.

### Example configuration

```yaml
benchmark_name: tool_version_update_benchmark
runs:
  - run_identifier: run_identifier_1
    results_dir: /path/to/results_dir_1
    phenopacket_dir: /path/to/phenopacket_dir
    gene_analysis: true
    variant_analysis: false
    disease_analysis: true
    threshold:
    score_order: descending
  - run_identifier: run_identifier_2
    results_dir: /path/to/results_dir_2
    phenopacket_dir: /path/to/phenopacket_dir
    gene_analysis: true
    variant_analysis: true
    disease_analysis: true
    threshold:
    score_order: descending
plot_customisation:
  gene_plots:
    plot_type: bar_cumulative
    rank_plot_title:
    roc_curve_title:
    precision_recall_title:
  disease_plots:
    plot_type: bar_cumulative
    rank_plot_title:
    roc_curve_title:
    precision_recall_title:
  variant_plots:
    plot_type: bar_cumulative
    rank_plot_title:
    roc_curve_title:
    precision_recall_title:
```

The `benchmark_name` is used to name the DuckDB database that stores benchmarking statistics.
It should not contain whitespace or special characters.

---

## Runs section

Each entry in the `runs` list specifies a completed run to include in the benchmark.

Required fields:

- `run_identifier` → A human-readable identifier used in tables and plots.
- `results_dir` → Path to the directory containing `pheval_gene_results`, `pheval_variant_results`, and/or `pheval_disease_results`.
- `phenopacket_dir`  →Path to the phenopacket directory used during runner execution.
- `gene_analysis`, `variant_analysis`, `disease_analysis` →Boolean flags indicating which analyses to include.

Optional fields:

- `threshold` → Score threshold for result inclusion.
- `score_order` → Ranking order (`ascending` or `descending`).

---

## Plot customisation

The `plot_customisation` section allows optional control over plot appearance.

Available options:

- `plot_type` → One of `bar_cumulative`, `bar_non_cumulative`, or `bar_stacked`.
- `rank_plot_title` → Custom title for ranking summary plots.
- `roc_curve_title` → Custom title for ROC plots.
- `precision_recall_title` → Custom title for precision–recall plots.

If left unspecified, default titles and plot types are used.

---

## Executing the benchmark

Once the configuration file is prepared, benchmarking can be executed with:

```bash
pheval-utils benchmark --run-yaml benchmarking_config.yaml
```

> !!! note "**Command Note:**"
    As of `pheval` version **0.5.0** onwards, the command is `benchmark`.  
    In earlier versions, the equivalent command was `generate-benchmark-stats`.
    See the [v0.5.1 release notes](https://github.com/monarch-initiative/pheval/releases/tag/0.5.1) for more details.


---

## Outputs and interpretation

Benchmarking produces:

- A DuckDB database containing computed statistics, and comparisons between runs
- Rank-based and binary classification plots

These outputs can be used to compare tools, configurations, and experimental conditions in a reproducible manner.
