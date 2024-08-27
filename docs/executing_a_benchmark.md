# Executing a Benchmark

PhEval is designed for benchmarking algorithms across various datasets. To execute a benchmark using PhEval, you need to: 

1. Execute your runner; generating the PhEval standardised TSV outputs for gene/variant/disease prioritisation.
2. Configure the benchmarking parameters.
3. Run the benchmark.

PhEval will generate various performance reports, allowing you to easily compare the effectiveness of different algorithms.

## After the Runner Execution

After executing a run, you may be left with an output directory structure like so:

```tree
.
├── pheval_disease_results
│   ├── patient_1-pheval_disease_result.tsv
├── pheval_gene_results
│   ├── patient_1-pheval_gene_result.tsv
├── pheval_variant_results
│   ├── patient_1-pheval_variant_result.tsv
├── raw_results
│   ├── patient_1.json
├── results.yml
└── tool_input_commands
    └── tool_input_commands.txt
```
Whether you have populated `pheval_disease_results`, `pheval_gene_results`, and `pheval_variant_results` directories will depend on what is specified in the `config.yaml` for the runner execution. It is the results in these directories that are consumed in the benchmarking to produce the statistical comparison reports.

## Benchmarking Configuration File

To configure the benchmarking parameters, a YAML configuration file should be created and supplied to the CLI command.

An outline of the configuration file structure follows below:

```yaml
benchmark_name: exomiser_14_benchmark
runs:
  - run_identifier: run_identifier_1
    results_dir: /path/to/results_dir_1
    phenopacket_dir: /path/to/phenopacket_dir
    gene_analysis: True
    variant_analysis: False
    disease_analysis: True
    threshold:
    score_order: descending
  - run_identifier: run_identifier_2
    results_dir: /path/to/results_dir_2
    phenopacket_dir: /path/to/phenopacket_dir
    gene_analysis: True
    variant_analysis: True
    disease_analysis: True
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

The `benchmark_name` is what will be used to name the duckdb database that will contain all the ranking and binary statistics as well as comparisons between runs. The name provided should not have any whitespace or special characters.

### Runs section

The `runs` section specifies which run configurations should be included in the benchmarking. For each run configuration you will need to populate the following parameters:

- `run_identifier`: The identifier associated with the run - this should be meaningful as it will be used in the naming in tables and plots. 
- `results_dir`: The full path to the root directory where the directories `pheval_gene_results`/`pheval_variant_results`/`pheval_disease_results` can be found.
- `phenopacket_dir`: The full path to the phenopacket directory used during the runner execution.
- `gene_analysis`: Boolean specifying whether to perform benchmarking for gene prioritisation analysis.
- `variant_analysis`: Boolean specifying whether to perform benchmarking for variant prioritisation analysis
- `disease_analysis`: Boolean specifying whether to perform benchmarking for disease prioritisation analysis
- `threshold`: OPTIONAL score threshold to consider for inclusion of results. 
- `score_order`: Ordering of results for ranking. Either ascending or descending.

### Plot customisation section

The `plot_customisation` section specifies any additional customisation to the plots output from the benchmarking. Here you can specify title names for all the plots output, as well as the plot type for displaying the summary ranking stats. This section is split by the plots output from the gene, variant and disease prioritisation benchmarking. The parameters in this section do not need to be populated - however, if left blank it will default to generic titles. The parameters as follows are:

- `plot_type`: The plot type output for the summary rank stats plot. This can be either, bar_cumulative, bar_non_cumulative or bar_stacked.
- `rank_plot_title`: The customised title for the summary rank stats plot.
- `roc_curve_title`: The customised title for the ROC curve plot.
- `precision_recall_title` The customised title for the precision-recall curve plot.

## Executing the benchmark

After configuring the benchmarking YAML, executing the benchmark is relatively simple.

```bash
pheval-utils generate-benchmark-stats --run-yaml benchmarking_config.yaml
```


