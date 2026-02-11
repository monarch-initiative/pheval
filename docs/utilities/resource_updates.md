# Resource updates

This page documents utilities used to **download and update shared resources** required by PhEval workflows.

---

## Purpose

Resource update utilities are used to:

- Ensure identifier mappings are up to date
- Maintain consistency across benchmarking runs
- Reduce errors caused by outdated or missing reference data

Keeping resources updated is recommended, particularly when running new experiments or comparing results across time.

---

## Updating shared resources

The `update` command downloads and refreshes shared resources used by PhEval and its plugins.

This includes:

- The MONDO SSSOM mapping file from the Monarch Initia
- The HGNC complete gene set from the HGNC download site

### Example: update all resources

```bash
pheval update
```

The command will download the latest versions of the supported resources and store them in PhEval’s configured data directory.

---

## When to run resource updates

You should consider running `pheval update` when:

- Installing PhEval for the first time
- Starting a new benchmarking experiment
- Updating PhEval or plugin versions

Running updates explicitly helps ensure clarity about which resources are being used.

---

## Resource provenance and reproducibility

When tools are executed via `pheval run`, information about shared resources is **automatically recorded** in the run metadata.

Each run produces a `results.yml` file that captures, among other details:

- The tool and tool version
- The execution timestamp
- The corpus used
- The download dates of shared resources

Example:

```yaml
tool: <TOOL_NAME>
tool_version: <TOOL_VERSION>
config: <TOOL_CONFIG>
run_timestamp: <TIMESTAMP>
corpus: <CORPUS_NAME>
mondo_download_date: <MONDO_DOWNLOAD_DATE>
hgnc_download_date: <HGNC_DOWNLOAD_DATE>
tool_specific_configuration_options: null
```

By recording resource download timestamps alongside each run, PhEval enables:

- Tracing which ontology and mapping versions were in use
- Comparison of results across runs performed at different times
- Transparent reporting and reproducibility

As long as `results.yml` files are preserved alongside benchmarking outputs, manual tracking of resource versions is not required.

Updating resources using `pheval update` should therefore be treated as an **explicit experimental choice** and interpreted in the context of the recorded run metadata.

---

## How resource updates fit into a workflow

A typical workflow involving resource updates might look like:

1. Install PhEval
2. Update shared resources using `pheval update`
3. Prepare input data and corpora
4. Run tools via plugin-provided runners
5. Benchmark and analyse results

