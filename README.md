# PhEval - Phenotypic Inference Evaluation Framework

[![PyPI](https://img.shields.io/pypi/v/pheval)](https://pypi.org/project/pheval/)
![Build Status](https://img.shields.io/github/actions/workflow/status/monarch-initiative/pheval/pypi-publish.yml?branch=main)
![License](https://img.shields.io/github/license/monarch-initiative/pheval)
![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)
![Issues](https://img.shields.io/github/issues/monarch-initiative/pheval)

PhEval (Phenotypic Inference Evaluation Framework) is a **modular, reproducible benchmarking framework** for evaluating **phenotype-driven prioritisation tools**, such as gene, variant, and disease prioritisation algorithms.

It is designed to support **fair comparison across tools, tool versions, datasets, and knowledge updates**, addressing a long-standing gap in standardised evaluation for phenotype-based methods.

📖 **Full documentation:** https://monarch-initiative.github.io/pheval/
---

## Why PhEval?

Evaluating phenotype-driven prioritisation tools is challenging because performance depends on many moving parts, including:

- Phenotype representations and noise
- Ontology structure and versioning
- Gene and disease mappings
- Tool-specific scoring and ranking strategies
- Input cohorts and simulation approaches

PhEval provides a framework that makes these factors **explicit, controlled, and comparable**.

Key features:

- **Standardised outputs** across tools
- **Reproducible benchmarking** with recorded metadata
- **Plugin-based architecture** for extensibility
- **Separation of execution and evaluation**
- Support for **gene, variant, and disease prioritisation**

---

## Installation

PhEval requires **Python 3.10 or later**.

Install from PyPI:

```bash
pip install pheval
```

This installs:

* The core pheval CLI (for running tools via plugins)
* `pheval-utils` (for data preparation, benchmarking, and analysis)

Verify installation:

```bash
pheval --help
pheval-utils --help
```

## How PhEval is used

PhEval workflows typically consist of three phases:

1.	Prepare data 
    Prepare and manipulate phenopackets and related inputs (e.g. VCFs).
2. Run tools 
   Execute phenotype-driven prioritisation tools via plugin-provided runners using:
   ```bash
   pheval run --runner <runner_name> ...
   ```
3. Benchmark and analyse
   Compare results across runs using standardised metrics and plots.

Each phase is documented in detail in the user documentation.

## Plugins and runners

PhEval itself is tool-agnostic.

Support for specific tools is provided via plugins, which implement runners responsible for:

* Preparing tool inputs 
* Executing the tool 
* Converting raw outputs into PhEval standardised results

A list of available plugins is maintained in the documentation:

Plugins: https://monarch-initiative.github.io/pheval/plugins/

Each plugin repository contains tool-specific installation instructions and examples.

## Documentation

The PhEval documentation is organised by audience and task:
* Getting started: installation and first steps 
* Using PhEval: running tools, plugins, and workflows 
* Utilities: data preparation, phenopacket manipulation, simulations 
* Benchmarking: executing benchmarks, metrics, and plots 
* Developer documentation: plugin development and API reference

Start here: https://monarch-initiative.github.io/pheval/

## Contributions

Contributions are welcome across:

* Code 
* Documentation 
* Testing 
* Plugins and integrations

## Citation

If you use **PhEval** in your research, please cite the following publication:

> **Bridges, Y., Souza, V. d., Cortes, K. G., et al.**  
> *Towards a standard benchmark for phenotype-driven variant and gene prioritisation algorithms: PhEval – Phenotypic Inference Evaluation Framework.*  
> **BMC Bioinformatics** 26, 87 (2025).  
> https://doi.org/10.1186/s12859-025-06105-4


