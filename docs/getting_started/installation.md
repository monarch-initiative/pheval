# Installation

This page explains how to install **PhEval** and verify that it is available on your system.  
Tool-specific execution details are handled by **plugins** and documented in their respective repositories.

---

## Prerequisites

Before installing PhEval, ensure you have:

- **Python 3.10 or newer**
- A working Python environment (virtualenv, conda, uv, or similar)
- `pip` or a compatible Python package manager

It is strongly recommended to install PhEval in an isolated environment.

---

## Create a virtual environment (recommended)

Using `venv`:

```bash
python -m venv pheval-env
source pheval-env/bin/activate
```

On macOS/Linux, you should now see your environment name in the shell prompt.

## Install PhEval

Install the latest released version from PyPI:

```bash
pip install pheval
```

This installs:

* The PhEval command-line interface 
* The `pheval-utils` command-line utilities
* Shared utilities used by plugins and runners

## Verify the installation

Check that PhEval is installed correctly:

```bash
pheval --help
```

You should see output similar to:

```bash
Usage: pheval [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  run     Execute a phenotype-driven tool via a runner
  update  Download or update required mapping resources
```

Verify that the utility commands are also available:

```bash
pheval-utils --help
```

You should see a list of available utility commands, including data preparation and benchmarking.

If either command is not found, ensure your virtual environment is activated.

## Install a plugin

PhEval does not execute tools directly.
Instead, plugins provide runners that implement tool-specific preparation, execution, and post-processing.

After installing PhEval, install one or more plugins corresponding to the tools you want to evaluate.

Each plugin:

* Documents its own installation requirements 
* Exposes one or more runners 
* Explains how to invoke those runners using `pheval run`

See the Plugins section for a list of available plugins and links to their documentation.

## Update mapping resources (recommended)

Some workflows require up-to-date ontology and identifier mappings.

You can download or refresh these resources using:

```bash
pheval update
```

This step is recommended before running benchmarks, particularly when working with gene or disease identifiers.

## Next steps

Once PhEval is installed:

* Learn how plugins and runners fit into the execution model
→ [Plugins and runners](../using_pheval/plugins_and_runners.md)￼
* Prepare input data and corpora
→ [Utilities](../utilities/index.md)
* Benchmark and analyse results
→ [Benchmarking](../benchmarking/index.md)
* Write your own runner (for developers)
→ [Developing a PhEval plugin](../resources_for_contributors/developing_a_pheval_plugin.md)