# Plugins and runners

This page defines **how execution works in PhEval**.

It is the single authoritative description of plugins, runners, and how phenotype-driven tools are run under the PhEval framework.

## Execution model

PhEval provides a general execution framework and command-line interface.
It does not embed logic for running individual tools.

Instead:

- **Plugins** extend PhEval with tool-specific functionality
- **Runners** (provided by plugins) implement complete execution workflows

PhEval is responsible for orchestration.
Runners are responsible for execution semantics.

## Plugins

A plugin is a Python package that integrates a specific tool with PhEval.

A plugin typically:

- Depends on PhEval
- Wraps a specific phenotype-driven tool
- Registers one or more runners
- Documents tool-specific configuration and usage

Plugins are installed separately from PhEval.
Once installed, their runners are automatically discovered.

## Runners

A runner is the unit of execution within PhEval.

Each runner is responsible for the full workflow for a given tool and configuration:

1. Preparing inputs (phenopackets, variants, resources)
2. Executing the tool
3. Post-processing raw outputs into the PhEval standardised results format

Runners encapsulate tool-specific assumptions while conforming to shared PhEval interfaces.

## Running a runner

All tool execution is performed via the PhEval CLI.

The general pattern is:

    pheval run --runner <runner_name> [options]

Where:

- `<runner_name>` identifies a runner exposed by an installed plugin

PhEval manages discovery and orchestration.
The runner controls execution logic and output generation.

## Multiple runners

A plugin may expose multiple runners.

This allows support for:

- Different tool modes
- Different input types
- Alternative workflows or configurations

Each runner is invoked explicitly by name.

## Outputs and standardisation

All runners are expected to produce outputs in the **PhEval standardised results format**.

This standardisation allows results to be:

- Benchmarked using shared metrics
- Compared across tools and versions
- Analysed using common plotting and reporting utilities

This is a core design principle of PhEval.

## Where execution details live

Tool-specific instructions do not live in the main PhEval documentation.

They are documented in:

- The plugin README
- Tool-specific configuration guides
- Examples provided by plugin authors

This avoids duplication and keeps framework documentation focused and stable.

