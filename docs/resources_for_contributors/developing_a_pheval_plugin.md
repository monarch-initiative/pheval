# Developing a PhEval plugin

This guide explains how to develop a **PhEval plugin** that exposes a **runner** and produces **PhEval standardised results** that can be benchmarked consistently.

## Video walkthrough

If you prefer a guided walkthrough, start here:

[Write your own PhEval runner](https://www.youtube.com/watch?v=GMYzQO4OcfU)

> !!! abstract "**Key takeaways**"
    1. A runner must implement all `PhEvalRunner` methods (`prepare`, `run`, `post_process`).
    2. Your runner must write **standardised result files** with the required columns for the benchmark type.
    3. **Result filenames must match phenopacket filenames** (file stem matching) so PhEval can align outputs to cases.

---

## Standardised result schemas (required)

PhEval benchmarking operates on **standardised result files**.
Each result file must conform exactly to the required schema for the type of prioritisation being produced.

Schemas are **validated** during post-processing.  
Missing or incorrectly named columns will cause validation to fail.

---

### Gene prioritisation results

Each gene result must contain the following columns:

| Column name       | Type       | Description                  |
|------------------|------------|------------------------------|
| `gene_symbol`     | `pl.String`  | Gene symbol                  |
| `gene_identifier` | `pl.String`  | Gene identifier              |
| `score`           | `pl.Float64` | Tool-specific score          |
| `grouping_id`     | `pl.Utf8`    | Optional grouping identifier |

---

### Variant prioritisation results

Each variant result must contain the following columns:

| Column name | Type       | Description |
|------------|------------|-------------|
| `chrom`     | `pl.String`  | Chromosome |
| `start`     | `pl.Int64`   | Start position |
| `end`       | `pl.Int64`   | End position |
| `ref`       | `pl.String`  | Reference allele |
| `alt`       | `pl.String`  | Alternate allele |
| `score`     | `pl.Float64` | Tool-specific score |
| `grouping_id` | `pl.Utf8`    | Optional grouping identifier |

---

### Disease prioritisation results

Each disease result must contain the following columns:

| Column name          | Type       | Description                          |
|---------------------|------------|--------------------------------------|
| `disease_identifier` | `pl.String`  | Disease identifier  |
| `score`              | `pl.Float64` | Tool-specific score                  |

### The `grouping_id` column (optional but important)

`grouping_id` is optional and enables **joint ranking** of entities that should be treated as a single unit without penalty.

Typical examples include:

- Compound heterozygous variants (multiple variants contributing together)
- Grouped variant representations within the same gene
- Polygenic or grouped signals where multiple items should be evaluated jointly

**How to use it**

- Variants in the same group share the same `grouping_id`
- Variants not in any group should each have a unique `grouping_id`

This preserves ranking semantics when benchmarking.

---

## Result file naming (required)

PhEval aligns result files to cases using **filename stem matching**.

> !!! danger "**Rule:**" 
    The **result filename stem must exactly match the phenopacket filename stem**.

Example:

- Phenopacket: `patient_001.json`
- Result filename: `patient_001-exomiser.json`
- Processed result filename passed to PhEval: `patient_001.json`

If the stems do not match, PhEval cannot reliably associate results with
phenopackets, and benchmarking may be incomplete or incorrect.

> !!! tip "**Recommendation:**"
    Always derive result filenames programmatically from the phenopacket stem.
---

## Step-by-step plugin development

PhEval plugins are typically derived from the runner template and standardised tooling.
The recommended approach uses the PhEval runner template, MkDocs, tox, and uv.

The template is available [here](https://github.com/monarch-initiative/pheval-runner-template)

---

### 1. Scaffold a new plugin

Install `cruft` (used to create projects from the template and keep them up to date):

```bash
pip install cruft
```

Create a project using the template:

```bash
cruft create https://github.com/monarch-initiative/pheval-runner-template
```

---

### 2. Environment and dependencies

Install `uv` (if you do not already use it):

```bash
pip install uv
```

Install dependencies and activate the environment:

```bash
uv sync
source .venv/bin/activate
```

Run the test suite to confirm the setup:

```bash
uv run tox
```

> !!! note 
    The template uses `uv` by default, but this is not required.
    You may use any packaging/dependency manager.
    PhEval only requires a valid `pheval.plugins` entry point.

---

### 3. Implement your custom runner

In the generated template, implement your runner in `runner.py` (under `src/`).

At minimum, implement `prepare`, `run`, and `post_process`:

```python
"""Runner."""

from dataclasses import dataclass
from pathlib import Path

from pheval.runners.runner import PhEvalRunner


@dataclass
class CustomRunner(PhEvalRunner):
    """Runner class implementation."""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        """Prepare inputs."""
        print("preparing")

    def run(self):
        """Execute the tool."""
        print("running")

    def post_process(self):
        """Convert raw outputs to PhEval standardised results."""
        print("post processing")
```

---

### 4. Register the runner entry point

The template populates your `pyproject.toml` entry points.
If you rename the runner class or move files, update this accordingly:

```toml
[project.entry-points."pheval.plugins"]
customrunner = "pheval_plugin_example.runner:CustomRunner"
```

> !!! tip
    The module path and class name are case-sensitive.

---

## Tool-specific configuration (config.yaml)

For `pheval run` to execute, the input directory must contain a `config.yaml`:

```yaml
tool:
tool_version:
variant_analysis:
gene_analysis:
disease_analysis:
tool_specific_configuration_options:
```

- `variant_analysis`, `gene_analysis`, `disease_analysis` must be booleans (`true` / `false`)
- `tool_specific_configuration_options` is optional and may include plugin-specific configuration

### Parsing tool-specific configuration (recommended)

Using `pydantic` can simplify parsing:

```python
from pydantic import BaseModel, Field

class CustomisedConfigurations(BaseModel):
    environment: str = Field(...)
```

Then parse in your runner:

```python
config = CustomisedConfigurations.parse_obj(
    self.input_dir_config.tool_specific_configuration_options
)
environment = config.environment
```

---

## Post-processing: generating standardised results

PhEval can handle ranking and writing result files in the correct locations.
Your runner’s post-processing must:

1. Read tool-specific raw outputs
2. Extract the required fields
3. Construct a Polars DataFrame with the required schema
4. Call the appropriate PhEval helper method to write standardised results

### Result generation helpers

!!! warning "Breaking change (v0.5.0)"

    `generate_pheval_result` was replaced with:

    - `generate_gene_result`
    - `generate_variant_result`
    - `generate_disease_result`

#### Generating gene result files

Use `generate_gene_result` to write PhEval-standardised gene results
from a Polars DataFrame.

```python
from pheval.post_processing.post_processing import (
    generate_gene_result,
    SortOrder,
)

generate_gene_result(
    results=pheval_gene_result,      # Polars DataFrame (gene schema)
    sort_order=SortOrder.DESCENDING, # or SortOrder.ASCENDING
    output_dir=output_directory,     # typically self.output_dir
    result_path=result_path,         # path to raw tool output, stem MUST match phenopacket stem exactly
    phenopacket_dir=phenopacket_dir, # directory containing phenopackets
)
```

#### Generating variant result files

Use `generate_variant_result` to write PhEval-standardised variant results.

```python
from pheval.post_processing.post_processing import (
    generate_variant_result,
    SortOrder,
)

generate_variant_result(
    results=pheval_variant_result,   # Polars DataFrame (variant schema)
    sort_order=SortOrder.DESCENDING,
    output_dir=output_directory,
    result_path=result_path,         # stem must match phenopacket stem
    phenopacket_dir=phenopacket_dir,
)
```

#### Generating disease result files

Use `generate_disease_result` to write PhEval-standardised disease results.

```python
from pheval.post_processing.post_processing import (
    generate_disease_result,
    SortOrder,
)

generate_disease_result(
    results=pheval_disease_result,   # Polars DataFrame (disease schema)
    sort_order=SortOrder.DESCENDING,
    output_dir=output_directory,
    result_path=result_path,         # stem must match phenopacket stem
    phenopacket_dir=phenopacket_dir,
)
```

> !!! important 
    The stem of `result_path` must exactly match the phenopacket stem.
    This often requires stripping tool-specific suffixes from raw output filenames.

---

## Adding metadata to results.yml (optional)

PhEval writes a `results.yml` file to the output directory by default.
You can add customised metadata by overriding `construct_meta_data()`.

Example dataclass:

```python
from dataclasses import dataclass

@dataclass
class CustomisedMetaData:
    customised_field: str
```

Runner implementation:

```python
def construct_meta_data(self):
    self.meta_data.tool_specific_configuration_options = CustomisedMetaData(
        customised_field="customised_value"
    )
    return self.meta_data
```

---

## Helper utilities (optional)

PhEval provides helper methods that can simplify runner implementations.

### PhenopacketUtil

Useful for extracting observed phenotypes when tools do not accept phenopackets directly:

::: src.pheval.utils.phenopacket_utils.PhenopacketUtil
    handler: python
    options:
      members:
        - PhenopacketUtil
      show_root_heading: false
      show_source: true

Example usage:

```python
from pheval.utils.phenopacket_utils import phenopacket_reader, PhenopacketUtil

phenopacket = phenopacket_reader("/path/to/phenopacket.json")
phenopacket_util = PhenopacketUtil(phenopacket)

observed_phenotypes = phenopacket_util.observed_phenotypic_features()
observed_phenotypes_hpo_ids = [p.type.id for p in observed_phenotypes]
```

---

## Testing your runner

Install dependencies:

```bash
uv sync
```

Run PhEval using your custom runner:

```bash
pheval run -i ./input_dir -t ./test_data_dir -r customrunner -o output_dir
```

Notes:

- the `-r/--runner` value must match the entry point name (lowercase)
- confirm that standardised result files are produced and validate correctly
- confirm that result file stems match the phenopacket file stems

---

## Checklist before release

- Runner implements `prepare`, `run`, `post_process`
- Entry point registered under `pheval.plugins`
- Standardised results conform to required schema(s)
- Result filenames use **phenopacket stem matching**
- Optional: `grouping_id` correctly set for grouped ranking scenarios
- Optional: `results.yml` metadata populated where useful
