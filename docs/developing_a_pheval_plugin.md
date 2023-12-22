# Developing a PhEval Plugin

## Description

Plugin development allows PhEval to be extensible, as we have designed it.
The plugin goal is to be flexible through custom runner implementations. This plugin development enhances the PhEval functionality. You can build one quickly using this step-by-step process.

==_All custom Runners implementations must implement all_ **PhevalRunner** _methods_==

::: src.pheval.runners.runner.PhEvalRunner
    handler: python
    options:
      members:
        - PhEvalRunner
      show_root_heading: false
      show_source: true
---

## Step-by-Step Plugin Development Process

The plugin structure is derived from a [cookiecutter](https://cookiecutter.readthedocs.io/en/stable/) template, [Sphintoxetry-cookiecutter](https://github.com/hrshdhgd/sphintoxetry-cookiecutter), and it uses [Sphinx](https://www.sphinx-doc.org/en/master/), [tox](https://tox.wiki/en/latest/) and [poetry](https://python-poetry.org) as core dependencies.
This allows PhEval extensibility to be standardized in terms of documentation and dependency management.

### 1. Sphintoxetry-cookiecutter scaffold

First, install the cruft package. Cruft enables keeping projects up-to-date with future updates made to this original template.

Install cruft from pip

```bash
pip install cruft
```

Next, create a project using the sphintoxetry-cookiecutter template.

```
cruft create https://github.com/monarch-initiative/monarch-project-template
```

### 2. Further setup

#### Install poetry if you haven't already.

```
pip install poetry
```

#### Install dependencies

```
poetry install
```

#### Add PhEval dependency

```
poetry add git+https://github.com/monarch-initiative/pheval.git
```

#### Run tox to see if the setup works

```
poetry run tox
```

### 3. Implement PhEval Custom Runner

_The runner name is arbitrary and custom Runner name was chose by demonstrative purposes_

Create a runner file inside the plugin project, e.g:

![a](./imgs/plugin_example_folder_structure.png)

```python
"""Custom Pheval Runner."""
from dataclasses import dataclass
from pathlib import Path
from pheval.runners.runner import PhEvalRunner


@dataclass
class CustomPhevalRunner(PhEvalRunner):
    """CustomPhevalRunner Class."""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        """prepare method."""
        print("preparing")

    def run(self):
        """run method."""
        print("running with custom pheval runner")

    def post_process(self):
        """post_process method."""
        print("post processing")

```

### 4. Add PhEval Plugins section to the pyproject.toml file

```toml
[tool.poetry.plugins."pheval.plugins"]
customrunner = "pheval_plugin_example.runner:CustomPhevalRunner"
```

==Replace the value above with the path to your custom runner plugin==

### 5. Implementing PhEval helper methods

Streamlining the creation of your custom PhEval runner can be facilitated by leveraging PhEval's versatile helper methods, where applicable.

Within PhEval, numerous public methods have been designed to assist in your runner methods. The utilisation of these helper methods is optional, yet they are crafted to enhance the overall implementation process.

#### Utility methods

The `PhenopacketUtil` class is designed to aid in the collection of specific data from a Phenopacket.

::: src.pheval.utils.phenopacket_utils.PhenopacketUtil
    handler: python
    options:
      members:
        - PhenopacketUtil
      show_root_heading: false
      show_source: true
---

`PhenopacketUtil` proves particularly beneficial in scenarios where the tool for which you're crafting a runner implementation does not directly accept Phenopackets as inputs. Instead, it might require elements—such as HPO IDs— via the command-line interface (CLI). In this context, leveraging PhenopacketUtil within the runner's preparation phase enables the extraction of observed phenotypic features from the Phenopacket input, facilitating seamless processing.

An example of how this could be implemented is outlined here:

```python
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil

phenopacket = phenopacket_reader("/path/to/phenopacket.json")
phenopacket_util = PhenopacketUtil(phenopacket)
# To return a list of all observed phenotypes for a phenopacket
observed_phenotypes = phenopacket_util.observed_phenotypic_features()
# To extract just the HPO ID as a list
observed_phenotypes_hpo_ids = [
    observed_phenotype.id for observed_phenotype in observed_phenotypes
]
```
#### Additional tool-specific configurations

For the `pheval run` command to execute successfully, a `config.yaml` should be found within the input directory supplied on the CLI.

```yaml
tool: 
tool_version: 
variant_analysis: 
gene_analysis: 
disease_analysis: 
tool_specific_configuration_options:
```

The `tool_specific_configuration_options` is an optional field that can be populated with any variables specific to your runner implementation that is required for the running of your tool.

All other fields are required to be filled in. The `variant_analysis`, `gene_analysis`, and `disease_analysis` are set as booleans and are for specifying what type of analysis/prioritisation the tool outputs.

To populate the `tool_specific_configurations_options` with customised data, we suggest using the `pydantic` package as it can easily parse the data from the yaml structure.

e.g.,

_Define a `BaseModel` class with the fields that will populate the `tool_specific_configuration_options`_

```python
from pydantic import BaseModel, Field

class CustomisedConfigurations(BaseModel):
    """
    Class for defining the customised configurations in tool_specific_configurations field,
    within the input_dir config.yaml
    Args:
        environment (str): Environment to run
    """
    environment: str = Field(...)
```

_Within your runner parse the field into an object._

```python
from dataclasses import dataclass
from pheval.runners.runner import PhEvalRunner
from pathlib import Path

@dataclass
class CustomPhevalRunner(PhEvalRunner):
    """CustomPhevalRunner Class."""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        """prepare method."""
        print("preparing")
        config = CustomisedConfigurations.parse_obj(
            self.input_dir_config.tool_specific_configuration_options
        )
        environment = config.environment
        
    def run(self):
        """run method."""
        print("running with custom pheval runner")

    def post_process(self):
        """post_process method."""
        print("post processing")
        

```

#### Post-processing methods

PhEval currently supports the benchmarking of gene, variant, and disease prioritisation results. 

To benchmark these result types, PhEval TSV result files need to be generated. 

PhEval can deal with the ranking and generation of these files to the correct location. However, the runner implementation must handle the extraction of essential data from the tool-specific raw results. This involves transforming them into a list comprising PhEval data classes, with each instance representing a result entry.

The dataclasses representing essential information extracted from tool-specific output for gene, variant, and disease prioritisation are defined as follows:

::: src.pheval.post_processing.post_processing.PhEvalGeneResult
    handler: python
    options:
      members:
        - PhEvalGeneResult
      show_root_heading: false
      show_source: true
---

::: src.pheval.post_processing.post_processing.PhEvalVariantResult
    handler: python
    options:
      members:
        - PhEvalVariantResult
      show_root_heading: false
      show_source: true
---

::: src.pheval.post_processing.post_processing.PhEvalDiseaseResult
    handler: python
    options:
      members:
        - PhEvalDiseaseResult
      show_root_heading: false
      show_source: true
---

The `generate_pheval_result()` can be implemented in your runner to write out the PhEval TSV results.

An example of how the method can be called is outlined here:

```python
from pheval.post_processing.post_processing import generate_pheval_result

generate_pheval_result(
    pheval_result=pheval_gene_result, # this is the list of extracted PhEval result requirements
    sort_order_str="descending", # or can be ascending - this determines in which order the scores will be ranked
    output_dir=output_directory, # this can be accessed from the runner instance e.g., self.output_dir
    tool_result_path=tool_result_json # this is the path to the tool-specific raw results file
)
```

#### Adding metadata to the results.yml

By default, PhEval will write a `results.yml` to the output directory supplied on the CLI. 

The `results.yml` contains basic metadata regarding the run configuration, however, there is also the option to add customised run metadata to the `results.yml` in the `tool_specific_configuration_options` field.

To achieve this, you'll need to create a `construct_meta_data()` method within your runner implementation. This method is responsible for appending customised metadata to the metadata object in the form of a defined dataclass. It should return the entire metadata object once the addition is completed.

e.g.,

_Defined customised metadata dataclass:_

```python
from dataclasses import dataclass

@dataclass
class CustomisedMetaData:
    customised_field: str
```

_Example of implementation in the runner._

```python
from dataclasses import dataclass
from pheval.runners.runner import PhEvalRunner
from pathlib import Path

@dataclass
class CustomPhevalRunner(PhEvalRunner):
    """CustomPhevalRunner Class."""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        """prepare method."""
        print("preparing")

    def run(self):
        """run method."""
        print("running with custom pheval runner")

    def post_process(self):
        """post_process method."""
        print("post processing")
        
    def construct_meta_data(self):
        """Add metadata."""
        self.meta_data.tool_specific_configuration_options = CustomisedMetaData(customised_field="customised_value")
        return self.meta_data

```

### 6. Test it.

To update your custom pheval runner implementation, you must first install the package

```
poetry install
```

Now you have to be able to run PhEval passing your custom runner as parameter. e.g.,

```
pheval run -i ./input_dir -t ./test_data_dir -r 'customphevalrunner' -o output_dir
```

The `-r` parameter stands for your plugin runner class name, and it must be entirely lowercase.

Output:

```
preparing
running with custom pheval Runner
post processing
```

Pay attention to "_==running with custom pheval Runner==_" line, this is exactly what we had implemented in the **CustomPhevalRunner** Example
