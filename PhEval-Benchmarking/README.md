# PheVal Benchmarking

Assesses the percentage of the top, top3 and top5 genes from Exomiser runs - retrieving the data from phenopackets containing the known variant and gene variants. 


## Set up
```

git clone --branch PheVal-Benchmarking https://github.com/monarch-initiative/pheval.git

cd pheval/Pheval-Benchmarking

poetry install

poetry shell

```

## Basic usage

Requires a list of results directories output from Exomiser as well as the full path to the phenopackets used as the input.

```
python assess_prioritisation.py -d /path/to/directory_list -ppackets /path/to/phenopackets/ -o output.txt

```

For more details into the inputs:

```
python assess_prioritisation.py --help

```
