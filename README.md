# Contact Matrix Generator

This Python 3.9 package can generate C-alpha contact matrices by parsing an
input mmCIF file, and writing out matrices per PDB chains in .npy format.

## Quick start

These instructions will get you a copy of the project up and running on your
local machine for development and testing purposes.

## Installing

Get the repository

```shell
git clone https://github.com/mvaradi/contact-matrix-generator
cd contact-matrix-generator
```

Install the dependencies defined in [requirements.txt](requirements.txt)

```shell
pip install -r requirements.txt
```

The package relies on the following:
* numpy - For calculations and data structures
* gemmi - For parsing mmCIF files
* pytest - For testing
* codecov - For testing coverage
* pytest-cov - For testing coverage

## Usage

The tool can be used as a CLI:

```shell
python .\main.py -i ./input.cif -o /path/to/output/files
```

OR

```shell
python .\main.py --input_path=./input.cif --output_path=/path/to/output/files
```

## Development

The repository has a `pre-commit` configuration that performs a number of
check on the code, before allowing a commit.

### Using pre-commit

```shell
pip install
pre-commit
pre-commit install
```
