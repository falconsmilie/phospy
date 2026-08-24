# Installation

PhosPy supports Python 3.11 and 3.12.

## Install From PyPI

```bash
pip install phospy
```

For Parquet input and output, install the optional extra:

```bash
pip install "phospy[parquet]"
```

Confirm the installation:

```bash
python -c "import phospy; print(phospy.__version__)"
```

## Install for Development

From a local clone:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
```

For documentation maintenance only, install the docs extra and run the
standalone strict build:

```bash
pip install -c constraints/ci.txt -e ".[docs]"
make docs-build
```

## Continue

- [Prepare a dataset](api/dataset-build-workflow.md)
- [Run your first analysis](quickstart.md)
- [Choose a workflow](index.md#choose-a-workflow)
