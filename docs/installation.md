# Installation

PhosPy requires Python 3.11 or 3.12.

## Install from PyPI

```bash
pip install phospy
```

Install the Parquet extra if you need `.parquet` input or output:

```bash
pip install "phospy[parquet]"
```

## Development install

From a local clone:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet,docs]"
```

The project-supported documentation build is:

```bash
mkdocs build --strict
```

## Next steps

- [Prepare a dataset](api/dataset-build-workflow.md)
- [Run your first analysis](quickstart.md)
- [Choose a workflow](index.md#workflow-map)
