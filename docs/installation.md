# Installation

PhosPy supports Python 3.11 and 3.12.

For a first analysis, install the standard package from PyPI. The development
installation is only needed when you plan to change PhosPy itself.

## Install From PyPI

```bash
pip install phospy
```

For Parquet input and output, install the optional extra:

```bash
pip install "phospy[parquet]"
```

Confirm the installation:

```python
import phospy

print(phospy.__version__)
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

- New to PhosPy? Follow the [Quickstart](quickstart.md).
- Preparing your own tables? Read [Build an Analysis-Ready Dataset](api/dataset-build-workflow.md).
- Choosing an analysis? See [Choose and Follow an Analysis Path](index.md#choose-and-follow-an-analysis-path).
