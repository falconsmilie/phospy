# PhosPy

PhosPy is a Python package for selected phosphoproteomics workflows inspired by
PhosR. It is aimed at scientists who want a clear Python lane from phosphosite
intensity tables to kinase scoring, kinase prediction, and optional signalome
analysis.

PhosPy does **not** provide HTTP endpoints or a web service. The supported user
interfaces are:

- the Python API
- the `phospy` command-line interface

## Install

PhosPy requires Python 3.10 or newer.

```bash
pip install phospy
```

For `.parquet` input or output support:

```bash
pip install "phospy[parquet]"
```

For local development from a clone:

```bash
pip install -e ".[dev]"
pyright
pytest -m "not parity"
```

For reproducible scientific/regression runs aligned to CI:

```bash
pip install -c constraints/ci.txt -e ".[dev,test]"
pytest tests/parity -m parity -s
```

## Beginner Lane

The smallest supported lane is:

1. build a dataset with `organism=Organism.RAT`
2. run kinase with `references=ReferencePreset.AUTO`
3. run signalome only when `site_metadata.protein_id` is present

Bundled runtime references in `1.5.0` are rat-only. For human or mouse work,
create and pass an explicit `ReferenceBundle` in Python instead of using
`ReferencePreset.AUTO`.

## Minimum Input Shape

`phospho` is a numeric site-by-sample table. Its index should use standard
PhosPy site IDs such as `TSC2;S939;`.

`site_metadata` must align to `phospho.index` and include:

- `gene_symbol`
- `site`
- `site_sequence` when you want motif-based scoring
- `protein_id` when you want signalome analysis

## Minimal Python Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SignalomeWorkflowRequest,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70],
        "sample_b": [1.10, 0.80],
        "sample_c": [0.95, 0.75],
    },
    index=["TSC2;S939;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "_______MSGRPRTTSFAESCKPVQQPSAFG",
        ],
        "protein_id": ["TSC2", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        activity_config=None,  # keep this tiny two-site example in the safe lane
    )
)

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)

print(dataset.phospho.shape)
print(kinase_result.prediction_result.pred_mat.round(4))
print(signalome_result.signalome_modules.table)
```

With the tables above you should get a strict `AnalysisReadyPhosphoDataset`, a
non-empty kinase prediction matrix, and signalome tables because every site has
an explicit `protein_id`.

## Minimal CLI Example

```bash
phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --skip-activity \
  --outdir ./out
```

The CLI writes a dataset directory, a kinase directory, and a short list of file
paths written.

## Import Contract

Use top-level `phospy` for the four main entrypoints only:

```python
from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy import KinaseWorkflow, SignalomeWorkflow
```

Use `phospy.api` for requests, configs, results, enums, references, and public
exceptions.

## Documentation

1. [Quickstart](https://phospy.com/docs/quickstart/)
2. [API Guide](https://phospy.com/docs/api/)
3. [CLI Guide](https://phospy.com/docs/cli/)
4. [Validation Guide](https://phospy.com/docs/validation/)

## Citation

If you use PhosPy in scientific work, cite this software release using
[`CITATION.cff`](CITATION.cff) and also cite the upstream PhosR project and
publications described in [`NOTICE.md`](NOTICE.md).
