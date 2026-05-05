# Quickstart

This page gives you one small, tested workflow. Keep it small first; add custom
references, total-protein correction, and larger signalome runs only after this
works.

## 1. Install

```bash
pip install phospy
```

For local development:

```bash
pip install -e ".[dev]"
```

## 2. Prepare Two Tables

`phospho` is numeric, with phosphosite IDs as the index. `site_metadata` uses the
same index and describes each site.

Required `site_metadata` columns for this lane:

- `gene_symbol`
- `site`
- `site_sequence`
- `protein_id` for signalome

A site ID should look like `TSC2;S939;`.

## 3. Run the Python Workflow

```python title="Analysis Ready Dataset, Kinase Workflow, and Signalome Workflow"
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
    SignalomeWorkflow
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SignalomeConfig,
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
        preprocessing_config=DatasetPreprocessingConfig.from_raw_phosphosite_table(),
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig.default(),
        prediction_config=KinasePredictionConfig.deterministic(),
        activity_config=None,
        site_sequence_conflict_policy="prefer_reference",
    )
)

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig.sampled_candidate_scoring(),
    )
)

print(
    "dataset shape",
    dataset.phospho.shape
)
print(
    "prediction shape",
    kinase_result.prediction_result.pred_mat.shape
)
print(
    "signalome modules",
    signalome_result.signalome_modules.table.shape
)

scale_guard = signalome_result.provenance.workflow_parameters["scale_guard"]
print("tree generation mode", scale_guard["tree_generation_mode"])
print("tree generation is approximate", scale_guard["tree_generation_is_approximate"])
print("candidate scoring mode", scale_guard["candidate_scoring_mode"])
print(
    "candidate scoring is approximate",
    scale_guard["candidate_scoring_is_approximate"],
)
```

`candidate_scoring_policy="sampled"` only approximates candidate module-count
scoring. It does not make tree generation approximate.

Why `activity_config=None`? The example has only two sites. The activity stage is
more useful on larger data and defaults to a higher substrate-support threshold.
For real datasets, you can remove that line or configure `KinaseActivityConfig`.

## 4. Run the Same Lane With the CLI

Write `phospho.csv` and `site_metadata.csv` with the first column as the index,
then run:

```bash
phospy signalome \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --skip-activity \
  --outdir ./out
```

Supported table formats are `.csv`, `.tsv`, `.txt` as tab-separated text, and
`.parquet` when optional parquet dependencies are installed.

## 5. If It Fails

| Symptom                                          | Most likely fix                                                                                                      |
|--------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|
| `ReferencePreset.AUTO` cannot resolve references | Use `organism=Organism.RAT` with bundled references, or pass an explicit `ReferenceBundle`.                          |
| Signalome fails on `protein_id`                  | Add a non-empty `protein_id` for every interpreted site. Gene symbols are not used as a protein-identity substitute. |
| Missing-value error                              | Start with a complete matrix, or configure row-median imputation deliberately.                                       |
| Site metadata does not align                     | Make `site_metadata.index` exactly match `phospho.index`.                                                            |
| File input fails                                 | Check that the first CSV/TSV column is the row index and that the suffix is supported.                               |
| The tiny example fails after enabling activity   | Use more sites or lower the activity thresholds deliberately.                                                        |
