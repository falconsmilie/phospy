# PhosPy

PhosPy currently exposes a focused rewrite public surface.

## Package Boundary

```text
src/phospy/         # supported rewrite package
legacy_archive/phospy_legacy/  # migration reference only (not installed package content)
```

## Supported Today

- Dataset builder route: supported
- Kinase workflow route: supported
- Signalome workflow route: first real vertical slice (module assignments, modules, network)

## Rewrite Contract

All public execution entry points use `run(request)`.

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    KinaseWorkflow,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 0.7], "sample_b": [1.2, 0.8]},
            index=["MAPK14;Y182;", "GSK3B;S9;"],
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "GSK3B"],
                "site": ["Y182", "S9"],
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                ],
            },
            index=["MAPK14;Y182;", "GSK3B;S9;"],
        ),
        organism=Organism.RAT,
    )
)

references = ReferenceBundle(
    organism=Organism.RAT,
    kinase_substrate_map=pd.DataFrame(
        {
            "kinase": ["MAP2K6", "MAP2K6"],
            "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
        }
    ),
    site_sequences=pd.DataFrame(
        {
            "site_sequence": dataset.site_metadata.loc[:, "site_sequence"],
        },
        index=pd.Index(dataset.site_metadata.index, name="site_id"),
    ),
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
    )
)

profile_scores = result.scoring_result.profile_scores
pred_mat = result.prediction_result.pred_mat
weighted_activity = result.activity_result.weighted_activity
```

## Current Limits

- `DatasetBuildRequest` accepts pandas `DataFrame` inputs or supported table file paths.
- `ReferencePreset.AUTO` requires `dataset.organism`.
- Bundled references are currently available only for the rat lane.
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` remain valid enum values, but
  they are not bundled runtime lanes in this release. Use an explicit
  `ReferenceBundle` for non-rat organisms.
- Kinase scoring enforces a two-substrate scientific floor
  (`scoring_config.min_substrates >= 2`). Single-site kinase profiles are rejected.
- The supported kinase workflow currently performs profile-based scoring only.
  `result.scoring_result.profile_scores` is the scientific score matrix used by
  downstream workflow stages; `motif_scores`, `combined_scores`, and `weights`
  remain optional compatibility fields and are `None` in this route.
- `SignalomeWorkflowResult.expanded_signalome` remains optional and is currently `None`.

## Examples

- [`examples/dataset_builder_demo.py`](examples/dataset_builder_demo.py)
- [`examples/simple_workflow_demo.py`](examples/simple_workflow_demo.py)
- [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py)

## Docs

- [`docs/api.md`](docs/api.md)
- [`docs/cli.md`](docs/cli.md)
- [`docs/output_bundles.md`](docs/output_bundles.md)
- [`docs/validation.md`](docs/validation.md)
- [`docs/roadmap.md`](docs/roadmap.md)
