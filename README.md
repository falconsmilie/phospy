# PhosPy

PhosPy currently exposes a focused rewrite public surface.

## Package Boundary

```text
src/phospy/         # supported rewrite package
legacy_archive/phospy_legacy/  # migration reference only (not installed package content)
```

## Supported Today

- Dataset builder route: supported
- Kinase workflow route: supported (including activity-stage outputs when enabled)
- Signalome workflow route: first real vertical slice (module assignments, modules, network)

## Rewrite Contract

All public execution entry points use `run(request)`.

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseActivityConfig,
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
        activity_config=KinaseActivityConfig(
            enabled=True,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=20,
        ),
    )
)

profile_scores = result.scoring_result.profile_scores
pred_mat = result.prediction_result.pred_mat
if result.activity_result is not None:
    weighted_activity = result.activity_result.weighted_activity
```

## Current Limits

- `DatasetBuildRequest` accepts pandas `DataFrame` inputs or supported table file paths.
- Builder convention handling is explicit: supported site-metadata aliases are
  `gene`/`gene_name` -> `gene_symbol`, `residue`/`phosphosite`/`site_position` -> `site`,
  and `centralized_sequence` -> `site_sequence`.
- Ambiguous legacy names like `sequence` and `protein` are not auto-guessed; rename
  them to explicit supported fields (`site_sequence`, `protein_id`) for reliable ingestion.
- `ReferencePreset.AUTO` requires `dataset.organism`.
- Bundled references are currently available only for the rat lane.
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` remain valid enum values, but
  they are not bundled runtime lanes in this release. Use an explicit
  `ReferenceBundle` for non-rat organisms.
- Kinase scoring enforces a two-substrate scientific floor
  (`scoring_config.min_substrates >= 2`). Single-site kinase profiles are rejected.
- The supported kinase scoring route now emits profile, motif, and
  profile/motif-combined score tables (`profile_scores`, `motif_scores`,
  `combined_scores`, `weights`) for scientific audit and parity tracking.
- Prediction in this release remains profile-driven in the workflow executor,
  with motif/combined outputs carried in the scoring result contract.
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
- [`docs/architecture/activity_science_port_review.md`](docs/architecture/activity_science_port_review.md)
