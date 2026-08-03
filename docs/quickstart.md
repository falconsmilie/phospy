# Quick Start

This page gives you one small, tested workflow. Keep it small first; add custom
references, total-protein correction, and larger signalome runs only after this
works.

## 1. Installation

```bash
pip install phospy
```

If this is your first Python project, create and activate a virtual
environment first:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1
```

For local development:

```bash
pip install -e ".[dev]"
```

## 2. Prepare Two Tables

`phospho` is numeric. Builder input may use display labels such as
`TSC2;S939;` as the index when `site_metadata` provides enough protein context
to derive `site_key`. The built analysis-ready dataset uses `site_key` as its
row index and preserves the display label in `display_id`.

Required `site_metadata` columns for this lane:

- `gene_symbol`
- `site`
- `site_sequence`
- `localisation_confidence`
- `organism`
- `protein_namespace`
- `protein_identifier`
- `protein_group_id` for signalome grouping

Builder input may omit `site_key` only when those protein-context fields are
available for deterministic derivation. The direct
`AnalysisReadyPhosphoDataset` constructor raises immediately. Advanced
trusted reconstruction of already prepared analysis-ready tables must use
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with typed evidence or
explicit waivers for identity, intensity scale, quantitative meaning,
localisation, sequence, and reference context, plus non-waivable aligned-table
structure evidence. Supplied trusted provenance must match the actual table
fingerprints. A display label should look like `TSC2;S939;`; it is not
unique row identity and may repeat when distinct `site_key` rows preserve
distinct protein context. Rows that resolve to the
same `site_key` are duplicate scientific evidence for the same analysis-ready
site and fail by default unless you intentionally choose and audit a non-error
duplicate-site preprocessing policy.

## 3. Run the Python Workflow

```python title="Analysis Ready Dataset, Kinase Workflow, and Signalome Workflow"
from dataclasses import replace

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
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceContextCompatibilityPolicy,
    ReferencePreset,
    SignalomeConfig,
    SignalomeValidationConfig,
    SignalomeWorkflowRequest,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70, 0.85],
        "sample_b": [1.10, 0.80, 0.88],
        "sample_c": [0.95, 0.75, 0.92],
    },
    index=["TSC2;S939;", "GSK3B;S9;", "MAPK14;Y182;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B", "MAPK14"],
        "site": ["S939", "S9", "Y182"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "ATMSGRPRTTSFAESSSPVQQPSAFGQAAAL",
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
        ],
        "display_id": ["TSC2;S939;", "GSK3B;S9;", "MAPK14;Y182;"],
        "organism": ["rat", "rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["TSC2", "GSK3B", "MAPK14"],
        "protein_group_id": ["TSC2", "GSK3B", "MAPK14"],
        "localisation_confidence": [0.96, 0.93, 0.95],
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
        scoring_config=KinaseScoringConfig(
            reliability_profile=KinaseReliabilityProfile.CUSTOM,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            )
        ),
        prediction_config=KinasePredictionConfig.deterministic(),
        activity_config=None,
        site_sequence_conflict_policy="prefer_reference",
    )
)

signalome_config = SignalomeConfig.compatibility()
signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=replace(
            signalome_config,
            validation=replace(
                signalome_config.validation,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                )
            ),
        ),
    )
)

print(
    "dataset shape",
    dataset.phospho.shape
)
print(
    dataset.site_metadata.loc[
        :,
        [
            "site_key",
            "display_id",
            "gene_symbol",
            "site",
            "protein_namespace",
            "protein_identifier",
            "protein_group_id",
        ],
    ]
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

print(
    "tree generation mode", 
    scale_guard["tree_generation_mode"]
)
print(
    "tree generation is approximate", 
    scale_guard["tree_generation_is_approximate"]
)
print(
    "candidate scoring mode", 
    scale_guard["candidate_scoring_mode"]
)
print(
    "candidate scoring is approximate",
    scale_guard["candidate_scoring_is_approximate"],
)
```

This tiny example opts into `SignalomeConfig.compatibility()` because it has
fewer than the five retained sites required by the production network
paired-observation rule. Use `SignalomeConfig.production()` for real runs once
your retained site set is large enough.

Why `activity_config=None`? Activity execution is opt-in and this example has
only two sites. The activity stage is more useful on larger data and defaults
to a higher substrate-support threshold. For real datasets, provide an explicit
`KinaseActivityConfig` when you want activity-like score outputs.

Why no `input_intensity_scale` declaration in this example? The
`from_raw_phosphosite_table()` preset applies a log2 transform, so the builder
can establish intensity scale from the configured preprocessing path.

Supported file-backed table formats for API-driven loading are `.csv`, `.tsv`,
`.txt` as tab-separated text, and `.parquet` when optional parquet dependencies
are installed.

## 4. If It Fails

| Symptom                                          | Most likely fix                                                                                                      |
|--------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|
| `ReferencePreset.AUTO` cannot resolve references | Use `organism=Organism.RAT` with bundled references, or pass an explicit `ReferenceBundle`.                          |
| Signalome fails on `protein_group_id`            | Add a non-empty signalome grouping `protein_group_id` for every interpreted site. Legacy `protein_id` is accepted only as a migration alias. Keep core protein identity in `protein_namespace` and `protein_identifier`; do not use gene symbols or `display_id` labels as fallbacks. |
| Missing-value error                              | Start with a complete matrix, or configure row-median imputation deliberately.                                       |
| Site metadata does not align                     | For builder input, make `site_metadata.index` match `phospho.index`; for trusted `from_trusted_tables(...)` reconstruction, use matching `site_key` indexes. |
| File input fails                                 | Check that the first CSV/TSV column is the row index and that the suffix is supported.                               |
| The tiny example fails after enabling activity   | Use more sites or lower the activity thresholds deliberately.                                                        |
