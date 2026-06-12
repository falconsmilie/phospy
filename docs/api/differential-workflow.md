# Differential Workflow

This page explains the public differential workflow API.

## Purpose

`DifferentialAnalysisWorkflow` is a first-class PhosPy workflow entrypoint
exposed from top-level `phospy`. It runs moderated differential analysis from an
`AnalysisReadyPhosphoDataset` plus explicit experimental design and contrast
definitions.

The supported public import route is:

```python
from phospy import DifferentialAnalysisWorkflow
```

`from phospy import DifferentialAnalysis` and
`from phospy.api import DifferentialAnalysis` are not supported public routes.

Differential analysis is separate from kinase and signalome workflows.

```python
differential_result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
)
```

## Imports

```python
from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    SampleDesignRecord,
    TechnicalReplicatePolicy,
)
```

## Input Contract

`DifferentialAnalysisWorkflow.run(...)` accepts `DifferentialAnalysisRequest`.
`DifferentialAnalysisRequest` is a lightweight command payload: constructing it
does not prove the dataset/design/contrast combination is valid. The workflow
validator runs first inside `DifferentialAnalysisWorkflow.run(...)` and rejects
invalid request fields or scientific states before interpretation or
statistical execution.

Required inputs:

- `dataset`: `AnalysisReadyPhosphoDataset`
- `design`: `ExperimentalDesign`
- `contrasts`: `tuple[Contrast, ...]`

`dataset.sample_metadata` may exist, but differential workflow design semantics
are owned by the explicit `design` and `contrasts` request fields.

Optional inputs:

- `config` (`DifferentialAnalysisConfig()` by default), including:
  - `technical_replicate_policy` (`TechnicalReplicatePolicy.REJECT`)
  - `allow_design_subset` (`False`)
  - `minimum_condition_replicates` (`2`)
  - `empirical_bayes` (`EmpiricalBayesConfig()`)

### Matrix Shape Expectations

- `dataset.phospho` must be a numeric phosphosite-by-sample matrix.
- Rows are phosphosite features keyed by `site_key`; columns are sample IDs.
- Values are expected to be finite and analysis-ready (no hidden preprocessing
  is run inside the differential workflow).

### Design Matrix Expectations

- `ExperimentalDesign.samples[].sample_id` must align to dataset sample IDs.
- `condition` labels are required and cannot be empty.
- Duplicate `sample_id` values are rejected.
- By default, all dataset samples must appear in the design.
- With `config.allow_design_subset=True`, design may define a strict sample subset.
- Conditions, replicate identity, batch fields, and block fields are taken from
  `ExperimentalDesign`, not inferred from `dataset.sample_metadata`.
- Declared fixed-effect covariates can be included in the model as ordinary
  fixed terms: `BatchCovariate`, `CategoricalCovariate`, and
  `ContinuousCovariate`.

### Contrast Expectations

- Contrasts are explicit condition-vs-condition definitions (`numerator` and
  `denominator`).
- Contrast conditions must exist in the design.
- Minimum replicate requirements are enforced per condition.
- Invalid contrast definitions fail before statistical execution.

### Sample Alignment Rules

- Alignment is label-based on sample IDs, not positional.
- Reordered sample columns are allowed when labels remain consistent.
- Design/sample mismatch fails with explicit validation errors.

### Intensity Scale and Missing Values

- Differential analysis requires an established log2 phospho intensity scale
  at the `AnalysisReadyPhosphoDataset` boundary before it can emit `logFC`.
- Differential analysis assumes valid upstream preprocessing.
- It does not build datasets, perform sequence resolution, localisation
  filtering, imputation, normalisation, or batch correction.
- Missing-value policy is inherited from the analysis-ready dataset boundary.

## Design and Contrast Requirements

- Conditions are never inferred from sample names.
- Conditions are not inferred from `dataset.sample_metadata` columns.
- `ExperimentalDesign.samples[].sample_id` must align to dataset sample IDs.
- By default, every dataset sample must appear in design.
- Duplicate design sample IDs are rejected.
- Empty condition labels are rejected.
- Contrast conditions must exist in the design.
- Each contrast must satisfy minimum replicate counts for numerator and
  denominator conditions.
- Use at least two biological replicates per condition for meaningful
  differential examples and interpretation.
- Fixed-effect batch, categorical covariate, and continuous covariate terms are
  executable as ordinary design covariates when the resolved design is full
  rank and requested contrasts are estimable.
- Batch fixed effects are not batch correction. They do not implement ComBat,
  RUV, `removeBatchEffect`, or any data-cleaning/removal step.
- Paired/blocking, repeated-measure, `duplicateCorrelation`-style correlated
  replicate, and mixed-effects designs remain rejected unless explicitly
  implemented in a future release.

## Empirical-Bayes Settings

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `method` | `str` | `"standard"` | `"standard"`, `"robust"` | `standard` applies limma-style moderation. `robust` applies winsorized robust hyperparameter estimation. |
| `trend` | `bool` | `False` | `True`, `False` | Enables mean-intensity trend fitting (`limma-trend` style). |
| `winsor_tail_p` | `tuple[float, float]` | `(0.05, 0.1)` | Two probabilities in `[0, 1)` summing to `< 1.0` | Tail clipping parameters used by robust moderation. |

## Public Workflow Shape

`DifferentialAnalysisWorkflow` is a thin public shell. Internal execution
follows the
same stage pattern as other PhosPy workflows:

1. `validated = validator.run(request)`
2. `interpreted = interpreter.run(validated)`
3. `result = executor.run(interpreted)`

The design contract is validated before statistical execution. The interpreter
receives resolved matrix-ready design/contrast structures.

The lower-level science executor owns numeric model computation only. Its
internal output is a stat-only computation payload, not a public scientific
result object. The workflow executor attaches dataset identity metadata from
`AnalysisReadyPhosphoDataset.site_metadata` and only then constructs the public
`DifferentialAnalysisResult`.

## Output Model

`DifferentialAnalysisWorkflow.run(...)` returns `DifferentialAnalysisResult`.

Common outputs include:

- `result.table_for("B_vs_A")` with columns:
  - `site_key`
  - `display_id`
  - `organism`
  - `protein_namespace`
  - `protein_identifier`
  - `gene_symbol`
  - `site`
  - optional workflow-relevant protein metadata such as `protein_id`,
    `protein_accession`, and `isoform_id`
  - `logFC`
  - `t`
  - `P.Value`
  - `adj.P.Val`
- `result.prior_diagnostics`
- `result.mean_variance_trend_diagnostics` (when trend is enabled)

Each contrast result table is indexed by the input `site_key` values. The
`site_key` column must exactly match the index. The minimum public identity
columns are `site_key`, `display_id`, `organism`, `protein_namespace`,
`protein_identifier`, `gene_symbol`, and `site`. Workflow-created results
preserve that required protein-scoped context from `dataset.site_metadata` and
also preserve optional workflow-relevant protein metadata such as `protein_id`
when present. `display_id` remains a human-readable label, and repeated
`display_id` values remain distinct rows when their `site_key` values differ.

For a dataset containing two protein-scoped rows that both display as
`MAPK14;Y182;`, the differential result remains keyed by the two input
`site_key` values:

```python
table = result.table_for("B_vs_A")
duplicate_rows = table.loc[
    table["display_id"] == "MAPK14;Y182;",
    ["site_key", "display_id", "protein_identifier", "logFC", "adj.P.Val"],
]

assert duplicate_rows.index.name == "site_key"
assert duplicate_rows["site_key"].tolist() == duplicate_rows.index.tolist()
assert duplicate_rows["site_key"].is_unique
assert duplicate_rows["display_id"].tolist() == ["MAPK14;Y182;", "MAPK14;Y182;"]
```

The duplicated display label is carried for readability, but statistics are
attached to the analysis-ready `site_key` rows. Differential analysis does not
merge result rows by `display_id`, `gene_symbol`, or `site`.

Direct public `DifferentialAnalysisResult` construction follows the same
identity contract. Display-indexed, stat-only, `GENE;SITE;`-keyed, and
arbitrary non-encoded contrast tables are rejected in public construction.
Result validation does not derive `site_key` from `display_id`, infer protein
identity from `gene_symbol`, decode missing protein context for users, or repair
weak identity metadata.
It is safe to construct directly only when the caller already has complete
public result tables that satisfy this identity-bearing contract.

Stat-only contrast tables belong only to the internal statistical computation
payload. They are not valid `DifferentialAnalysisResult` tables, even when their
row index happens to contain encoded `site_key` values.

### Statistical Method and Multiple Testing

- Effect estimates are OLS contrast estimates (`logFC`).
- Test statistics are moderated t-statistics (`t`) using empirical-Bayes
  variance moderation.
- Raw p-values are reported in `P.Value`.
- Multiple-testing-adjusted p-values are reported in `adj.P.Val`.
- Current multiple-testing policy supports Benjamini-Hochberg adjustment.

### Result Structure and Metadata

- Per-contrast result tables (`result.table_for(contrast_name)`)
- Prior diagnostics (`result.prior_diagnostics`)
- Method metadata:
  - `result.empirical_bayes_method`
  - `result.empirical_bayes_robust`
  - `result.empirical_bayes_trend`
- Trend diagnostics when enabled (`result.mean_variance_trend_diagnostics`)
- Structured policy provenance (`result.policy_provenance`) including:
  - design formula/description, condition columns, fixed-effect covariate
    columns, and fixed-effect covariate kinds
  - rank and contrast-estimability validation status
  - typed contrast definitions and contrast vectors
  - replicate and technical-replicate policy details
  - empirical-Bayes settings
  - p-value and adjusted p-value methods
  - missing-value handling policy
  - unsupported-design rejection policy and intentionally rejected unsupported
    design features

## Limitations and Non-goals

- Assumes valid upstream preprocessing and quantitative inputs.
- Assumes valid design matrix, contrast definitions, and replicate structure.
- Does not resolve peptide/site ambiguity.
- Does not perform localisation filtering unless explicitly implemented upstream.
- Does not perform missing-value imputation unless explicitly implemented
  upstream.
- Does not perform batch correction unless explicitly implemented upstream.
  A batch term in `ExperimentalDesign.fixed_effects` is a fixed model covariate,
  not ComBat, RUV, `removeBatchEffect`, `duplicateCorrelation`, or mixed-effects
  modelling.
- Statistical interpretation depends on design, contrast specification, and
  replicate structure.

## Worked Example

This site-level example sets localisation policy at dataset build time so
low-confidence phosphosite assignments fail fast before differential statistics.
That avoids reporting `logFC` on ambiguously localised sites.

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    IntensityScaleKind,
    Organism,
    SampleDesignRecord,
)

phospho = pd.DataFrame(
    {
        "control_rep1": [8200.0, 9100.0],
        "control_rep2": [8000.0, 9000.0],
        "treatment_rep1": [16200.0, 9150.0],
        "treatment_rep2": [15800.0, 9050.0],
    },
    index=["MAPK14;Y182;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3B"],
        "site": ["Y182", "S9"],
        "site_sequence": [
            "MPRKSLVGTPYWMNQYAVNQKQTLRDLKQEN",
            "ATMSGRPRTTSFAESSKPVQQPSAFGQAAAL",
        ],
        "display_id": ["MAPK14;Y182;", "GSK3B;S9;"],
        "organism": ["rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "GSK3B"],
        "protein_id": ["MAPK14", "GSK3B"],
        "localisation_confidence": [0.95, 0.92],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            ),
        ),
    )
)

assert dataset.intensity_scale_state.kind is IntensityScaleKind.LOG2
assert dataset.intensity_scale_state.is_established

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord(
            sample_id="control_rep1",
            condition="control",
            biological_replicate_id="control_r1",
        ),
        SampleDesignRecord(
            sample_id="control_rep2",
            condition="control",
            biological_replicate_id="control_r2",
        ),
        SampleDesignRecord(
            sample_id="treatment_rep1",
            condition="treatment",
            biological_replicate_id="treatment_r1",
        ),
        SampleDesignRecord(
            sample_id="treatment_rep2",
            condition="treatment",
            biological_replicate_id="treatment_r2",
        ),
    )
)
contrasts = (
    Contrast(
        name="treatment_vs_control",
        numerator_condition="treatment",
        denominator_condition="control",
    ),
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
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
            "protein_id",
        ],
    ]
)
print(result.table_for("treatment_vs_control").loc[:, ["logFC", "adj.P.Val"]])
```

Interpretation notes for this tiny synthetic matrix:

- The `MAPK14;Y182;` display label has higher treatment intensity than control,
  so `logFC` should be positive for `treatment_vs_control`.
- The `GSK3B;S9;` display label is approximately unchanged across conditions,
  so its `logFC` should be near zero.
- The example demonstrates workflow contracts and mechanics, not biological
  discovery or study-level statistical power.
