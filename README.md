# PhosPy

PhosPy is a Python package for selected phosphoproteomics workflows inspired by
PhosR. It is aimed at scientists who want a clear Python lane from phosphosite
intensity tables to differential phosphorylation analysis, offline
over-representation enrichment, kinase scoring and prediction, and optional
signalome analysis. Enrichment support is offline ORA over caller-supplied
selected identifiers, local set collections, and an explicit background
universe.

"PhosR-inspired" in PhosPy docs means scoped, feature-level comparison lanes. It
does not imply full PhosR package parity or full PhosR API compatibility.

PhosPy does **not** provide HTTP endpoints or a web service. The supported user
interface is the Python API.

## Recommended Reading
You can view the full documentation here: [PhosPy Docs](https://phospy.com/docs)

## Installation

PhosPy requires Python 3.10, 3.11, or 3.12.

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
python scripts/run_pyright.py
pytest -m "not parity"
```

For reproducible scientific/regression runs aligned to CI:

```bash
pip install -c constraints/ci.txt -e ".[dev,test]"
pytest tests/parity -m parity -s
```

For full release-gate validation (unit/integration, reproducibility goldens,
parity, and performance):

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
make test-release-gate
```

## Quick Start

1. Build an analysis-ready phosphoproteomics dataset.
2. Run the supported workflow lane you need: differential, enrichment, kinase,
   or signalome.
3. Explore full API workflow documentation:
   - [Dataset building](docs/api/dataset-build-workflow.md)
   - [Differential workflow](docs/api/differential-analysis.md)
   - [Enrichment workflow](docs/api/enrichment.md)
   - [Kinase workflow](docs/api/kinase.md)
   - [Signalome workflow](docs/api/signalome.md)

Bundled runtime references in the current release are rat-only. For human or
mouse work, create and pass an explicit `ReferenceBundle` in Python instead of
using `ReferencePreset.AUTO`.

The default kinase `scoring_mode="phosr_rank_weighted"` is PhosR-inspired
rank-weighted scoring implemented by PhosPy. It combines available profile and
motif support under PhosPy's support rules; it is not an exact PhosR
implementation and is not intended to provide numerical parity with PhosR.
Optional kinase activity outputs are exploratory kinase activity scores or
activity-like substrate summaries. They depend on substrate coverage and
reference evidence; they are not direct proof of kinase activation or causal
pathway activity.

Enrichment ORA results are overlap statistics under the caller-supplied
background universe. They do not prove pathway activation, regulation, or
biological causality, and PhosPy does not imply GSEA or PTM-SEA support.

Scientific scope categories and parity/open-gap status are maintained in
[`docs/scientific-coverage.md`](docs/scientific-coverage.md). Parity fixture
evidence lives in [`docs/parity.md`](docs/parity.md).
Future coverage direction is tracked in
[`ADR-0025`](docs/adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md);
that roadmap is not a current feature-support claim.
Native SPS/RUV-style batch-correction prerequisites are recorded in
[`ADR-0029`](docs/adr/adr_0029_native_sps_ruv_style_batch_correction_prerequisites.md);
those prerequisites are implemented for the supported
`SpsRuvBatchCorrectionConfig` preprocessing lane. This is a native PhosPy
implementation, not PhosR-equivalent SPS/RUV-III parity.

The `batch_correction` preprocessing group exposes two explicit lanes:
`linear_residualize_batch`, a limited fixed-effect residualisation step, and
native SPS/RUV-style correction through `SpsRuvBatchCorrectionConfig`.
`linear_residualize_batch` preserves condition effects by design and rejects
confounded batch/condition metadata; it is not ComBat, not RUV, not limma
`removeBatchEffect` parity, not native SPS/RUV-style correction, not
PhosR-equivalent batch correction, and not mixed-effects modelling. Native
SPS/RUV-style correction is executable only through explicit structured
preprocessing config and requires caller-supplied controls, batch and protected
condition metadata, a missingness policy, diagnostics, and provenance. RUV-III
style correction is not executable unless a future feature implements
replicate-aware semantics. Any `ruv_readiness` diagnostics are report-only
readiness signals and do not apply correction.

## Kinase Workflow Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    IntensityScaleKind,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

# Tiny synthetic example for workflow mechanics only (not biological discovery).
phospho = pd.DataFrame(
    {
        "control_rep1": [8200.0, 9100.0, 6000.0],
        "control_rep2": [8000.0, 9000.0, 5900.0],
        "treatment_rep1": [16200.0, 9150.0, 13000.0],
        "treatment_rep2": [15800.0, 9050.0, 12800.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A", "TSC2"],
        "site": ["Y182", "S21", "S939"],
        "site_sequence": [
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
        ],
        "display_id": ["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
        "organism": ["rat", "rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        # Signalome has a separate explicit protein_id requirement.
        "protein_id": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index.copy(),
)
sample_metadata = pd.DataFrame(
    {
        "comparison_group": ["control", "control", "treatment", "treatment"],
    },
    index=phospho.columns.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LINEAR,
        preprocessing_config=DatasetPreprocessingConfig(
            # Site-level workflows should fail fast when localisation is missing
            # or below threshold, because ambiguous site assignment can
            # mis-state kinase/substrate interpretation.
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
)

# Builder input may be display-indexed when enough protein context is present.
# The analysis-ready dataset itself is indexed by site_key; direct
# AnalysisReadyPhosphoDataset construction must already use site_key indexes.
print(
    dataset.site_metadata.loc[
        :,
        [
            "site_key",
            "display_id",
            "gene_symbol",
            "site",
            "organism",
            "protein_namespace",
            "protein_identifier",
            "protein_id",
            "site_sequence",
        ],
    ]
)
# sample_metadata is descriptive/alignment metadata on the dataset.
# Differential workflow design is provided separately via ExperimentalDesign.

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        # Safe in this example because organism=rat and bundled runtime
        # references in this release are rat-only.
        references=ReferencePreset.AUTO,
        activity_config=None,
    )
)

print(kinase_result.prediction_result.pred_mat.round(3).iloc[:3, :5])
if kinase_result.prediction_result.substrate_list is not None:
    print(kinase_result.prediction_result.substrate_list.head(5))
```

`site_key` is the true analysis-ready phosphosite row identity. `display_id` is
only a human-readable label and may repeat when different `site_key` values
preserve distinct protein context. Rows that resolve to the same `site_key` are
a scientific ambiguity; the builder fails by default, and any non-error
duplicate-site policy should be chosen deliberately and audited in the
preprocessing report.

Differential result tables use strict protein-scoped identity. Public
`DifferentialAnalysisResult` tables must be indexed by encoded `site_key` values
and include `site_key`, `display_id`, `gene_symbol`, and `site`. Workflow-created
results preserve available protein context such as `organism`,
`protein_namespace`, `protein_identifier`, and `protein_id`. Display-indexed or
stat-only result tables are not valid public inputs.

## Import Contract

Use top-level `phospy` for the dataset, differential, kinase, and signalome
convenience entrypoints:

```python
from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
```

Use `phospy.api` for requests, configs, results, enums, references, and public
exceptions. `EnrichmentWorkflow` is a supported public workflow from
`phospy.api` and `phospy.workflows`:

```python
from phospy.api import EnrichmentConfig, EnrichmentWorkflow
from phospy.api import EnrichmentWorkflowRequest, GeneSetCollection
```

## Documentation

1. [Quickstart](https://phospy.com/docs/quickstart/)
2. [API Guide](https://phospy.com/docs/api/)
3. [Workflow Contracts](https://phospy.com/docs/workflow_contracts/)
4. [Validation Guide](https://phospy.com/docs/validation/)
5. [Scientific Coverage Matrix](https://phospy.com/docs/scientific-coverage/)

## Citation

If you use PhosPy in scientific work, cite this software release using
[`CITATION.cff`](CITATION.cff) and also cite the upstream PhosR project and
publications described in [`NOTICE.md`](NOTICE.md).
