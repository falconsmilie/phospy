# PhosPy

[![PyPI version](https://img.shields.io/pypi/v/phospy.svg)](https://pypi.org/project/phospy/)
[![Python versions](https://img.shields.io/pypi/pyversions/phospy.svg)](https://pypi.org/project/phospy/)
[![Tests](https://github.com/falconsmilie/phospy/actions/workflows/ci.yml/badge.svg)](https://github.com/falconsmilie/phospy/actions/workflows/ci.yml)
[![License](https://img.shields.io/pypi/l/phospy.svg)](https://github.com/falconsmilie/phospy/blob/main/LICENSE)
[![Downloads](https://api.pepy.tech/badge/phospy)](https://pepy.tech/project/phospy)

PhosPy is a Python package for selected phosphoproteomics workflows inspired by
PhosR. It is aimed at scientists who want a clear Python lane from phosphosite
intensity tables to differential phosphorylation analysis, offline
over-representation enrichment, kinase scoring and prediction, and optional
signalome analysis. Enrichment support is offline ORA over caller-supplied
selected identifiers, local set collections, and an explicit background
universe.

"PhosR-inspired" in PhosPy docs means scoped, feature-level comparison lanes. It
does not imply full PhosR package parity or full PhosR API compatibility.
Current differential analysis is scoped to tested design and contrast envelopes;
it is not full limma or PhosR parity.

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
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

For public release checks, the maintainer command is `make release-check`. It
runs the normal lint, type, unit, blocking parity, performance,
release/golden/reproducibility, checked-in reference, and distribution build
checks:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
make release-check
```

This process provides normal CI/build confidence, not formal
exact-source/exact-artifact attestation. Partial local passes, such as default
tests, parity-only tests, or performance-only tests, are useful development
checks but are insufficient for publishing. Default pytest `testpaths` omit
`tests/release`, `tests/golden`, and `tests/performance`; the
`test-release-gates` Make target selects release/golden checks explicitly.
`make build` starts from an empty `dist/`, builds one wheel and one sdist, runs
metadata checks, and validates the packaged reference manifests and declared
file hashes in both archives.

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
Packaged references are governed by a manifest `redistribution_status`: only
`approved` references with typed upstream-package license evidence for the exact
snapshot and packaged files are release-eligible, `unresolved` bundled
references block release, and `external_only` references are caller-supplied
local data that must not be shipped as bundled data. The raw
`redistribution_allowed` manifest value is only a compatibility mirror and must
not contradict `redistribution_status`.

The default kinase `scoring_mode="phosr_rank_weighted"` is PhosR-inspired
rank-weighted scoring implemented by PhosPy. It combines available profile and
motif support under PhosPy's support rules; it is not an exact PhosR
implementation and is not intended to provide numerical parity with PhosR.
Kinase scoring and prediction outputs are relative support within a run, not
calibrated probabilities.
Optional kinase activity outputs are exploratory kinase activity scores or
activity-like substrate summaries. They depend on substrate coverage and
reference evidence; they are not direct proof of kinase activation or causal
pathway activity.

Enrichment ORA results are overlap statistics under the caller-supplied
background universe. They do not prove pathway activation, regulation, or
biological causality, and PhosPy does not imply GSEA or PTM-SEA support.

Scientific scope categories and parity/open-gap status are maintained in
[`docs/scientific-coverage.md`](docs/scientific-coverage.md). Parity fixture
evidence lives in [`docs/parity.md`](docs/parity.md). Parity claims are
fixture-scoped; they do not transfer to untested fixtures, broader PhosR
surfaces, or artifacts that did not pass the maintainer release checks.
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
condition metadata, a missingness policy, diagnostics, and provenance. The
native PhosPy SPS/RUV-style preprocessing correction estimates unwanted factors
from eligible control-site residuals after protected-design handling. Batch
terms are resolved for validation and diagnostics, including
batch-associated-variance summaries; they are not directly residualized as
fixed effects by the native correction. Replicate metadata, when supplied, is
validated and recorded for provenance and diagnostics only. Supplied replicate
labels must not be all the same, all unique, perfectly confounded with batch,
or perfectly confounded with protected condition metadata. Replicate metadata
is not used for numerical unwanted-factor estimation and does not enable
RUV-III or replicate-aware RUV-III semantics. RUV-III style correction is not
executable unless a future feature implements replicate-aware semantics. Any
`ruv_readiness` diagnostics are report-only readiness signals and do not apply
correction. Native SPS/RUV-style correction rejects correction-stage matrices
with actual missing values before executor invocation: temporary imputation
followed by restored missing values cannot produce analysis-ready corrected
output. Run missing-data preprocessing first, or provide a complete
upstream-imputed matrix with an observation mask.

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
from phospy.api.configs import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
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
        scoring_config=KinaseScoringConfig(
            reliability_profile=KinaseReliabilityProfile.CUSTOM,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            )
        ),
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
from phospy import AnalysisReadyDatasetBuilder
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
```

`AnalysisReadyPhosphoDataset` remains exported as a public result/domain type.
The direct `AnalysisReadyPhosphoDataset(...)` constructor raises immediately; it
does not provide a warning-based compatibility route and does not create
provenance for direct calls. Ordinary user construction should go
through `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`.
Advanced callers who already own fully prepared `site_key`-indexed
analysis-ready tables must use
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with complete
`TrustedDatasetConstructionAssertions`. The trusted factory runs private dataset
structural validation, including required `site_sequence`, but trusted
assertions and provenance are caller claims and do not prove biological
correctness. Supplied provenance must fingerprint the actual represented tables;
mismatched fingerprints are rejected.

Use `phospy.api` for the stable request, workflow, primary result, reference,
enum, and common exception names documented in the API guide. The aggregate
facade is intentionally smaller than the implementation modules; validators,
workflow executors, processing-state internals, nested diagnostic records, and
compatibility constants are not stable public API.

`EnrichmentWorkflow` is a supported public workflow from `phospy.api`, not a
top-level `phospy` convenience export:

```python
from phospy.api import EnrichmentConfig, EnrichmentWorkflow
from phospy.api import EnrichmentWorkflowRequest, GeneSetCollection
```

Specialized configuration and inspection helpers that are not part of the
stable facade are documented as advanced supported API. Import lower-level
constants and nested result models from explicit submodules such as
`phospy.api.configs` or `phospy.api.results` only when you need that advanced
surface.

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
