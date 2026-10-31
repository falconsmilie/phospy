# PhosPy Release Notes

## Version 1.8.1

Publication date: 2026-10-31.

Compare the complete change set from
[v1.8.0 to v1.8.1](https://github.com/falconsmilie/phospy/compare/v1.8.0...v1.8.1),
or read the [Version 1.8.0 release notes](release-notes-1.8.0.md).

## Release Overview

PhosPy 1.8.1 completes the supported API and owner-namespace consolidation
started in Version 1.8.0. Ordinary scientific imports now come from
`phospy.api`, specialist imports come from `phospy.api.advanced`, and the
package root exposes version metadata only. The release also adds a governed
route for supplying caller-attested output from an external batch correction
and makes kinase and signalome bundles preserve the complete kinase-scoring
state.

This release does not add or change a scientific estimator, alter analysis
defaults, expand an external parity claim, or strengthen biological
interpretation. Python support remains 3.11 and 3.12, with
`Requires-Python: >=3.11,<3.13`.

## Compatibility Warning

Although numbered as a patch release, Version 1.8.1 contains breaking import
and reloadable-bundle changes. This follows the repository's accepted beta API
reset policy, but applications that import from the Version 1.8.0 package root
or `phospy.advanced` must update before upgrading. No forwarding modules or
compatibility shims are provided.

Install and test the candidate in a clean environment before updating a
production analysis environment:

```bash
python -m pip install --upgrade "phospy==1.8.1"
```

## Update Supported Imports

Move the three former package-root conveniences to the ordinary aggregate:

```python
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

Move specialist configuration, diagnostics, reference helpers, and publishers
from `phospy.advanced` to the nested specialist aggregate:

```python
from phospy.api.advanced import (
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
    publish_kinase_workflow,
)
```

The ordinary `phospy.api` inventory retains its 66 supported names. The
specialist inventory retains its Version 1.8.0 names and adds
`CorrectedPreprocessingOutput` and
`create_external_corrected_preprocessing_output`, for 126 names in total.

Reusable implementation packages such as the former `phospy.errors`,
`phospy.references`, `phospy.quantitative`, `phospy.provenance`, and
`phospy.validation` locations moved beneath `phospy.core`. Those owner modules
remain unsupported external import targets. Import a documented retained name
from `phospy.api` or `phospy.api.advanced`; internal-only helpers have no
supported replacement.

## Supply Externally Corrected Preprocessing Output

Version 1.8.1 adds an advanced typed route for a supported SPS/RUV-style
correction executed outside PhosPy:

```python
from phospy.api import AnalysisReadyDatasetBuilder, DatasetBuildRequest
from phospy.api.advanced import create_external_corrected_preprocessing_output

corrected_output = create_external_corrected_preprocessing_output(
    pre_correction_matrix=phospho,
    corrected_matrix=externally_corrected,
    method="sps_ruv_style",
    source_software="external-correction-tool",
    source_version="2.1",
    control_site_set=controls,
    n_unwanted_factors=1,
    batch_column="batch",
    condition_columns=("condition",),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=organism,
        input_intensity_scale="log2",
        preprocessing_config=preprocessing_config,
        corrected_preprocessing_output=corrected_output,
    )
)
```

The factory supports `sps_ruv_style` and `ruv_iii_style`; the latter also
requires `replicate_column`. It requires finite, axis-aligned input and output,
eligible controls, a positive unwanted-factor count, a batch column, protected
condition columns, and preservation of condition effects. PhosPy validates and
fingerprints the supplied state, but records the execution as caller-attested
and does not claim to have verified the external method or environment.

Externally corrected output must be the only matrix-changing preprocessing
input. It cannot be combined with configured total-protein correction,
site-matrix construction, normalisation, comparison building, or another batch
correction. See [Prepare a Dataset](api/dataset-build-workflow.md) for the full
validated example and rejection boundaries.

## Regenerate Reloadable Bundles

Current internal bundle writers emit kinase manifest version 4 and signalome
manifest version 3. Both formats now store the complete typed kinase-scoring
state, including all optional score tables, score-source data, diagnostics,
validation results, policies, dtypes, and `None` versus empty distinctions.
Loaders reconstruct the original `KinaseScoringResult` through its validating
constructor instead of substituting defaults.

Readers accept only the current schema. Kinase bundle versions 1 through 3 and
signalome versions 1 and 2 are rejected because they do not contain enough
state for faithful reconstruction. Regenerate reloadable bundles with Version
1.8.1 before relying on them. User-facing `publish_*` table directories remain
the supported output route and are distinct from these internal reloadable
codecs. See [Output Bundles](output_bundles.md) for details.

## Runtime Version and Release Verification

`pyproject.toml` remains the single package-version authority. An installed
distribution now takes precedence when resolving `phospy.__version__`; an
uninstalled fixed-layout source checkout falls back only to its own verified
project file and otherwise reports `0+unknown`. Reinstall editable metadata
after changing the project version.

Release tooling now:

- verifies architecture inventory generation through its package module entry;
- checks all supported public aggregates during installed wheel and source
  distribution verification;
- enforces namespace-aware logical ownership after the `phospy.core` moves; and
- keeps strict Pyright coverage aligned with every supported aggregate and
  owner package.

These are release and architecture checks, not new scientific capabilities or
performance claims.

## Scientific and Environment Continuity

The directory and API changes do not change the supported workflow
request-and-`run()` shape, scientific calculations, defaults, dataset
invariants, dependency lower bounds, or external validation scope. Existing
PhosR, limma, RUV-III, and other parity claims remain limited to their
documented fixtures and semantics. The external-correction handoff records
caller declarations and validation evidence; it does not create a new native
estimator or a broader PhosR/RUV equivalence claim.

Review [Scientific Coverage](scientific-coverage.md), [Parity](parity.md), and
the workflow-specific caveats before interpreting results.

## Continue With the Current Guides

After updating imports and regenerating any reloadable bundles, run the
[Quickstart](quickstart.md) or the relevant workflow guide. The
[Public Python API](api/guide.md) is the supported import authority, and the
[changelog](https://github.com/falconsmilie/phospy/blob/main/CHANGELOG.md)
preserves the complete release history.
