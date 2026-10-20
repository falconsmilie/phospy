# PhosPy Release Notes

## Version 1.8.0

Provisional release date: 2026-10-20.

Reconfirm this date on the publication day if the release is published later.
These notes prepare Version 1.8.0 for release; they do not indicate that the
package or tag has already been published. They describe changes since Version
1.7.6.

## Release Overview

PhosPy 1.8.0 completes a domain-first package reorganisation and makes a clean
public API break. Supported external imports now come from the aggregate
`phospy`, `phospy.api`, and `phospy.advanced` modules. Existing code that uses
facade submodules, workflow packages, implementation modules, or deprecated
aliases must be updated before upgrading.

The release does not add a scientific method, change analysis defaults, broaden
an external parity claim, or strengthen biological interpretation. Supported
algorithms, dataset invariants, and versioned table and bundle formats continue
across the reorganisation.

## Prepare for the Upgrade

Audit imports and deprecated calls while still using Version 1.7.6. After
Version 1.8.0 is published, install the release in a clean environment and run
representative workflows and bundle round trips before updating production
environments:

```bash
python -m pip install --upgrade "phospy==1.8.0"
```

PhosPy continues to support Python 3.11 and 3.12. The declared runtime
dependency lower bounds are unchanged.

## Use the Aggregate Public Facades

The package root remains a deliberately small convenience surface:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

Use `phospy.api` for ordinary requests, configurations, results, errors,
references, and workflows:

```python
from phospy.api import (
    AnalysisReadyPhosphoDataset,
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    KinaseWorkflowRequest,
    SampleDesignRecord,
)
```

Use `phospy.advanced` only for documented specialist configuration,
diagnostics, resources, and publishing helpers:

```python
from phospy.advanced import (
    DatasetIntensityTransformConfig,
    KinaseScoringConfig,
    SignalomeConfig,
    publish_dataset,
)
```

The [Public Python API](api/guide.md) is the current import authority. A class
or function available from an owner module is still unsupported for external
callers unless it is exported by one of these aggregate facades.

## Migrate Imports

The following table covers common Version 1.7.6 imports. A replacement applies
only when the named object remains present in the destination aggregate.

| Version 1.7.6 code | Version 1.8.0 code or action |
| --- | --- |
| `from phospy import AnalysisReadyPhosphoDataset` | `from phospy.api import AnalysisReadyPhosphoDataset` |
| `from phospy import DifferentialAnalysisWorkflow` | `from phospy.api import DifferentialAnalysisWorkflow` |
| `from phospy.api.<submodule> import Name` | `from phospy.api import Name`; if `Name` is absent, it has no supported ordinary route |
| `from phospy.advanced.configs import Name` or `from phospy.advanced.results import Name` | `from phospy.advanced import Name`; if absent, it has no supported advanced route |
| `phospy.workflows.*`, `phospy.tables.*`, or direct `phospy.science.*` imports | Import a documented retained object from `phospy.api` or `phospy.advanced`; internal-only helpers have no public replacement |
| `DesignMatrix` or `ContrastMatrix` | Build an `ExperimentalDesign` and one or more `Contrast` values, then pass them through `DifferentialAnalysisRequest` |
| `DifferentialAnalysis(...)` | `DifferentialAnalysisWorkflow().run(DifferentialAnalysisRequest(...))` from `phospy.api` |
| Former semi-public preprocessing, clustering, or rank-fusion helpers | Use the documented builder or workflow boundary; there is no public low-level extension replacement in this release |

The aggregate modules re-export the owner-defined objects. Changing the import
does not require a new request shape for a retained object.

## Replace Removed Runtime Behaviour

Version 1.8.0 removes deprecated behaviour immediately. Removed calls raise
normal import, attribute, type, or validation errors; they do not emit a
deprecation warning or forward to the new path.

| Removed behaviour | Supported migration |
| --- | --- |
| `KinaseScoringConfig.default()` | Use `KinaseScoringConfig.exploratory()` for the former permissive preset, `KinaseScoringConfig.production(...)` with study-specific retention floors, or explicit custom settings with `reliability_profile="custom"`. |
| `scoring_mode="kinase_library_motif"` | Use `scoring_mode="kinase_library_contextual_motif"`. |
| `KinaseActivityResult.activity_scores` or `.weighted_activity` | Use `KinaseActivityResult.activity_matrix`. |
| `KinaseActivityResult.legacy_condition_statistics_table_dataframe()` | Use `statistics_table_dataframe()` or `statistics_table`; `profile_id` is the canonical row identity. |
| `ActivityMethodSummary.kinase_condition_pairs_*` fields | Use the corresponding `kinase_profile_pairs_*` fields. Historical payloads containing the former field names are rejected. |
| Activity inputs without typed semantics, a raw `effect_matrix`, or activity results without `input_semantics` and `profile_metadata` | Use the public `KinaseWorkflow` boundary. Code constructing internal activity objects directly has no supported low-level compatibility route. |
| `TechnicalReplicateResolver` | Run `DifferentialAnalysisWorkflow` with the public request, design, and technical-replicate policy. The former resolver was an internal workflow wrapper. |
| `load_enrichment_sets_gmt`, `load_enrichment_sets_table`, `load_enrichment_sets_csv`, or `load_enrichment_sets_tsv` | Read and validate the local file in application code, then construct a public `GeneSetCollection`, `PtmSetCollection`, or `EnrichmentSetCollection`. |
| Bare string entries passed to internal motif-library construction helpers | Use the documented `KinaseLibraryResource` loading route from `phospy.advanced`; the former semi-public helper is not supported. |
| `PreprocessingPipeline(stage_metadata_registry=...)` and unknown-attribute delegation from validated request wrappers | Use documented builder and workflow configuration fields. These implementation hooks have no supported external replacement. |

Ordinary dataset creation still uses the builder:

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=organism,
        preprocessing_config=preprocessing_config,
    )
)
```

Direct `AnalysisReadyPhosphoDataset(...)` construction was already sealed in
Version 1.7.x and remains sealed. Use
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` only for the documented
advanced reconstruction case with complete trusted assertions.

## Review Inputs, Outputs, and Persistence

The API reset changes import locations and removes compatibility behaviour. It
does not change the supported request-and-`run()` workflow shape, default
analysis policies, required analysis-ready dataset invariants, or supported
Python versions.

No persistence or schema version was bumped for the package reorganisation.
Continue to use documented table and bundle formats and verify a representative
round trip after upgrading. Python pickle module paths, class `__module__`
values, and imports from implementation packages are not persistence
contracts; applications relying on them must migrate to the supported
table/bundle interfaces described in [Output Bundles](output_bundles.md).

The fixed sample-order profile-correlation reduction can produce different
last-bit floating-point values from Version 1.7.6 on some numerical backends.
The scoring method, ranking and decision contracts, and supported parity limits
remain unchanged. Applications should compare floating outputs with the
documented scientific tolerances rather than requiring backend-identical bytes.

## Portability, Performance, and Release Checks

Profile scoring now calculates its correlation numerator in a fixed sample
order instead of relying on a backend-dependent BLAS reduction. Adaptive
prediction concurrency is more stable, validated SPS fingerprints are reused,
and log2 conversion handles nullable values correctly. These changes target
portable and reliable execution; they do not establish a general workflow
speedup.

The release infrastructure adds generated architecture inventories,
repository-hygiene checks, distribution-to-committed-Git byte comparison, and
stronger clean wheel and source-distribution consumer checks. Absolute
wall-clock contracts have dedicated isolated lanes so they do not compete with
the broader source/runtime suite. The constrained documentation environment
uses `mkdocs-material` 9.6.23; this is a documentation-tool update rather than
a runtime dependency change.

## Scientific Continuity and Limits

The domain-first move changes code ownership and dependency direction. It does
not change supported estimators, dataset invariants, defaults, scientific
policy, or biological claim boundaries. In particular:

- Kinase outputs remain scores or associations rather than causal activity
  claims unless a documented method supplies eligible inferential output.
- Signalome networks remain descriptive, score-derived summaries rather than
  causal, directional, or experimental evidence.
- Enrichment remains offline over-representation analysis over an explicit
  foreground, background, and caller-supplied collection.
- Existing PhosR, limma, RUV-III, and other external parity claims remain
  limited to their documented fixtures and semantics.
- The profile-correlation portability fix does not introduce a new scientific
  method or broader numerical-equivalence claim.

Review [Scientific Coverage](scientific-coverage.md), [Parity](parity.md), and
the workflow-specific caveats before interpreting results.

## Continue With the Current Guides

After migrating imports, run the [Quickstart](quickstart.md) or the relevant
workflow guide. The [Public Python API](api/guide.md) shows the supported
surfaces, and the project [changelog](https://github.com/falconsmilie/phospy/blob/main/CHANGELOG.md)
preserves the complete history for Version 1.7.6 and earlier releases.
