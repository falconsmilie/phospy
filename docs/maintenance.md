# Maintenance

This page is the authoritative maintainer release-process page for PhosPy.

## Development Setup

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # optional parquet support
```

For CI-aligned dependency resolution:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,docs]"
```

For minimum supported dependency validation, use the dedicated lower-bound
constraint file rather than the current pinned CI stack:

```bash
pip install -c constraints/minimum.txt -e ".[test]"
python -m pip check
pytest -m "not parity"
pytest -o addopts= tests/release tests/golden -m "release_gate or golden or reproducibility"
```

For release checks, install the release extras first:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet,docs]"
```

Then run the final aggregate maintainer command:

```bash
make release-check
```

If `make release-check` fails with import errors for optional engines, install
the optional extras above and rerun.

## Common Checks

```bash
ruff check .
ruff format --check .
python tools/testing/pyright_strict_coverage.py --check
python scripts/run_pyright.py
pytest -m "not parity"
```

Run blocking parity tests when scientific logic, fixture data, reference
handling, or scoring behaviour changes:

```bash
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

Run performance checks when preprocessing, scoring, prediction, or signalome hot
paths change:

```bash
pytest tests/performance -m "performance or release_gate"
```

Run the optional local release-scale benchmark only when you need same-machine
scale profiling:

```bash
make benchmark-release-scale
```

This 50,000 x 48 builder+differential benchmark is informational,
machine-dependent, and excluded from `make test-performance`,
`make release-check`, and GitHub Actions.

Run release/golden checks when changing release policy, provenance fixtures,
reference manifests, or reproducibility contracts:

```bash
make test-release-gates
```

Run the maintainer release checks before tagging a release:

```bash
make release-check
```

Default `pytest` or `pytest -m "not parity"` is a fast local development check,
not sufficient for publishing. Ordinary CI/build success provides normal
development confidence: it shows the source-tree checks and packaging checks
selected by CI are healthy for that commit. Final release verification is
different. `make release-check` is the authoritative aggregate command, and it
performs release-owned exact-source and exact-artefact verification by
validating staged reference-bundle bytes from a Git-backed checkout, building
fresh wheel and sdist artefacts, checking packaged reference manifests and
hashes, and executing installed wheel/sdist probes outside the checkout.
Because staged-byte verification reads the Git index, release checks must run
from a Git-backed checkout. A successful source-tree test run alone is not
proof that the built wheel and sdist are valid.

The configured default pytest `testpaths` omit `tests/release`, `tests/golden`,
and `tests/performance`. Blocking parity tests, performance contracts,
release/golden/reproducibility tests, reference validation, metadata checks,
archive-level packaged-reference checks, and the installed wheel/sdist verifier
are not optional for public release decisions.

The publish pipeline (`.github/workflows/publish.yml`) runs
`make release-check` once on the checked-out tag, uploads the fresh `dist/`
directory, verifies the uploaded wheel and sdist on Python 3.11 and 3.12, and
publishes those artifacts through trusted publishing only after that matrix
passes.

Do not treat a partial local pass, a parity-only pass, a performance-only pass,
ordinary CI/build confidence from a different commit, or a check pass from a
different distribution as sufficient evidence for public release.

Release-blocking coverage in `make release-check` is:

| Gate | Command selector |
| --- | --- |
| Lint | `ruff check .` |
| Type checking | `python scripts/run_pyright.py` |
| Default non-parity suite | `pytest -m "not parity"` |
| External-consumer public API contract | `pytest -o addopts= tests/contract` |
| Threshold-bearing parity | `pytest tests/parity -m "parity and not parity_diagnostic" -s` |
| Performance release contracts | `pytest tests/performance -m "performance or release_gate"` |
| Checked-in reference bundles | `python scripts/validate_reference_bundle_index.py --repo-root .` |
| Release/golden/reproducibility gates | `pytest -o addopts= tests/release tests/golden -m "release_gate or golden or reproducibility"` |
| Distribution build and packaged-reference checks | `make build` |
| Installed wheel/sdist verification | `make verify-installed-distributions` |

The optional `docs-build` and `benchmark-release-scale` targets are deliberately
absent from this table. `docs-build` remains the explicit local documentation
validation command. Benchmark runtime and memory observations are local
benchmark data, not release-blocking evidence.

`parity_diagnostic` checks are intentionally excluded from the blocking parity
target unless a maintainer deliberately promotes them into the release selector
and updates this policy. `make release-check` is the authoritative aggregate
command. The release policy tests use
`tools/testing/release_selector_coverage.py` to collect pytest node IDs and
effective markers, then fail if any release-blocking node is absent from every
authoritative release target.

## Type Checking

Pyright is the configured type checker. The checked scope is listed in
`pyproject.toml` under `[tool.pyright]` and includes:

- `src/phospy/api`
- `src/phospy/contracts`
- `src/phospy/errors`
- `src/phospy/frames`
- `src/phospy/io`
- `src/phospy/policies`
- `src/phospy/provenance`
- `src/phospy/science`
- `src/phospy/tables`
- `src/phospy/validation`
- `src/phospy/workflows`

Strict checking is enabled for selected stable scientific/core modules listed
under `[tool.pyright].strict`. This strict list currently includes
`src/phospy/science/datasets/models.py` and
`src/phospy/contracts/results`; those paths are already strict-checked, not
future targets. Strict scope can be expanded incrementally. The strict coverage
policy check fails if any declared strict path is missing or outside the
configured Pyright include scope.

Avoid suppressions by default. Use them only when Pyright cannot model correct
runtime behaviour. In strict files, every Pyright suppression must use this
inline format:

```python
# pyright: ignore[reportRuleName] - concrete technical rationale
```

The rule list must contain one or more explicit Pyright `report...` diagnostic
names. The rationale must explain the concrete typing limitation and why the
runtime operation is safe. Blanket `# pyright: ignore`, placeholder reasons such
as `TODO`, file-wide diagnostic downgrades such as
`# pyright: reportUnknownMemberType=false`, file-wide strictness downgrades, and
`[tool.pyright].ignore` entries that intersect effective strict files are
rejected by `python tools/testing/pyright_strict_coverage.py --check`.

## Fixture Policy

Active fixture roots:

- `tests/fixtures/`
- `tests/support/`
- `scripts/active/`

Regeneration scripts should be deterministic and should say which fixture family
they update. Generated benchmark reports belong in `benchmarks/reports/`, which
is ignored by git.

Release-validation fixtures added under `tests/fixtures/` must keep their
classification explicit:

- external parity fixtures must record the external implementation and version,
  generation command, seed, timestamp/source policy, generator SHA-256, and
  file SHA-256 hashes.
- PhosPy regression/golden/performance fixtures must not be described as
  external parity unless an external expected-output source is actually used and
  documented.
- The large differential trend parity fixture is regenerated with
  `make fixtures-large-differential-limma-trend`.
- Compact PhosPy release-validation regression fixtures are regenerated with
  `make fixtures-release-validation-regression`.

Manifest-governed text fixtures are exact-byte assets. Their project-standard byte
policy is UTF-8 with LF line endings and a final newline for CSV, JSON,
Markdown, and `MANIFEST.json` files. Generators must write bytes explicitly
rather than relying on platform text-mode newline translation, and manifests
must hash the final bytes written. Tests must validate raw checked-in bytes; do
not normalize line endings before hashing. The repository pins
`tests/fixtures/release_validation_regression/**`,
`tests/fixtures/rewrite_parity/differential_limma_trend_large/**`, and the
active generators that record source hashes to LF via `.gitattributes`.

After changing a manifest-governed fixture generator, run the relevant
regeneration target and the release fixture reproducibility gate:

```bash
make fixtures-release-validation-regression
make fixtures-large-differential-limma-trend
pytest -o addopts= tests/release/test_manifest_fixture_byte_reproducibility.py
```

## Source and Release Archive Hygiene

Source and release archives should be built from the tagged source state. Use a
clean tree and `make build` for package distributions after
`make release-check` passes:

```bash
make build
```

`make build` clears stale wheel/sdist files, builds one wheel and one sdist,
runs metadata checks, and validates both built archives against their packaged
reference manifests and declared file hashes. It does not require Git metadata
and can run from a copied source tree. This is an archive-level check; it is
complementary to installed execution.

Run the installed-distribution verifier through `make release-check`, or invoke
the verifier target directly to rebuild the artifacts and execute the installed
checks:

```bash
make verify-installed-distributions
```

The verifier creates separate temporary environments outside the checkout,
installs the wheel and sdist, runs Python with isolation enabled, checks that
`phospy.__file__` resolves inside the installed environment, verifies bundled
rat reference manifest files and SHA-256 values, and runs representative
dataset, differential, and kinase public workflow contracts without importing
repository tests, fixtures, or `conftest.py`.

Do not broaden a release claim from one artifact to another. If the tagged
source tree, source archive, wheel, dependency constraints, or bundled fixture
set changes after a release-check run, rerun the checks and rebuild `dist/`.

Release confidence is based on the ordinary CI matrix, the maintainer
release-check command, and the freshly built wheel/sdist artifacts.

Do not include generated artefacts in source/release archives. Exclude build
outputs, documentation sites, cache directories, previous archive files, and
local benchmark reports, including `build/`, `dist/`, `site/`, `.pytest_cache/`,
`.ruff_cache/`, `.hypothesis/`, `__pycache__/`, `benchmarks/reports/`, and
files such as `app-src.zip`.

Before publishing or handing off an archive, verify that it contains source,
tests, docs, configuration, constraints, examples, and required metadata, but
not generated output directories or prior release artefacts.

## Documentation Policy

Docs should stay flat, beginner-friendly, and tested against the code. Prefer one
clear beginner path over several overlapping overview pages. Keep examples small
and runnable.

Docs subdirectories are intentional:

- `docs/api/` for workflow API references
- `docs/adr/` for architecture decision records
- `docs/testing/` for testing-audit and consolidation material

## Config Ownership

Config concepts must have one implementation owner. Public transport DTOs that
describe request shape or optional user intent belong under
`phospy.contracts.configs` and stay available through `phospy.api.configs`.
Science-owned algorithm policies belong in the owning science domain, and
contract/API routes may re-export them only by exact object identity.

Workflow interpreters own contextual translation from public config DTOs into
distinct resolved execution models. Numerical science modules must not import
`phospy.contracts` or accept public workflow config DTOs. The architecture guard
for this is `tests/architecture/test_config_ownership.py`.

### Scientific Claims Checklist

Before merging docs or examples that describe scientific results, check that the
wording:

- describes kinase outputs as scores or associations, not causal kinase activity
- avoids claiming exact PhosR equivalence unless scoped parity evidence supports it
- treats enrichment as statistical evidence, not pathway activation or biological
  proof
- names the background universe, test, correction, threshold, and assumptions
  before calling a result significant
- keeps `site_key` as row identity and `display_id` as a display label
- describes score-derived network edges as inferred, not experimentally proven
  relationships
- respects the `AnalysisReadyPhosphoDataset` boundary, private dataset
  validation, and validator -> interpreter -> executor responsibilities

## Deprecation and Removal Ledger

This ledger is the maintainer checkpoint for deprecated APIs, retained
backwards-compatibility paths, and removals that tests intentionally keep absent.
Add a row here before adding a new compatibility path or removing one.

### Warning-Based Deprecations

| Surface | Module | Status | Replacement or recommended route | Warning | Decision and removal target | Intentional tests |
| --- | --- | --- | --- | --- | --- | --- |
| `AnalysisReadyPhosphoDataset(...)` | `src/phospy/science/datasets/models.py` | Sealed public constructor; the class remains importable as a result/domain type. | `AnalysisReadyDatasetBuilder.run(DatasetBuildRequest(...))` for ordinary construction, or `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with seven required `TrustedDatasetConstructionAssertions` dimensions plus an optional typed numeric-semantic-domain waiver when needed for trusted caller-owned tables. | Raises `TypeError` immediately. | No warning-based compatibility path remains. | `tests/unit/science/datasets/test_direct_construction_assertions.py::test_direct_constructor_fails_immediately_and_names_supported_paths`. |
| `DifferentialAnalysis` | `src/phospy/science/differential/public.py` | Internal compatibility path; not exported from `phospy` or `phospy.api`. | `DifferentialAnalysisWorkflow.run(DifferentialAnalysisRequest)` from `phospy` or `phospy.api`. | `PhosPyDeprecationWarning` on construction. | Retain for now; remove in a future release. | `tests/unit/test_public_contract_import_routes.py::test_deprecated_differential_analysis_shell_warns_and_delegates`, `tests/unit/test_public_contract_import_routes.py::test_differential_analysis_is_not_supported_from_phospy_api_namespace`, `tests/unit/test_public_documentation_examples.py::test_readme_differential_import_example_matches_supported_route`. |
| `TechnicalReplicateResolver` | `src/phospy/workflows/differential/replicates.py` | Internal workflow compatibility wrapper. | `TechnicalReplicateAggregationPlanner` plus `TechnicalReplicateAggregator`, or normal `DifferentialAnalysisWorkflow.run(...)`. | `PhosPyDeprecationWarning` on construction and `run()`. | Retain for now; remove in a future release. | `tests/unit/test_differential_technical_replicate_policy.py::test_technical_replicate_resolver_warns_and_preserves_wrapper_behaviour`. |
| `load_enrichment_sets_gmt`, `load_enrichment_sets_table`, `load_enrichment_sets_csv`, `load_enrichment_sets_tsv` | `src/phospy/io/readers/enrichment_sets.py` | Public reader aliases. | Matching `read_enrichment_sets_*` function. | `PhosPyDeprecationWarning` when called. | Retain for now; remove in a future release. | `tests/unit/test_reader_enrichment_sets.py::test_reader_enrichment_load_aliases_warn_and_forward`. |
| Post-hoc peptide-to-site differential estimate combination | `src/phospy/science/differential/aggregation/experimental.py`, `src/phospy/science/differential/aggregation/models.py`, `src/phospy/science/differential/aggregation/executor.py` | Withdrawn from public support under `unsupported_withdrawn_posthoc_estimate_combination_v1`. | Resolve peptide evidence at sample-intensity level before `DifferentialAnalysisWorkflow`. Future public support requires executable mapping semantics, a coherent combined estimand, an inferential result, dependence handling, multiple-testing semantics, provenance semantics, docs, and tests. | Yes. The retained compatibility shell fails closed before calculation because coherent combined effect/inference and mapping semantics are not implemented. Mapping policies such as equal splitting or statistical-model exclusion must not silently execute as ordinary evidence. | Retain science-owned source only as internal future-work material; do not export it as production functionality or delegate public calls to calculation. | `tests/unit/test_peptide_to_site_aggregation.py`. |
| `KinaseActivityResult.activity_scores`, `KinaseActivityResult.weighted_activity` | `src/phospy/science/activities/results.py`; `src/phospy/science/activities/models.py` is an identity-preserving import route | Public result aliases. | `KinaseActivityResult.activity_matrix`. | `PhosPyDeprecationWarning` on property access. | Retain for now; remove in a future release after bundle/result consumers have migrated. | `tests/unit/test_activity_science.py::test_activity_scores_compatibility_alias_matches_activity_matrix`, `tests/unit/test_activity_science.py::test_weighted_activity_compatibility_alias_matches_activity_matrix`, `tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_kinase_bundle_integration.py`. |
| `KinaseActivityResult.legacy_condition_statistics_table_dataframe()` | `src/phospy/science/activities/results.py`; `src/phospy/science/activities/models.py` is an identity-preserving import route | Compatibility adapter for old condition-shaped activity statistics tables. | `KinaseActivityResult.statistics_table_dataframe()` or `KinaseActivityResult.statistics_table` with `profile_id` as the row identity. | `PhosPyDeprecationWarning` on method call. | Retain for now; remove after downstream statistics-table consumers have migrated to `profile_id`. | `tests/unit/test_activity_science.py::test_legacy_condition_statistics_table_adapter_is_deprecated_and_defensive`, `tests/integration/test_kinase_bundle_integration.py::test_kinase_bundle_round_trip_preserves_exact_activity_semantics`. |
| `KinaseScoringConfig.default()` | `src/phospy/contracts/configs/kinase.py` | Public compatibility constructor alias for the historical exploratory kinase scoring defaults. | `KinaseScoringConfig.exploratory()` for the old permissive behavior, or `KinaseScoringConfig.production(...)` with explicit study-specific attrition floors. | `PhosPyDeprecationWarning` when called. | Retain for now; remove in a future release after public examples and downstream callers have migrated. | `tests/unit/test_public_config_self_validation.py::test_kinase_scoring_default_is_deprecated_exploratory_alias`. |
| Bare string entries in `build_motif_library_from_sequences(motif_sequences=...)` | `src/phospy/science/prediction/motif_scoring/library_validation.py` | Semi-public motif-library input compatibility path. | `ExplicitMotifSequence` values or mapping entries with `reference_id`, optional `site_id`, `kinase`, and `sequence`. | `PhosPyDeprecationWarning` when any bare string entry is present. | Retain for now; reject in a future release. | `tests/unit/test_prediction_sequence_validation.py::test_structured_explicit_sequences_do_not_emit_bare_string_warning`, `tests/unit/test_prediction_sequence_validation.py::test_bare_explicit_sequences_warn_and_remain_supported_without_site_mismatch_claims`, `tests/unit/test_prediction_sequence_validation.py::test_bare_explicit_invalid_sequence_is_rejected_and_reported`. |

### Retained Compatibility Paths Without Deprecation Warnings

| Surface | Module | Status | Replacement or recommended route | Warning | Decision and removal target | Intentional tests |
| --- | --- | --- | --- | --- | --- | --- |
| `phospy.api.configs.*`, `phospy.api.requests`, `phospy.api.results` wrapper modules | `src/phospy/api/` | Public API wrapper routes. | Use `phospy.api` for public request, config, result, enum, reference, workflow, and exception contracts. | No. | Retain; no removal target. | `tests/unit/test_public_contract_workflows.py`, `tests/unit/test_public_contract_results.py`, `tests/unit/test_public_contract_import_routes.py`. |
| Dataset request constants for supported analysis-ready policies (`reject`, `exclude_from_sequence_scoring`, `split`) and `DATASET_SITE_RESOLUTION_MODE_*` | `src/phospy/api/requests.py`, `src/phospy/contracts/requests.py`, `src/phospy/science/evidence/dataset_resolution/` | Public request compatibility constants for policies that can produce strict site-level rows. | Use the documented `DatasetBuildRequest` policy fields; `keep_joint` is no longer an analysis-ready builder policy and fails request validation with migration guidance. | No. | Retain the supported constants only. | `tests/unit/test_public_contract_workflows.py::test_request_compatibility_constants_are_public_exports`, `tests/unit/test_dataset_peptide_evidence_resolution.py`. |
| `KinaseEligibilityReport`, `KinaseWorkflowPreprocessingAttritionSummary`, `KinaseWorkflowScoringAttritionSummary`, `KinaseWorkflowSiteAttritionSummary` | `src/phospy/api/results.py`, `src/phospy/contracts/results.py` | Public result compatibility aliases. | Prefer accessing these through `KinaseWorkflowResult` fields unless type names are needed. | No. | Retain; no removal target. | `tests/unit/test_public_contract_results.py::test_result_compatibility_aliases_are_public_exports`, `tests/unit/test_kinase_workflow_components.py`. |
| Semi-public science routes from ADR-0028: `PreprocessingStageMetadata`, clustering protocol/facade names, and `fuse_profile_and_motif_scores_by_rank_weight` | `src/phospy/science/datasets/preprocessing/stage_registry.py`, `src/phospy/science/signalomes/clustering/`, `src/phospy/science/prediction/scoring.py` | Semi-public compatibility routes for extension, parity, performance, and backend-contract use. | Use only the route-specific names listed in ADR-0028; do not promote them into `phospy` or `phospy.api`. | No. | Retain until a normal deprecation/removal process replaces ADR-0028. | `tests/unit/test_public_contract_import_routes.py`, `tests/parity/test_signalome_clustering_backend_parity.py`, `tests/unit/test_prediction_science.py`. |
| `GeneSetCollection`/`PtmSetCollection` mapping-style `sets` input and `.sets` mapping view | `src/phospy/science/enrichment/models.py` | Public enrichment collection convenience compatibility. | `EnrichmentSetCollection` with explicit `EnrichmentSet` objects; use `members_by_set_id` for a mapping snapshot. | No. | Retain; no removal target. | `tests/unit/test_enrichment_collections.py::test_enrichment_collection_legacy_gene_wrapper_preserves_mapping_view`, `tests/unit/test_result_snapshot_helpers.py::test_enrichment_collection_mappings_are_fresh_snapshots`. |
| `EnrichmentWorkflowResult.result_table` | `src/phospy/contracts/results.py` | Public result-table alias. | `EnrichmentWorkflowResult.table` or `to_dataframe()`. | No. | Retain; no removal target. | `tests/unit/test_result_snapshot_helpers.py::test_enrichment_result_table_aliases_return_equal_independent_snapshots`. |
| Builder display-indexed input plus narrow site-metadata aliases (`gene_name`, `centralized_sequence`, localisation spelling variants) | `src/phospy/science/datasets/builders/site_metadata_normalizer.py` | Public builder-ingestion compatibility only; direct `AnalysisReadyPhosphoDataset` construction remains strict. | Encoded `site_key` identity and explicit metadata columns such as `gene_symbol`, `site_sequence`, and `localisation_confidence`. | No. | Retain at the builder boundary; no removal target. | `tests/unit/test_dataset_builders_site_key_indexing.py`, `tests/unit/test_dataset_normalizer_identity.py`, `tests/unit/test_dataset_normalizer_collaborators.py`, `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_supports_documented_alias_and_index_derivation_conventions`. |
| Reference source-column aliases accepted by `ReferenceBundleBuilder` | `src/phospy/science/references/builder.py` | Public reference-builder input compatibility. | Normalised output columns `kinase`, `substrate_site`, `site_sequence`, `display_id`, organism, and optional protein metadata. | No. | Retain while error messages list accepted aliases; no removal target. | `tests/unit/test_reference_bundle_builder.py`. |
| Mode-only replacement semantics for old site-sequence resolution plans | `src/phospy/science/datasets/preprocessing/site_sequence/reference_loader.py` | Internal preprocessing plan path. | Explicit `site_sequence_resolution.conflict_policy="replace_existing"`. | No. | Retain for old internal plans; no removal target. | `tests/unit/test_dataset_site_sequence_resolution.py`. |
| Missing `row_medians_used`, missing `site_sequence_resolution`, missing report-only `ruv_readiness`, and older site-sequence fields in processing-state payloads | `src/phospy/science/datasets/_processing_state/`, `src/phospy/io/bundles/_shared/processing_state.py` | Internal saved-payload path. | Current `DatasetProcessingState` payload with versioned diagnostics and full site-sequence/readiness fields. | No. | Retain to avoid breaking saved-result readability; no removal target. | `tests/unit/test_processing_state_bundle_payload.py`. |
| Reference provenance payloads without `identifier_normalisation` | `src/phospy/provenance/serialization.py` | Internal saved-provenance path. | Current reference provenance payload including `identifier_normalisation` when available. | No. | Retain; no removal target. | `tests/unit/test_reference_provenance.py::test_reference_provenance_from_payload_supports_legacy_missing_identifier_normalisation`. |
| `MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY` / `"first_site_compatibility"` | `src/phospy/science/evidence/multi_site.py` | Internal explicit opt-in compatibility policy. | Prefer workflow-specific multi-site policies such as keeping joint sites, excluding ambiguous sequence scoring, or split policies. | No. | Retain only as explicit opt-in; no removal target. | `tests/unit/test_evidence_multi_site_handling.py::test_first_site_compatibility_is_only_enabled_explicitly`. |
| Activity sidecars `thresholded_substrate_mean_activity`, `thresholded_substrate_counts`, and `activity_substrate_counts` | `src/phospy/science/activities/results.py`; `src/phospy/science/activities/models.py` is an identity-preserving import route; `src/phospy/io/bundles/_kinase/`, `src/phospy/io/bundles/_signalome/` | Public result and bundle payload compatibility. | Prefer method-neutral `activity_matrix` and `substrate_count_matrix` for primary activity score outputs. | No. | Retain for current activity methods and bundle round trips; no removal target. | `tests/unit/test_activity_science.py`, `tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_kinase_bundle_integration.py`. |

### Already Removed or Explicitly Rejected Paths

| Surface | Module | Status | Replacement or recommended route | Warning | Decision and removal target | Intentional tests |
| --- | --- | --- | --- | --- | --- | --- |
| Legacy provenance schemas and table/stage hash aliases | `src/phospy/provenance/serialization.py` | Removed; payloads fail clearly. | Regenerate provenance with the current schema. | No; raises `PhosPyInputError`. | Already removed in 1.6.0. | `tests/unit/test_provenance_hashing.py`, `tests/unit/test_dataset_run_provenance.py`. |
| Legacy kinase and signalome bundle schemas, including old diagnostic fields | `src/phospy/io/bundles/_kinase/`, `src/phospy/io/bundles/_signalome/` | Removed; bundle load fails clearly. | Regenerate bundles with the current PhosPy version. | No; raises `PhosPyInputError`. | Already removed in 1.6.0. | `tests/integration/test_kinase_bundle_integration.py`, `tests/integration/test_signalome_bundle_integration.py`, `tests/unit/test_signalome_bundle_schema.py`, `tests/unit/test_signalome_clustering_diagnostics_helpers.py`. |
| Removed public/config aliases: `ratio_to_total`, `ensemble_size`, `signalome_cutoff`, `max_exact_clustering_sites`, `kinase_network_policy`, `tree_engine`, and old backend-style signalome names | Dataset, kinase, and signalome config/snapshot modules | Removed; construction or snapshot load rejects them. | Current names such as `subtract_log_total`, `n_iterations`, `scientific.substrate_support_cutoff`, `performance.max_exact_tree_sites`, `output.network_policy`, and `clustering.clustering_engine`. | No; raises `TypeError` or `PhosPyInputError`. | Already removed in 1.6.0. | `tests/unit/test_public_config_self_validation.py`, `tests/unit/test_kinase_bundle_snapshots.py`, `tests/unit/test_signalome_bundle_schema.py`, `tests/unit/test_validator_boundaries.py`, `tests/integration/test_dataset_builder_integration.py`. |
| Unsupported historical site-metadata aliases (`sequence`, `protein`, `gene`, `residue`, `phosphosite`) used as substitutes for strict builder columns | `src/phospy/science/datasets/builders/site_metadata_normalizer.py` | Removed/rejected. | Rename to supported columns such as `site_sequence`, `protein_id`, `gene_symbol`, or `site`. | No; raises `UnsupportedInputFormatError`. | Already removed in 1.6.0. | `tests/unit/test_domain_boundaries.py`, `tests/unit/test_dataset_normalizer_identity.py`, `tests/unit/test_dataset_normalizer_collaborators.py`. |
| Removed kinase/signalome workflow component re-export modules | `src/phospy/workflows/` | Removed import routes. | Import workflow classes from `phospy`, `phospy.api.workflows`, or `phospy.workflows.<workflow>`. | No; import raises `ModuleNotFoundError`. | Already removed. | `tests/unit/test_workflow_public_import_contract.py::test_removed_workflow_compatibility_imports_fail`, `tests/unit/test_workflow_public_import_contract.py::test_docs_and_examples_do_not_reference_removed_workflow_compatibility_paths`. |

Review rules:

- Warning-based deprecations must name a replacement, emit `PhosPyDeprecationWarning`,
  and have a test that intentionally preserves the warning during the window.
- Retained compatibility paths must stay thin and must not own new workflow or
  scientific behaviour.
- Removed paths must have absence/rejection tests before the row is closed out.
- Removal planning must not weaken provenance, bundle reconstruction, or
  saved-result readability.

## Frame Ownership Policy

PhosPy treats DataFrames as owned mutable state internally.

Input DataFrames are copied when accepted into validated dataset/table objects.
Workflow internals may pass owned DataFrames without repeated defensive copies.
Public result/table access returns caller-writable defensive snapshots; internal
borrow helpers are read-only by contract and must not be exposed as public API.

Frame ownership helpers must not set or restore pandas process-global options.
Internal borrowed views are mutation-isolated locally: NumPy-backed pandas
frames use shallow read-only borrowed blocks where possible, pandas runtimes
with native copy-on-write may detach locally, and unsupported pandas internals
fall back to deep copies. Code that needs to mutate scientific state must use an
owned frame, not a borrowed frame.

Provenance fingerprints describe the owned internal state at creation time.

Exposure categories:

- `owned_internal`: DataFrames stored in dataset/result/table dataclass fields.
- `safe_public_copy`: `to_dataframe(...)`, `to_pandas(...)`, and
  `*_dataframe(...)` helpers (always caller-writable defensive snapshots).
- `borrowed_internal_view`: private/internal helpers only (`_borrow_dataframe`,
  `_borrow_optional_dataframe`); read-only by contract, and writes may raise or
  detach locally depending on pandas, but must not mutate the owner.
- `export_snapshot`: persisted outputs and provenance fingerprints.

## Release Notes

- Current release notes: [PhosPy Release Notes](release-notes.md)
- Changelog: [`CHANGELOG.md`](https://github.com/falconsmilie/phospy/blob/main/CHANGELOG.md)
- Citation metadata: [`CITATION.cff`](https://github.com/falconsmilie/phospy/blob/main/CITATION.cff)

## ADRs

Architecture and governance decisions live in [ADR Index](adr/index.md). ADRs are
advanced maintainer documents; day-to-day users should start with the
[Quickstart](quickstart.md).
