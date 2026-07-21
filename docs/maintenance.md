# Maintenance

This page describes the maintainer material.

## Development Setup

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # optional parquet support
```

For CI-aligned dependency resolution:

```bash
pip install -c constraints/ci.txt -e ".[dev,test]"
```

For full release checks, install the release extras first. The maintainer
release command is `make release-check`; it runs normal lint, type, unit,
parity, performance, checked-in reference, metadata, packaged-reference, and
build checks:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
```

If `make release-check` fails with import errors for optional engines, install
the optional extras above and rerun.

## Common Checks

```bash
ruff check .
ruff format --check .
pyright
pytest -m "not parity"
```

Run parity tests when scientific logic, fixture data, reference handling, or
scoring behaviour changes:

```bash
pytest tests/parity -m parity -s
```

Run performance checks when preprocessing, scoring, prediction, or signalome hot
paths change:

```bash
pytest tests/performance -m performance
```

Run the maintainer release checks before tagging a release:

```bash
make release-check
```

Default `pytest` or `pytest -m "not parity"` is a fast local development check,
not sufficient for publishing. Parity tests, performance contracts, reference
validation, metadata checks, packaged-reference checks, and the wheel smoke test
are not optional for public release decisions.

The publish pipeline (`.github/workflows/publish.yml`) runs `make
release-check` once on the checked-out tag, uploads the fresh `dist/` directory,
and publishes those artifacts through trusted publishing.

This process provides normal CI/build confidence, not formal
exact-source/exact-artifact attestation. Do not treat a partial local pass, a
parity-only pass, a performance-only pass, or a check pass from a different
commit/distribution as sufficient evidence for public release.

Release-blocking coverage in `make release-check` is:

| Gate | Command selector |
| --- | --- |
| Lint | `ruff check .` |
| Type checking | `python scripts/run_pyright.py` |
| Default non-parity suite | `pytest -m "not parity"` |
| Threshold-bearing parity | `pytest tests/parity -m parity -s` |
| Performance release contracts | `pytest tests/performance -m "performance or release_gate"` |
| Checked-in reference bundles | `python scripts/validate_reference_bundle_index.py --repo-root .` |
| Distribution build and packaged-reference checks | `make build` |

`parity_diagnostic` checks are intentionally non-blocking diagnostics unless a
maintainer deliberately promotes them into the release selector and updates this
policy.

## Type Checking

Pyright is the configured type checker. The checked scope is listed in
`pyproject.toml` under `[tool.pyright]` and includes:

- `src/phospy/api`
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
`src/phospy/science/datasets/models.py`; that file is already strict-checked,
not a future target. Strict scope can be expanded incrementally.

Avoid suppressions by default. Use them only when Pyright cannot model correct
runtime behaviour. Every suppression must be narrow, error-code-specific,
commented, and justified by tests where practical.

## Fixture Policy

Active fixture roots:

- `tests/fixtures/`
- `tests/support/`
- `scripts/active/`

Regeneration scripts should be deterministic and should say which fixture family
they update. Generated benchmark reports belong in `benchmarks/reports/`, which
is ignored by git.

## Source and Release Archive Hygiene

Source and release archives should be built from the tagged source state. Use a
clean tree and `make build` for package distributions after `make
release-check` passes:

```bash
make build
```

`make build` clears stale wheel/sdist files, builds one wheel and one sdist,
runs metadata checks, and validates both built archives against their packaged
reference manifests and declared file hashes. It does not require Git metadata
and can run from a copied source tree.

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
| `AnalysisReadyPhosphoDataset(...)` | `src/phospy/science/datasets/models.py` | Public compatibility constructor for trusted advanced/internal analysis-ready tables. | `AnalysisReadyDatasetBuilder.run(DatasetBuildRequest(...))` for ordinary construction, or `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with seven-dimension `TrustedDatasetConstructionAssertions` for trusted caller-owned tables. | `DeprecationWarning` on direct public construction. | Retain for compatibility; restrict new advanced usage to `from_trusted_tables(...)`. | `tests/unit/science/datasets/test_direct_construction_assertions.py::test_direct_constructor_emits_deprecation_warning`. |
| `DifferentialAnalysis` | `src/phospy/science/differential/public.py` | Internal compatibility path; not exported from `phospy` or `phospy.api`. | `DifferentialAnalysisWorkflow.run(DifferentialAnalysisRequest)` from `phospy` or `phospy.api`. | `DeprecationWarning` on construction. | Retain for now; remove in a future release. | `tests/unit/test_public_contract_import_routes.py::test_deprecated_differential_analysis_shell_warns_and_delegates`, `tests/unit/test_public_contract_import_routes.py::test_differential_analysis_is_not_supported_from_phospy_api_namespace`, `tests/unit/test_public_documentation_examples.py::test_readme_differential_import_example_matches_supported_route`. |
| `TechnicalReplicateResolver` | `src/phospy/workflows/differential/replicates.py` | Internal workflow compatibility wrapper. | `TechnicalReplicateAggregationPlanner` plus `TechnicalReplicateAggregator`, or normal `DifferentialAnalysisWorkflow.run(...)`. | `DeprecationWarning` on construction and `run()`. | Retain for now; remove in a future release. | `tests/unit/test_differential_technical_replicate_policy.py::test_technical_replicate_resolver_warns_and_preserves_wrapper_behaviour`. |
| `load_enrichment_sets_gmt`, `load_enrichment_sets_table`, `load_enrichment_sets_csv`, `load_enrichment_sets_tsv` | `src/phospy/io/readers/enrichment_sets.py` | Public reader aliases. | Matching `read_enrichment_sets_*` function. | `DeprecationWarning` when called. | Retain for now; remove in a future release. | `tests/unit/test_reader_enrichment_sets.py::test_reader_enrichment_load_aliases_warn_and_forward`. |
| `PeptideToSiteAggregationConfig(strategy="compat_best_p_value")` | `src/phospy/science/differential/aggregation/models.py` | Public lower-level aggregation config value for historical reproduction. | `fixed_effect_meta`, `random_effect_meta`, or `stouffer_z`, chosen for the available peptide-level uncertainty. | `DeprecationWarning` at config construction. | Retain for now; remove in a future release. | `tests/unit/test_peptide_to_site_aggregation.py::test_compatibility_mode_reproduces_minimum_p_value_selection`, `tests/unit/test_peptide_to_site_aggregation.py::test_scientific_policy_metadata_warns_for_compatibility_min_p_mode`. |
| `KinaseActivityResult.activity_scores`, `KinaseActivityResult.weighted_activity` | `src/phospy/science/activities/models.py` | Public result aliases. | `KinaseActivityResult.activity_matrix`. | `DeprecationWarning` on property access. | Retain for now; remove in a future release after bundle/result consumers have migrated. | `tests/unit/test_activity_science.py::test_activity_scores_compatibility_alias_matches_activity_matrix`, `tests/unit/test_activity_science.py::test_weighted_activity_compatibility_alias_matches_activity_matrix`, `tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_kinase_bundle_integration.py`. |
| Bare string entries in `build_motif_library_from_sequences(motif_sequences=...)` | `src/phospy/science/prediction/motif_scoring/library_validation.py` | Semi-public motif-library input compatibility path. | `ExplicitMotifSequence` values or mapping entries with `reference_id`, optional `site_id`, `kinase`, and `sequence`. | `DeprecationWarning` when any bare string entry is present. | Retain for now; reject in a future release. | `tests/unit/test_prediction_sequence_validation.py::test_structured_explicit_sequences_do_not_emit_bare_string_warning`, `tests/unit/test_prediction_sequence_validation.py::test_bare_explicit_sequences_warn_and_remain_supported_without_site_mismatch_claims`, `tests/unit/test_prediction_sequence_validation.py::test_bare_explicit_invalid_sequence_is_rejected_and_reported`. |

### Retained Compatibility Paths Without Deprecation Warnings

| Surface | Module | Status | Replacement or recommended route | Warning | Decision and removal target | Intentional tests |
| --- | --- | --- | --- | --- | --- | --- |
| `phospy.api.configs.*`, `phospy.api.requests`, `phospy.api.results` wrapper modules | `src/phospy/api/` | Public API wrapper routes. | Use `phospy.api` for public request, config, result, enum, reference, workflow, and exception contracts. | No. | Retain; no removal target. | `tests/unit/test_public_contract_workflows.py`, `tests/unit/test_public_contract_results.py`, `tests/unit/test_public_contract_import_routes.py`. |
| Dataset request constants `DATASET_MULTI_SITE_POLICY_*` and `DATASET_SITE_RESOLUTION_MODE_*` | `src/phospy/api/requests.py`, `src/phospy/contracts/requests.py`, `src/phospy/science/evidence/dataset_resolution.py` | Public request compatibility constants. | Use the documented `DatasetBuildRequest` policy fields; constants remain valid for explicit imports. | No. | Retain; no removal target. | `tests/unit/test_public_contract_workflows.py::test_request_compatibility_constants_are_public_exports`, `tests/unit/test_dataset_peptide_evidence_resolution.py`. |
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
| Activity sidecars `thresholded_substrate_mean_activity`, `thresholded_substrate_counts`, and `activity_substrate_counts` | `src/phospy/science/activities/models.py`, `src/phospy/io/bundles/_kinase/`, `src/phospy/io/bundles/_signalome/` | Public result and bundle payload compatibility. | Prefer method-neutral `activity_matrix` and `substrate_count_matrix` for primary activity score outputs. | No. | Retain for current activity methods and bundle round trips; no removal target. | `tests/unit/test_activity_science.py`, `tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_kinase_bundle_integration.py`. |

### Already Removed or Explicitly Rejected Paths

| Surface | Module | Status | Replacement or recommended route | Warning | Decision and removal target | Intentional tests |
| --- | --- | --- | --- | --- | --- | --- |
| Legacy provenance schemas and table/stage hash aliases | `src/phospy/provenance/serialization.py` | Removed; payloads fail clearly. | Regenerate provenance with the current schema. | No; raises `PhosPyInputError`. | Already removed in 1.6.0. | `tests/unit/test_provenance_hashing.py`, `tests/unit/test_dataset_run_provenance.py`. |
| Legacy kinase and signalome bundle schemas, including old diagnostic fields | `src/phospy/io/bundles/_kinase/`, `src/phospy/io/bundles/_signalome/` | Removed; bundle load fails clearly. | Regenerate bundles with the current PhosPy version. | No; raises `PhosPyInputError`. | Already removed in 1.6.0. | `tests/integration/test_kinase_bundle_integration.py`, `tests/integration/test_signalome_bundle_integration.py`, `tests/unit/test_signalome_bundle_schema.py`, `tests/unit/test_signalome_clustering_diagnostics_helpers.py`. |
| Removed public/config aliases: `ratio_to_total`, `ensemble_size`, `signalome_cutoff`, `max_exact_clustering_sites`, `kinase_network_policy`, `tree_engine`, and old backend-style signalome names | Dataset, kinase, and signalome config/snapshot modules | Removed; construction or snapshot load rejects them. | Current names such as `subtract_log_total`, `n_iterations`, `scientific.substrate_support_cutoff`, `performance.max_exact_tree_sites`, `output.network_policy`, and `clustering.clustering_engine`. | No; raises `TypeError` or `PhosPyInputError`. | Already removed in 1.6.0. | `tests/unit/test_public_config_self_validation.py`, `tests/unit/test_kinase_bundle_snapshots.py`, `tests/unit/test_signalome_bundle_schema.py`, `tests/unit/test_validator_boundaries.py`, `tests/integration/test_dataset_builder_integration.py`. |
| Unsupported historical site-metadata aliases (`sequence`, `protein`, `gene`, `residue`, `phosphosite`) used as substitutes for strict builder columns | `src/phospy/science/datasets/builders/site_metadata_normalizer.py` | Removed/rejected. | Rename to supported columns such as `site_sequence`, `protein_id`, `gene_symbol`, or `site`. | No; raises `UnsupportedInputFormatError`. | Already removed in 1.6.0. | `tests/unit/test_domain_boundaries.py`, `tests/unit/test_dataset_normalizer_identity.py`, `tests/unit/test_dataset_normalizer_collaborators.py`. |
| Removed kinase/signalome workflow component re-export modules | `src/phospy/workflows/` | Removed import routes. | Import workflow classes from `phospy`, `phospy.api.workflows`, or `phospy.workflows.<workflow>`. | No; import raises `ModuleNotFoundError`. | Already removed. | `tests/unit/test_workflow_public_import_contract.py::test_removed_workflow_compatibility_imports_fail`, `tests/unit/test_workflow_public_import_contract.py::test_docs_and_examples_do_not_reference_removed_workflow_compatibility_paths`. |

Review rules:

- Warning-based deprecations must name a replacement, emit `DeprecationWarning`,
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
