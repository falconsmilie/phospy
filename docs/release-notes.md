# PhosPy Release Notes

## Version 1.7.0

Release date: 2026-08-14.

These notes describe the changes since Version 1.6.0.

## Release Overview

PhosPy 1.7.0 is a release-hardening and public-contract release. The package now
targets Python 3.11 and 3.12, separates stable and advanced public API routes,
seals ordinary dataset construction behind the builder, and records more of the
scientific contract in typed provenance. The release also expands native
preprocessing correction, differential helper utilities, kinase method
contracts, signalome diagnostics, result caveats, artifact verification, and
release-gate coverage.

## Kinase Scientific-Policy Versions

The current implementation owns these policy and schema versions:

| Policy | Implemented version |
| --- | ---: |
| KSEA activity policy | 5 |
| Membership-selection policy | 4 |
| Inferential policy | 4 |
| Membership payload schema | 2 |
| Membership-independence policy | 2 |

These versions govern the KSEA scientific contract recorded in provenance and
bundles: membership evidence, whether substrate membership was selected
independently of the tested matrix, whether ordinary KSEA p/q output is
eligible, and compatibility for persisted membership and provenance payloads.
They are compatibility and interpretation contract identifiers, not empirical
proof of scientific validity.

## Compatibility and Migration

- BREAKING: Python 3.10 support was removed. PhosPy 1.7.0 requires Python 3.11
  or Python 3.12.
- BREAKING: ordinary `AnalysisReadyPhosphoDataset(...)` construction now raises
  immediately. Use
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))` for normal
  inputs. Use `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` only for
  audited trusted-table replay with complete
  `TrustedDatasetConstructionAssertions`.
- BREAKING: kinase activity statistics tables use `profile_id` as authoritative
  row identity. The legacy condition-shaped statistics adapter remains only as a
  deprecated compatibility helper.
- `KinaseWorkflowRequest.activity_config` defaults to `None`, so activity-like
  summaries are disabled unless callers supply a `KinaseActivityConfig`.
- `KinaseScoringConfig` requires explicit reliability intent through the
  exploratory, production, or custom profile routes.
- Legacy kinase bundle manifest-version 2 payloads are rejected with regenerate
  instructions. Current bundles persist typed activity input semantics and
  profile metadata.
- The post-hoc peptide-to-site differential estimate-combination route is
  withdrawn from public support under
  `unsupported_withdrawn_posthoc_estimate_combination_v1`. The supported lane is
  peptide evidence resolution at sample-intensity level before
  `DifferentialAnalysisWorkflow`.
- `PreprocessingPipeline(stage_metadata_registry=...)` is deprecated; use
  `stage_contract_registry=...`.
- `load_enrichment_sets_gmt`, `load_enrichment_sets_table`,
  `load_enrichment_sets_csv`, and `load_enrichment_sets_tsv` remain deprecated
  in favour of the matching `read_enrichment_sets_*` functions.

## Major Additions

- Governed API stability tiers: `phospy.api` is the stable beta-user route,
  `phospy.advanced` is the supported route for specialist configuration,
  diagnostics, references, and publishing helpers, and implementation modules
  are unsupported import targets.
- Native PhosPy SPS/RUV-style `SpsRuvBatchCorrectionConfig` preprocessing
  correction, including explicit caller-supplied controls,
  missingness policy, factor feasibility checks, selected-control provenance,
  observation-mask provenance, and workflow orchestration.
- Dataset group-coverage filtering, configurable multiple-testing correction,
  differential contrast helpers, differential result filtering/ranking helpers,
  and stronger biological-replicate reliability policy.
- Explicit `paired_design_policy="duplicate_correlation"` support for blocked
  differential designs. This additive option estimates one REML consensus
  compound-symmetry within-block correlation and refits eligible features by
  GLS; `fixed_block` and `reject` remain valid policies, and the default policy
  remains `reject`.
- Common structured `ResultCaveat` records surfaced by differential, kinase,
  signalome, and enrichment workflow results.
- Typed row-attrition provenance across dataset building, differential,
  kinase, signalome, and enrichment workflows, with causal site-row attrition
  separated from compatibility metrics such as site/kinase-pair loss.
- Reference-context compatibility records attached to datasets and workflow
  results, plus configurable handling for unknown or mismatched contexts.
- Typed intensity-transformation events, derived quantitative data provenance,
  quantitative-operation contracts, and evidence sidecars for preprocessing
  traces.
- Kinase scoring support for explicit reliability profiles, profile
  self-inclusion policy, leave-one-out profile scoring, true Kinase
  Library-style motif-only scoring, method-specific quantitative contracts, and
  optional substrate contribution tables.
- Signalome module-selection stability diagnostics, small-sample and
  fully-missing clustering guards, network paired-observation guards, and edge
  skip diagnostics.
- Importer quality reports, shared reader table-parsing helpers, fixed
  fixture-byte policy, and hardened MaxQuant/FragPipe/PTMProphet edge-case
  coverage.
- Reference-bundle validation reports and stricter bundled rat `l6_native`
  exact-snapshot metadata, attribution, manifest, and source-version handling.
- Release-science assets: large-feature R/limma trend parity fixture,
  PhosPy-owned release-validation regression fixture families, minimum
  dependency constraints, installed wheel/sdist verifier, and an opt-in local
  50,000-site by 48-sample benchmark.

## Changed

- Public configuration ownership was consolidated so each config concept has one
  owner. Stable `phospy.api` imports and supported advanced imports are
  preserved through facade modules rather than duplicate class definitions.
- Scientific table schemas were consolidated under `phospy.science.tables`,
  with shared `TableSchema` infrastructure under `phospy.frames`. Supported
  `phospy.tables.*` imports are identity-preserving compatibility re-exports.
- Dataset construction now produces validated aggregate state, recursive
  immutability for exported JSON-like payloads, explicit trusted-construction
  assertions, and stricter organism/reference identity coherence.
- Preprocessing stage construction now uses explicit stage contracts and typed
  resolved sections rather than reflective collaborator negotiation.
- Result bundles use explicit overwrite policy, transactional writes, and
  content-addressed integrity checks.
- Kinase bundle schema advanced to manifest version 3, preserving exact typed
  activity semantics, profile axes, quantitative semantics, identifiers, and
  condition-summary aggregation records.
- Kinase, differential, and signalome workflows reuse validator-owned immutable
  dataset views within a run to reduce redundant DataFrame copies while keeping
  separate workflow runs isolated.
- The release process follows ADR-0039's lightweight solo-maintainer path:
  normal CI/build confidence, fresh wheel/sdist builds, archive-level packaged
  reference validation, installed artifact execution, and trusted publishing
  replace the former formal retained-evidence attestation run.

## Fixes and Hardening

- KSEA ordinary normal-approximation p/q availability is derived and enforced
  from typed membership provenance. Adaptive, profile-derived, fused
  profile/motif, leave-one-out, unknown, incomplete, or tampered membership
  cannot reconstruct finite ordinary p/q output.
- ssGSEA-style substrate-enrichment activity is invariant to equal-valued row
  ordering by scoring ties as method-owned midrank expectation blocks.
- Quantitative-meaning relabelling now requires authority-gated semantic
  provenance. Caller declarations are limited to supported meanings, and
  derived total-protein meanings carry operation fingerprints.
- Preprocessing `QuantitativeOperationContract.required_evidence` is executable
  policy. Stages fail before trace acceptance when required evidence is missing,
  and processing-state reconstruction revalidates typed trace evidence.
- KNN imputation records column-mean fallback as typed scientific provenance and
  preserves explicit missingness/ordering semantics.
- Peptide evidence resolution rejects conflicting accessions or sequence
  contexts deterministically and keeps non-executable mapping policies from
  silently running as ordinary evidence.
- Reference-bundle release hardening validates exact source-tree and built-wheel
  bytes, rejects unknown manifest/evidence fields, rejects non-Boolean or JSON
  null raw `redistribution_allowed` values, requires explicit `verified_at`
  dates for approved bundled evidence, and treats hashes as integrity checks
  rather than redistribution approval.
- Release validation now installs both wheel and sdist in isolated environments
  outside the checkout, checks installed import origins, validates bundled rat
  reference resources against manifest SHA-256 values, and runs representative
  public dataset, differential, and kinase workflow contracts.
- The 50,000-site by 48-sample workload is no longer a release-blocking
  pytest/CI performance contract. It is an explicitly invoked local benchmark
  with runtime, output, missingness, fingerprint/provenance, and RSS reporting.
- `make release-check` now includes release/golden selector coverage,
  threshold-bearing parity, performance contracts, strict documentation build,
  archive validation, and installed wheel/sdist verification so release gates
  cannot silently omit required nodes.
- DataFrame ownership, immutable scientific containers, provenance hashing,
  package dependency DAG checks, strict Pyright suppressions, and public
  boundary adversarial tests were strengthened.

## Scientific Scope

Bundled runtime references remain rat-only for `ReferencePreset.AUTO`. The only
approved packaged runtime reference is the exact rat `l6_native` snapshot
derived from upstream PhosR 1.20.0 package data. Its approval is scoped only to
the exact packaged files in the committed PhosPy snapshot. It does not approve
future bundles, other rat bundles, other organisms, or future PhosR/PhosPy
snapshots, and it does not claim independent direct permission from
PhosphoSitePlus, PRIDE, Kinase Library, or another upstream scientific
database.

Human and mouse remain valid organisms, but workflows require an explicit
caller-supplied `ReferenceBundle` unless a future release adds approved
redistributable packaged data.

`site_key` is the analysis-ready row identity. `display_id` remains available
for interpretation and reporting, but it is not the row key and may repeat
across distinct protein contexts.

Kinase prediction, Kinase Library motif scoring, KSEA-style activity, and
ssGSEA-style substrate enrichment are explicit PhosPy workflow methods. They do
not claim calibrated causal kinase inference, full PhosR kinase-activity
equivalence, validated Kinase Library parity, or PTM-SEA support.

Enrichment support is offline ORA over caller-supplied collections and explicit
backgrounds. It does not bundle GO, KEGG, Reactome, PTM-SEA, PTMsigDB, Enrichr,
gseapy, clusterProfiler, GSEA, or online-service behaviour. Typed
selected/background identifier-set provenance is optional and does not change
enrichment statistics.

Fixed-effect batch/covariate/block differential designs, explicit
`duplicate_correlation` paired differential designs, linear batch
residualisation, and native SPS/RUV-style preprocessing correction are
executable PhosPy features. `duplicate_correlation` uses one REML-estimated
consensus compound-symmetry correlation and GLS; it is not a general
mixed-effects model or feature-specific random-effects fit. The preprocessing
batch-correction methods are not ComBat, PhosR-equivalent RUV/SPS/RUV-III, or
limma `removeBatchEffect` parity claims.

Committed differential limma parity fixtures are implementation evidence for
the exact fixture-scoped model envelopes they cover. They are not independent
scientific validation and do not imply general limma equivalence.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
