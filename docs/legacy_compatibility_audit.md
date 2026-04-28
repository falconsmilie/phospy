# Legacy/Compatibility Audit

Date: 2026-04-27

Scope: repository search for legacy/rewrite/migration/compatibility/deprecated/alias/canonical/backwards/prerelease markers, then triage into keep/remove/replace/follow-up.

## Keep

- `docs/adr/**`: historical design records.
- `CHANGELOG.md` and `docs/release_notes/**`: historical release records.
- Domain-meaningful `compatibility` wording in runtime validation (for example dataset/reference organism compatibility checks).

## Remove (implemented)

- `SignalomeConfig` deprecated alias fields:
  - removed `clustering_engine`
  - removed `max_exact_clustering_sites`
- Signalome deprecated alias exports:
  - removed clustering-backend alias constants/types from `phospy.api.configs`.
- CLI deprecated flags:
  - removed `--prediction-ensemble-size`
  - removed `--clustering-backend`
  - removed `--max-exact-clustering-sites`
- Kinase prediction alias:
  - removed constructor alias `ensemble_size`.
- Dataset preprocessing alias:
  - removed `ratio_to_total` policy alias and resolver plumbing.
- Signalome bundle prerelease/internal payload fallbacks:
  - removed `signalome_cutoff`
  - removed `kinase_network_policy`
  - removed `clustering_engine`
  - removed `max_exact_clustering_sites`
  - removed deprecated alias fields from emitted snapshot/provenance payloads.
- Old compatibility tests that enforced removed behavior:
  - replaced with explicit rejection tests.
- `LEGACY-001` bundle schema hardening:
  - removed unsupported old-manifest tolerances and fallback defaults.
  - manifest loading now requires explicit `provenance`, strict v1 markers, explicit activity enabled markers, and full declared table keys.
  - signalome bundle loading now requires declared `kinase_network_candidate_correlations` marker (path or `null`).
  - config snapshot and signalome diagnostics payload decoders no longer fill missing fields from defaults.
  - unsupported old/partial manifests fail clearly and instruct users to regenerate bundles with the current PhosPy version.

## Replace (implemented)

- User-facing namespace wording:
  - replaced "canonical namespace" with "primary namespace".
- User-facing docs wording:
  - removed deprecated/legacy framing tied to removed aliases.
  - replaced "legacy aliases" with "historical aliases" where docs describe rejected inputs.

## LEGACY-001 Audit Classification

Bundle-loading tolerance findings were classified as follows:

- Missing top-level `provenance` (missing key or `null`):
  - unsupported prerelease/internal format; removed.
- Missing manifest contract markers (`bundle_type`, `manifest_version`, `table_format`, required sections):
  - unsupported prerelease/internal format; removed.
- Partial old metadata variants (for example missing activity `enabled` marker, missing declared optional table markers):
  - unsupported prerelease/internal format; removed.
- Snapshot/diagnostics fallback defaults used for partial payloads:
  - unsupported prerelease/internal format; removed.
- Required by real public artifacts:
  - none identified in repository fixtures or current integration evidence.

## Follow-Up (new tickets)

- `LEGACY-002` Historical docs cleanup:
  - decide whether to retire or archive rewrite/legacy-heavy maintainer docs (`docs/parity.md`, `docs/roadmap.md`, selected `docs/architecture/*`) outside user-facing lane.
- `LEGACY-003` Namespace cleanup:
  - consider renaming internal modules that still use compatibility-era names (for example `_signalome/compatibility.py`) now that compatibility fallbacks are removed.
