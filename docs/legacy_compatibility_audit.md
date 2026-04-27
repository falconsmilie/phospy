# Legacy/Compatibility Audit

Date: 2026-04-27

Scope: repository search for legacy/rewrite/migration/compatibility/deprecated/alias/canonical/backwards/prerelease markers, then triage into keep/remove/replace/follow-up.

## Keep

- `docs/adr/**`: historical design records.
- `CHANGELOG.md` and `docs/release_notes/**`: historical release records.
- Domain-meaningful `compatibility` wording in runtime validation (for example dataset/reference organism compatibility checks).

## Remove (implemented)

- `SignalomeConfig` deprecated alias fields:
  - removed `clustering_backend`
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
  - removed `clustering_backend`
  - removed `max_exact_clustering_sites`
  - removed deprecated alias fields from emitted snapshot/provenance payloads.
- Old compatibility tests that enforced removed behavior:
  - replaced with explicit rejection tests.

## Replace (implemented)

- User-facing namespace wording:
  - replaced “canonical namespace” with “primary namespace”.
- User-facing docs wording:
  - removed deprecated/legacy framing tied to removed aliases.
  - replaced “legacy aliases” with “historical aliases” where docs describe rejected inputs.

## Follow-Up (new tickets)

- `LEGACY-001` Bundle schema hardening:
  - decide whether to remove remaining old-manifest load tolerances (`manifest` without provenance, partial old metadata variants) after confirming real public artifact usage.
- `LEGACY-002` Historical docs cleanup:
  - decide whether to retire or archive rewrite/legacy-heavy maintainer docs (`docs/parity.md`, `docs/roadmap.md`, selected `docs/architecture/*`) outside user-facing lane.
- `LEGACY-003` Namespace cleanup:
  - consider renaming internal modules that still use compatibility-era names (for example `_signalome/compatibility.py`) now that compatibility fallbacks are removed.
