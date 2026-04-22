# Fixtures

This page is the fixture governance truth source for active, archival, and
removable material.

> Audience: contributors and maintainers working on parity, provenance, and
> fixture regeneration.

For parity intent and science coverage tiers, see [`parity.md`](parity.md).

## Classification Model

Every committed fixture/output family must be explicitly classified as one of:

- `active`: used by active tests or supported rewrite reference workflows
- `archival`: retained only for provenance/forensics, not part of routine
  maintainer lanes
- `removable`: dead weight with no active test/provenance role; should be
  deleted

## Active Fixture Lanes

Canonical active roots:

- `tests/fixtures/rewrite_parity/`
- `tests/fixtures/public_workflow_reference/`

Default supported maintainer generators:

- `scripts/active/generate_r_l6_fixtures.R` ->
  `tests/fixtures/rewrite_parity/r_reference_l6`
- `scripts/active/generate_signalome_public_workflow_reference.py` ->
  `tests/fixtures/public_workflow_reference`

Primary active bootstrap lane:

```bash
make fixtures-all
```

`fixtures-all` intentionally regenerates active rewrite/public workflow fixture
families only.

## Legacy-Script Output Classification Register (2026-04-22)

| Family | Legacy/archived script lineage | Classification | Current role |
| --- | --- | --- | --- |
| `tests/fixtures/rewrite_parity/r_reference_l6/` | activity tables were materialized with archived `scripts/archive/generate_activity_donor_snapshot.py` from promoted donor inputs | `active` | activity-stage parity gate (`tests/parity/test_activity_stage_parity.py`) |
| `tests/fixtures/rewrite_parity/r_reference_l6_prediction/` | promoted donor prediction outputs with rewrite refresh | `active` | L6 prediction parity gate (`tests/parity/test_l6_prediction_parity.py`) |
| `tests/fixtures/rewrite_parity/fragile_support_reference/` | regeneration utility archived at `scripts/archive/generate_fragile_support_reference.py` | `active` | prediction-science parity gate (`tests/parity/test_prediction_science_parity.py`) |
| `tests/fixtures/rewrite_parity/adaptive_sampling_edge/` | originated from archived synthetic seam tooling (`scripts/archive/generate_synthetic_adaptive_sampling_edge_fixtures.py`) | `active` | adaptive prediction parity gate (`tests/parity/test_adaptive_prediction_parity.py`) |
| `tests/fixtures/rewrite_parity/adaptive_sampling_replay/` | promoted donor replay traces from legacy L6 prediction trace outputs | `active` | adaptive replay parity gate (`tests/parity/test_adaptive_replay_parity.py`) |
| `tests/fixtures/rewrite_parity/protein_correction/` | promoted legacy donor preprocessing outputs | `active` | preprocessing parity gate (`tests/parity/test_preprocessing_science_parity.py`) |
| `tests/fixtures/rewrite_parity/site_matrix/` | promoted legacy donor site-matrix outputs | `active` | preprocessing parity gate (`tests/parity/test_preprocessing_science_parity.py`) |
| `tests/fixtures/rewrite_parity/comparison_building/` | promoted legacy donor comparison-building outputs | `active` | preprocessing parity gate (`tests/parity/test_preprocessing_science_parity.py`) |
| `tests/fixtures/public_workflow_reference/` | includes rewrite-owned public workflow outputs plus donor benchmark context tables | `active` | public workflow parity gates (`tests/parity/test_public_predmat_parity.py`, `tests/parity/test_signalome_workflow_parity.py`) |
| `tests/fixtures/archive/adaptive_sampling_edge_trace_debug/` | non-gated seam-debug trace tables moved out of active lane on 2026-04-22 | `archival` | provenance/forensics only |
| `tests/fixtures/archive/adaptive_sampling_replay_trace_debug/` | non-gated replay debug trace tables moved out of active lane on 2026-04-22 | `archival` | provenance/forensics only |
| `tests_legacy/fixtures/python_reference_l6/prediction_trace*/` | outputs from archived python prediction trace export workflows | `archival` | historical seam/probability-trace forensics only |
| `tests_legacy/fixtures/r_reference/` | historical small R fixture generator (`scripts/archive/generate_r_fixtures.R`) | `archival` | legacy parity forensics only |
| `tests_legacy/fixtures/r_reference_l6/` | historical donor L6 reference tree used as promotion source | `archival` | provenance source only; active tests consume promoted rewrite-owned copies |
| `tests_legacy/fixtures/r_reference_l6_seam_stress/` | historical seam-stress generator (`scripts/archive/generate_l6_seam_stress_reference.py`) | `archival` | seam-forensics only |
| `tests_legacy/fixtures/synthetic_adaptive_sampling_edge/` | historical synthetic adaptive seam-debug source tree | `archival` | provenance source only |
| `tests_legacy/fixtures/fragile_support_reference/` | historical fragile-support seam-debug source tree | `archival` | provenance source only |
| `tests_legacy/fixtures/public_workflow_reference/` | historical donor public-workflow outputs | `archival` | provenance source only |
| `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*_selected.csv` | older debug slice exports | `removable` | deleted on 2026-04-22 (not used by active tests/docs) |
| probability-diff report outputs (`scripts/archive/diff_prediction_trace_probabilities.py`) | debug report artifacts from archived trace-diff tooling | `removable` | no committed outputs are retained in active fixture roots |

## Maintainer Targets and Archive Policy

Repository script layout is status-first:

- `scripts/active/`: supported current maintainer generators
- `scripts/support/`: helpers used by active generators
- `scripts/archive/`: historical parity/debug tooling

Archived/niche target retained for explicit forensics:

- `make fixtures-r-small-archive`

This target is archival only and excluded from `make fixtures-all`.

## Archive Notes

Archival fixture trees:

- `tests/fixtures/archive/`
- `tests_legacy/fixtures/`

Archive material is retained for provenance and audit workflows, not routine
rewrite maintenance.

## Where Next

- Parity governance and confidence tiers: [Parity to PhosR](parity.md)
- Maintainer navigation hub: [Contributor and maintainer docs](contributor/index.md)
