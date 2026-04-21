# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. The rewrite does not
claim full package equivalence with PhosR.

This page is the project truth source for scientific regression confidence
tiers. In this repository, `implemented`, `supported`, `parity-gated`, and
`closed` are not interchangeable terms.

## Scope and Contract

Scoped parity passes in this document do not equal whole-package legacy parity.
Audit boundaries, explicit exclusions, and area-level evidence are tracked in:
`docs/architecture/legacy_science_gap_audit.md`.

## Coverage Tier Vocabulary (Normative)

Every scientific lane should be described with one of these coverage tiers:

- `PARITY_GATED_ACTIVE_SCIENCE`: rewrite-owned behavior guarded by active
  parity-focused tests in `tests/parity/` and treated as the highest regression
  confidence tier in this project. Promotion to this tier requires explicit
  `tests/parity/...` test evidence in the legacy-science audit inventory.
- `DONOR_BACKED_REWRITE_COVERAGE`: rewrite-implemented behavior supported by
  rewrite-owned unit/integration coverage and donor-informed fixtures/evidence,
  but not promoted to the same parity-gated tier.
- `CONTRACT_CHANGED_SUPPORTED_LANE`: rewrite behavior that is intentionally
  supported under a narrowed or reshaped contract relative to legacy behavior.
  This lane can still have strong tests (including parity tests), but should
  not be described as broad legacy-equivalent parity.
- `OPEN_SCIENTIFIC_GAP`: unresolved area where science decisions, parity
  decisions, or regression confidence are still insufficient for closure.

Legacy inventory status labels are a separate axis and remain:
`PORTED`, `CONTRACT_CHANGED`, `OPEN_GAP`, `INTENTIONALLY_RETIRED`.
Status labels describe governance state; coverage tiers describe confidence and
regression protection strength.

## Promotion Guardrails (Normative)

- Every new or changed science area must be added to the legacy-science
  inventory with both `Status` and `Coverage tier`.
- `implemented`, `PORTED`, and closed ticket labels are not enough to claim
  parity closure.
- `PARITY_GATED_ACTIVE_SCIENCE` claims require active rewrite-owned parity tests
  under `tests/parity/` plus rewrite-owned fixture/provenance evidence.
- If that gate evidence does not exist yet, classify the area as
  `DONOR_BACKED_REWRITE_COVERAGE`,
  `CONTRACT_CHANGED_SUPPORTED_LANE`, or `OPEN_SCIENTIFIC_GAP`.
- When parity gates are added, removed, or materially changed, update all of:
  `docs/parity.md`, `docs/architecture/legacy_science_gap_audit.md`, and
  `tests/support/legacy_donor_inventory.py` in the same change.

## What Parity Means Here

Parity in this repository is:

- seam-level
- selective
- tied to committed fixtures
- strongest where lanes are explicitly classified as
  `PARITY_GATED_ACTIVE_SCIENCE`

Parity here does not mean:

- every PhosR feature is implemented
- every Python path must numerically match PhosR
- every `PORTED` row has the same regression confidence tier

## Active Parity-Gated Science

The following areas currently run as active parity-focused science in supported
rewrite lanes (`PARITY_GATED_ACTIVE_SCIENCE`):

- prediction-science parity on the fragile-support rewrite fixture lane
- kinase workflow parity on the supported L6 rewrite lane
- full promoted L6 downstream prediction/scoring parity against rewrite-owned
  promoted reference tables (`profile_scores`, `combined_scores`, `weights`,
  candidate substrates, prediction-matrix ranking summaries, top-k export
  ranking summaries)
- release gates for the core kinase lane, including strict candidate overlap,
  prediction-matrix ranking agreement thresholds, top-k export agreement
  thresholds, and replay-surface agreement in rewrite-owned fixture lanes
- adaptive prediction parity from promoted adaptive-sampling fixtures, executed
  in both supported rewrite policy lanes:
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`
- adaptive replay-trace parity from promoted replay fixtures, including initial
  negative-set replay surfaces, per-iteration sample-membership replay
  surfaces, final ensemble probabilities, top-k replay summaries, and
  stable-vs-r_parity comparison metrics under fixed seed
- public end-to-end `predMat` benchmark parity for both supported adaptive
  policies
- public end-to-end `predMat` order-invariance parity on the stable/default
  lane
- activity-stage parity outputs from fixed `predMat` + phospho inputs
- full-table signalome regression contracts on the supported L6 lane:
  `module_assignments`, `signalome_modules`, `kinase_network.nodes`,
  `kinase_network.edges`, `expanded_signalome`
- preprocessing-science parity on supported builder lanes:
  `total_protein_correction.policy="ratio_to_total"`,
  `site_matrix.policy="build_from_metadata"`, and
  `comparisons.policy="sample_metadata_pairs"`

## Donor-Backed Rewrite Coverage (Not Parity-Gated)

No lanes in the current audited inventory snapshot dated `2026-04-21` are
classified in this tier.

## Contract-Changed Supported Lanes

These lanes are intentionally supported but not legacy-equivalent by contract:

- adaptive sampling public contract uses `adaptive_policy` naming (not legacy
  `svm_mode`)
- builder transformation establishment is intentionally narrow
  (`transformation_state.label == "linear"` via identity establishment path)
- signalome requires explicit `site_metadata.protein_id` (no legacy
  site-id-prefix fallback)
- signalome public entrypoint is contracted to
  `SignalomeWorkflowRequest(kinase_result=...)`
- motif sequence authority in supported kinase lane is
  `references.site_sequences` (not dataset-sequence fallback)

## Open Scientific Gaps

`OPEN_SCIENTIFIC_GAP` remains an active vocabulary tier and should be used
whenever evidence is incomplete or unresolved. In the current audited inventory
snapshot dated `2026-04-21`, no rows are classified `OPEN_GAP`.

This is not a blanket claim that every possible legacy-science surface is
closed; it applies only to the explicit inventory and boundaries documented
here and in the architecture audit.

## Core Kinase Lane Status (2026-04-20)

The central kinase scoring/prediction lane is treated as closed at the parity
tier only when both rewrite parity gates pass:

- `tests/parity/test_l6_prediction_parity.py`
- `tests/parity/test_adaptive_replay_parity.py`

These gates enforce downstream behavior (candidate selection, prediction-matrix
ranking agreement, top-k export agreement, and adaptive replay surfaces)
against promoted rewrite fixture references.

## Rewrite-Owned Parity Reporting

Parity chatter is emitted by default from `tests/parity/`. Reporting is
rewrite-owned (`tests/support/` + `tests/conftest.py`) and is not routed
through `tests_legacy/`.

When parity tests run, terminal output includes grouped scientific summaries
for prediction-science parity, kinase workflow parity, L6 core scoring/prediction,
adaptive prediction, adaptive replay-trace parity, public end-to-end predMat
parity, predMat order-invariance parity, activity-stage parity,
preprocessing-science parity, and signalome workflow parity.

No `PHOSPY_SHOW_*` environment variables are required.

## Active Parity Gate Files (2026-04-21)

Active parity-gated science currently resolves through the following
rewrite-owned parity files:

- `tests/parity/test_prediction_science_parity.py`
- `tests/parity/test_kinase_workflow_parity.py`
- `tests/parity/test_l6_prediction_parity.py`
- `tests/parity/test_adaptive_prediction_parity.py`
- `tests/parity/test_adaptive_replay_parity.py`
- `tests/parity/test_public_predmat_parity.py`
- `tests/parity/test_activity_stage_parity.py`
- `tests/parity/test_preprocessing_science_parity.py`
- `tests/parity/test_signalome_workflow_parity.py`

## Legacy Science Coverage Inventory

This table is the parity-governance truth source for tracked legacy-science
areas. `Status` and `Coverage tier` are intentionally separate columns so that
`PORTED` does not imply parity-gated closure.

| Legacy science area | Status | Coverage tier | Contract relation | Rewrite coverage summary |
| --- | --- | --- | --- | --- |
| profile policy behavior | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `strict` + `median_skipna` profile behavior is supported and parity-tested. |
| core kinase scoring/prediction lane | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | Candidate/ranking/replay behavior is parity-gated in rewrite-owned tests. |
| adaptive sampling / svm_mode | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (`adaptive_policy` replaces legacy `svm_mode` naming) | Adaptive science is implemented with parity evidence, but contract naming intentionally differs from legacy. |
| signalome clustering/module selection | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | Clustering and module-count diagnostics are implemented and parity-backed. |
| weighted-top assignment behavior | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | Weighted-top assignment and fractional support propagation are implemented and parity-backed. |
| network policy variants | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `positive_only`, `absolute_threshold`, and `signed` are implemented and parity-tested. |
| expanded signalome outputs | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `expanded_signalome` is materialized in the supported workflow path and parity-tested. |
| activity parity lock | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | Activity/KSEA science is rewrite-ported and guarded by parity CI gates. |
| preprocessing transformation establishment | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (narrow builder establishment policy) | Supported builder lane establishes only `linear` pass-through transformation state. |
| total/protein correction | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `total_protein_correction.policy="ratio_to_total"` is parity-gated in rewrite-owned fixture tests, including strict phospho/total alignment behavior. |
| site-matrix construction | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | Supported site-matrix construction (`build_from_metadata`) is parity-gated with rewrite-owned fixture expectations for row retention, site identity, and output matrix values. |
| comparison-building | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | Sample-metadata comparison construction is parity-gated for explicit and inferred pair policies, including pair identity/order and expected output values. |
| site-to-protein resolution fallback behavior | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (no legacy site-id-prefix fallback) | Signalome requires explicit `site_metadata.protein_id`. |
| signalome input route contraction | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (workflow entrypoint narrowed) | Supported signalome entrypoint is `SignalomeWorkflowRequest(kinase_result=...)`. |
| dataset-vs-reference sequence authority decisions | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (reference bundle is sequence authority) | Motif sequence authority in supported kinase lane is `references.site_sequences`. |

Open legacy-science areas in this inventory snapshot (`2026-04-21`):
- none currently classified `OPEN_GAP` in the audited list.

Rewrite-side visibility check:

- `tests/unit/test_legacy_donor_inventory.py`

## Fixture Locations

### Rewrite-owned parity inputs and expectations

- `tests/fixtures/rewrite_parity/r_reference_l6/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/r_reference_l6_prediction/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/r_reference_l6_prediction/PROVENANCE.md`
- bundled donor motif tables used by the supported rat L6 native bundle:
  `src/phospy/data/reference_bundles/rat/l6_native/motif_scores.csv`,
  `src/phospy/data/reference_bundles/rat/l6_native/motif_sizes.csv`
- `tests/fixtures/rewrite_parity/fragile_support_reference/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/fragile_support_reference/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/adaptive_sampling_edge/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/adaptive_sampling_edge/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/adaptive_sampling_replay/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/adaptive_sampling_replay/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/protein_correction/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/protein_correction/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/comparison_building/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/comparison_building/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/site_matrix/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/site_matrix/PROVENANCE.md`

These files are the normal source for active parity tests in `tests/parity/`
and helpers in `tests/support/rewrite_fixture_data.py`.

### Rewrite workflow regression expectations

- public predMat benchmark inputs and committed rewrite outputs:
  `tests/fixtures/public_workflow_reference/predmat_input_*.{csv,json}`,
  `tests/fixtures/public_workflow_reference/predmat_rewrite_*.csv`,
  `tests/fixtures/public_workflow_reference/predmat_rewrite_contract.json`
- provenance and promotion history:
  `tests/fixtures/public_workflow_reference/PROVENANCE.md`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*.csv`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json`

### Historical reference archive

- `tests_legacy/fixtures/` is retained for provenance and archival material.
- Active rewrite parity tests should not resolve fixtures from this path as
  their normal source.
- `tests_legacy/test_parity-with_metrics.py` is archival/provenance only and is
  not active reporting infrastructure for rewrite parity runs.

## Run the Parity Suite

```bash
pytest tests/parity -m parity -s
```

or:

```bash
make test-parity
```

Useful variants:

```bash
pytest tests/parity -m parity -rs -s
pytest -m parity -vv
pytest tests/parity/test_activity_stage_parity.py -m "parity and activity_parity" -vv
pytest tests/parity/test_signalome_workflow_parity.py -vv
```
