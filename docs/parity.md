# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. The rewrite does not
claim full package equivalence with PhosR.

Scoped parity passes in this document do not equal full legacy-science parity.
Legacy-science coverage status is tracked explicitly in the inventory below.

## What Parity Means Here

Parity in this repository is:

- seam-level
- selective
- tied to committed fixtures

Parity here does not mean:

- every PhosR feature is implemented
- every Python path must numerically match PhosR

## Active Parity Coverage

The parity suite currently protects rewrite-era parity families for:

- prediction-science parity on the fragile-support rewrite fixture lane
- kinase workflow parity on the supported L6 rewrite lane
- full promoted L6 downstream prediction/scoring parity against rewrite-owned
  promoted donor references (`profile_scores`, `combined_scores`, `weights`,
  candidate substrates, ranking/top-k summaries)
- legacy-grade release gates for the core kinase lane, including strict
  candidate overlap, ranking agreement thresholds, and replay-surface
  agreement (no permissive donor-similarity fallback in this lane)
- adaptive prediction parity from promoted adaptive-sampling fixtures, executed
  in both supported rewrite policy lanes:
  `adaptive_policy="stable"` (default lane) and
  `adaptive_policy="r_parity"`
- adaptive replay-trace parity from promoted replay fixtures, including:
  initial negative-set replay surfaces, per-iteration sample-membership replay
  surfaces, final ensemble probabilities, top-k replay summaries, and
  stable-vs-r_parity comparison metrics under fixed seed
- public end-to-end `predMat` benchmark parity on the rewrite workflow path for
  both supported adaptive policies:
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`
- public end-to-end `predMat` order-invariance parity on the stable/default lane
  (normalized equality under reference-map order perturbation)
- activity-stage outputs from fixed `predMat` + phospho inputs
- full-table signalome regression contracts on the supported L6 lane:
  `module_assignments`, `signalome_modules`, `kinase_network.nodes`,
  `kinase_network.edges`, `expanded_signalome`

## Analysis-Ready Transformation Coverage Status

Legacy donor preprocessing produced analysis-ready phosphosite matrices before
workflow execution. In the supported rewrite builder lane, transformation
establishment is currently a contract-changed, narrow pass-through policy:

- `AnalysisReadyDatasetBuilder.run(...)` preserves provided quantitative matrix values
- builder establishes `transformation_state.label == "linear"` through the
  supported identity transformer path
- no broader legacy-style transformation-selection API is currently supported

## Core Kinase Lane Status (2026-04-20)

The central kinase scoring/prediction lane is closed only when both rewrite
parity gates pass:

- `tests/parity/test_l6_prediction_parity.py`
- `tests/parity/test_adaptive_replay_parity.py`

These gates enforce legacy-grade downstream behavior (candidate selection,
ranking/top-k agreement, and adaptive replay surfaces) against promoted donor
fixtures.

## Rewrite-owned parity reporting

Parity chatter is emitted by default from the active rewrite suite under
`tests/parity/`. The reporting layer is rewrite-owned (`tests/support/` +
`tests/conftest.py`) and is not routed through `tests_legacy/`.

When parity tests run, terminal output includes grouped scientific summaries for:

- prediction-science parity
- kinase workflow parity
- L6 core kinase scoring/prediction parity
- adaptive prediction parity
- core kinase adaptive replay-trace parity
- public end-to-end predMat parity
- public predMat order-invariance parity
- activity-stage parity
- signalome workflow parity

No `PHOSPY_SHOW_*` environment variables are required.

Adaptive policy comparison is part of the active rewrite parity output:

- both supported rewrite policies execute in parity:
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`
- terminal chatter prints both lanes with clear policy labeling for review:
  `stable (default)` and `r_parity`
- side-by-side comparison metrics are printed in the adaptive parity section
- `svm_mode` remains archival naming and is not a rewrite public API field

Activity parity checks cover:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

Activity parity is a hard regression gate in rewrite CI:

- dedicated job: `activity-parity-gate`
- required marker selection: `parity and activity_parity`
- fixture source pinned to `tests/fixtures/rewrite_parity/r_reference_l6/`
  with provenance in
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`
- active parity assertions compare rewrite runtime outputs to committed
  rewrite-owned fixture expectations; no live `legacy_archive` execution is part
  of this gate

This lane is supported and parity-backed, not provisional.

## Legacy Science Coverage Inventory

This inventory is the parity-governance truth source for what legacy science is
ported, open, or contract-changed in supported rewrite lanes.

Status vocabulary:

- `PORTED`: implemented in supported rewrite lanes and guarded by rewrite-owned tests.
- `INTENTIONALLY_RETIRED`: intentionally unsupported legacy area.
- `OPEN_GAP`: not yet ported into supported rewrite lane.
- `CONTRACT_CHANGED`: rewrite intentionally narrows/reshapes behavior vs legacy.
- The current inventory has no `INTENTIONALLY_RETIRED` rows.

| Legacy science area | Status | Science-gap ticket | Rewrite coverage summary |
| --- | --- | --- | --- |
| profile policy behavior | PORTED | `SCI-GAP-01` | `strict` + `median_skipna` profile behavior is supported and parity-tested. |
| core kinase scoring/prediction lane | PORTED | `SCI-GAP-12` | Candidate/ranking/replay behavior is parity-gated in rewrite-owned tests. |
| adaptive sampling / svm_mode | CONTRACT_CHANGED | `SCI-GAP-05` | Adaptive science is ported, but public contract uses `adaptive_policy` instead of legacy `svm_mode` naming. |
| signalome clustering/module selection | PORTED | `SCI-GAP-06` | Clustering and module-count diagnostics are implemented and parity-backed. |
| weighted-top assignment behavior | PORTED | `SCI-GAP-08` | Weighted-top assignment and fractional support propagation are implemented. |
| network policy variants | PORTED | `SCI-GAP-09` | `positive_only`, `absolute_threshold`, and `signed` are implemented and tested. |
| expanded signalome outputs | PORTED | `SCI-GAP-10` | `expanded_signalome` is materialized in the supported workflow path. |
| activity parity lock | PORTED | `SCI-GAP-11` | Activity/KSEA science is rewrite-ported and guarded by parity CI gates. |
| preprocessing transformation establishment | CONTRACT_CHANGED | - | Builder preprocessing is intentionally narrow: pass-through linear transformation establishment plus limited missing-data policy. |
| total/protein correction | PORTED | - | `total_protein_correction.policy="ratio_to_total"` is supported in builder preprocessing with strict phospho/total matching checks. |
| site-matrix construction | OPEN_GAP | - | Legacy site-matrix construction policy surface is not fully ported. |
| comparison-building | PORTED | - | Builder preprocessing supports sample-metadata-based pairwise comparison construction with explicit or inferred pairs. |
| site-to-protein resolution fallback behavior | CONTRACT_CHANGED | - | Signalome requires explicit `site_metadata.protein_id` and does not apply legacy fallback to site-id prefix. |
| signalome input route contraction | CONTRACT_CHANGED | - | Supported signalome entrypoint is contracted to `SignalomeWorkflowRequest(kinase_result=...)`. |
| dataset-vs-reference sequence authority decisions | CONTRACT_CHANGED | - | Motif sequence authority in supported kinase lane is `references.site_sequences`, not dataset-sequence fallback. |

Open legacy-science areas in this inventory:

- `site-matrix construction`

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
- Active rewrite parity tests should not resolve fixtures from this path as their
  normal source.
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
