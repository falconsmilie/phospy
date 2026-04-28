# PhosR Parity

PhosPy’s PhosR parity scope is intentionally narrow and fixture-backed. It does
not claim whole-package equivalence with PhosR.

> Audience: advanced users and maintainers tracking scientific confidence and
> governance tiers.
>
> If you are onboarding, start with [Getting started](getting-started/index.md)
> and return here when you need PhosR comparison and governance detail.

This page is the project truth source for scientific regression confidence
tiers. In this repository, `implemented`, `supported`, PhosR parity-gated, and
`closed` are not interchangeable terms.

## Scope and Contract

Scoped PhosR parity passes in this document do not equal whole-package PhosR
equivalence.

Audit boundaries, explicit exclusions, and area-level evidence are tracked in:

- [Scientific coverage audit](architecture/science_gap_audit.md)

## Coverage Tier Vocabulary

Every scientific lane should be described with one of these coverage tiers:

- `PHOSR_PARITY_ACTIVE_SCIENCE`: PhosPy-owned behaviour guarded by active PhosR
  parity tests in `tests/parity/` and treated as the highest regression
  confidence tier in this project. Promotion to this tier requires explicit
  `tests/parity/...` test evidence in the science audit inventory.
- `PHOSPY_VALIDATED_COVERAGE`: PhosPy-implemented behaviour supported by
  PhosPy-owned unit/integration coverage and fixture evidence, but not promoted
  to the same PhosR parity tier.
- `CONTRACT_CHANGED_SUPPORTED_LANE`: PhosPy behaviour that is intentionally
  supported under a narrowed or reshaped contract relative to PhosR behaviour.
  This lane can still have strong tests, including PhosR parity tests, but
  should not be described as broad PhosR-equivalent behaviour.
- `OPEN_SCIENTIFIC_GAP`: unresolved area where science decisions, PhosR
  comparison decisions, or regression confidence are still insufficient for
  closure.

Inventory status labels are a separate axis and remain:

- `PHOSR_ALIGNED`
- `CONTRACT_CHANGED`
- `OPEN_GAP`
- `INTENTIONALLY_RETIRED`

Status labels describe governance state. Coverage tiers describe confidence and
regression protection strength.

## Promotion Guardrails

- Every new or changed science area must be added to the science inventory with
  both `Status` and `Coverage tier`.
- `implemented`, `PHOSR_ALIGNED`, and closed ticket labels are not enough to
  claim PhosR parity closure.
- `PHOSR_PARITY_ACTIVE_SCIENCE` claims require active PhosPy-owned PhosR parity
  tests under `tests/parity/` plus fixture/provenance evidence.
- If that gate evidence does not exist yet, classify the area as
  `PHOSPY_VALIDATED_COVERAGE`, `CONTRACT_CHANGED_SUPPORTED_LANE`, or
  `OPEN_SCIENTIFIC_GAP`.
- When PhosR parity gates are added, removed, or materially changed, update
  `docs/parity.md` and the science audit in the same change.

## What PhosR Parity Means Here

PhosR parity in this repository is:

- seam-level
- selective
- tied to committed fixtures
- strongest where lanes are explicitly classified as
  `PHOSR_PARITY_ACTIVE_SCIENCE`

PhosR parity here does not mean:

- every PhosR feature is implemented
- every Python path must numerically match PhosR
- every `PHOSR_ALIGNED` row has the same regression confidence tier

## Active PhosR Parity Science

The following areas currently run as active PhosR parity science in supported
PhosPy lanes (`PHOSR_PARITY_ACTIVE_SCIENCE`):

- prediction-science PhosR parity on committed prediction-science fixtures,
  including fragile-support reference tables
- kinase workflow PhosR parity on the supported L6 PhosPy lane
- adaptive prediction PhosR parity from promoted adaptive-sampling fixtures,
  executed in both supported PhosPy policy lanes:
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`, with
  policy-specific checks and stable-vs-`r_parity` divergence checks reported
  separately
- adaptive replay-trace PhosR parity from promoted replay fixtures, including
  initial negative-set replay surfaces, per-iteration sample-membership replay
  surfaces, final ensemble probabilities, top-k replay summaries, and
  stable-vs-`r_parity` comparison metrics under fixed seed
- public end-to-end `predMat` benchmark PhosR parity for both supported adaptive
  policies
- public end-to-end `predMat` order-invariance PhosR parity on the
  stable/default lane
- activity-stage PhosR parity outputs from fixed `predMat` + phospho inputs
- full-table signalome regression contracts on the supported L6 lane:
  `module_assignments`, `signalome_modules`, `kinase_network.nodes`,
  `kinase_network.edges`, `expanded_signalome`
- preprocessing-science PhosR parity on supported builder lanes:
  `total_protein_correction.policy="subtract_log_total"` with
  `intensity_transform.policy="log2"`,
  `site_matrix.policy="build_from_metadata"`, and
  `comparisons.policy="sample_metadata_pairs"`

## PhosPy Validated Coverage

No lanes are currently classified in this tier in the audited `2026-04-22`
snapshot.

The core kinase scoring/prediction lane (`SCI-GAP-12`) was promoted out of this
tier after ranking comparison-surface repair and restoration of ranking-focused
closure gates on explicit like-for-like surfaces.

## Contract-Changed Supported Lanes

These lanes are intentionally supported but not PhosR-equivalent by contract:

- adaptive sampling public contract uses `adaptive_policy` naming, not PhosR
  `svm_mode`
- core kinase scoring/prediction runtime defaults are intentionally reshaped:
  profile+motif combine keeps profile-only fallback, missing motif values
  preserve profile evidence, and workflow candidate filtering uses
  `score_threshold=0.0` with `inclusion=1`
- builder transformation establishment is intentionally narrow:
  `intensity_scale_state.label == "linear"` via identity establishment path
- signalome requires explicit `site_metadata.protein_id`
- signalome public entrypoint is contracted to
  `SignalomeWorkflowRequest(kinase_result=...)`
- motif sequence authority in the supported kinase lane is
  `references.site_sequences`, not dataset-sequence fallback

## Kinase Scoring/Prediction PhosPy-vs-PhosR Classification

Audit snapshot: `2026-04-22`

| Difference | PhosPy contract | PhosR baseline | Classification |
| --- | --- | --- | --- |
| Profile-only fallback in score combine | Enabled in supported workflow scoring path | PhosR default disabled | Intentional and supported |
| Missing motif-value handling | Profile score is preserved when motif value is missing for that kinase/site cell | No explicit profile-rescue path in PhosR combine | Intentional and supported |
| Workflow candidate filtering defaults | Candidate selection uses `score_threshold=0.0`, `inclusion=1` with caller-owned `top_k` | PhosR defaults used `score_threshold=0.8`, `inclusion=20`, `top=50` | Intentional and supported |
| Request/config knobs | PhosR-style knobs (`allow_profile_only_fallback`, `score_threshold`, `inclusion`, `min_motif_size`, `svm_mode`, `profile_policy`) are out of public contract | PhosR prediction config exposed those knobs | Intentional and supported |
| Ranking comparison surface | Ranking closure bars run on explicit like-for-like surfaces: `phospy(stable)` vs `promoted_reference(stable)` for prediction-matrix ranking and ranked top-k export ranking; cross-policy divergence (`stable` vs `r_parity`) is tracked separately | Earlier mixed source/policy reporting is not used for closure bars | Repaired and governance-locked |
| Ranking gate posture | Ranking gates are enforced as release bars in `tests/parity/test_l6_prediction_parity.py` using `tests/support/l6_prediction_parity_thresholds.py` | Earlier interim/loosened bars | Closure-grade gates active |

Resolved design decisions in this lane at this audit snapshot:

- Ticket 1 resolved: ranking comparison surfaces are now source-consistent and
  policy-consistent for closure decisions.
- Ticket 2 resolved: ranking-specific closure gates are restored on the
  repaired surface with explicit threshold configuration.

## Open Scientific Gaps

`OPEN_SCIENTIFIC_GAP` remains an active vocabulary tier and should be used
whenever evidence is incomplete or unresolved. In the current audited inventory
snapshot dated `2026-04-22`, no rows are classified `OPEN_GAP`.

This is not a blanket claim that every possible PhosR science surface is closed.
It applies only to the explicit inventory and boundaries documented here and in
the architecture audit.

## Core Kinase Lane Status

Audit snapshot: `2026-04-22`

The central kinase scoring/prediction lane is PhosR parity-gated for ranking
behaviour in this snapshot.

Active PhosR parity evidence still runs through:

- `tests/parity/test_l6_prediction_parity.py`
- `tests/parity/test_adaptive_replay_parity.py`

This evidence now gives closure-grade confidence for scoring numerics,
candidate selection, ranking behaviour, and replay behaviour in promoted PhosPy
fixture lanes.

This status does not override contract-changed defaults in this lane. Those
defaults are intentionally documented and tested as PhosPy contract.

## L6 Ranking Gate Policy

Audit snapshot: `2026-04-22`

For `tests/parity/test_l6_prediction_parity.py`, ranking thresholds are
release-governance gates configured in
`tests/support/l6_prediction_parity_thresholds.py`:

- prediction-matrix ranking gates:
  - mean Spearman rank correlation `>= 0.99`
  - mean top-20 overlap `>= 0.95`
  - mean top-30 overlap `>= 0.95`
  - kinases with top-10 overlap >= 70%: `>= 24`
- ranked top-k export ranking gates:
  - mean Spearman rank correlation `>= 0.99`
  - mean top-20 overlap `>= 0.95`
  - mean top-30 overlap `>= 0.95`
  - kinases with top-10 overlap >= 70%: `>= 24`

Current governance interpretation:

- these checks are closure-grade PhosR parity gates on repaired like-for-like
  surfaces
- cross-policy divergence (`stable` vs `r_parity`) remains separately reported
  and separately asserted, and is not used as a replacement for
  promoted-reference ranking closure bars

Informational diagnostics:

- `test_l6_prediction_parity_reporting_is_surface_explicit` is marked
  `parity_diagnostic` and records surface metrics for operator visibility
- this diagnostic reporting is separate from threshold assertions and does not
  itself define release bars

## PhosPy-Owned PhosR Parity Reporting

PhosR parity reporting is emitted by default from `tests/parity/`. Reporting is
PhosPy-owned through `tests/support/` and `tests/conftest.py`.

When PhosR parity tests run, terminal output includes grouped scientific
summaries for prediction science, kinase workflow, L6 core scoring/prediction,
adaptive prediction, adaptive replay trace, public end-to-end `predMat`,
`predMat` order invariance, activity stage, preprocessing science, and
signalome workflow coverage.

No `PHOSPY_SHOW_*` environment variables are required.

## Active PhosR Parity Gate Files

Audit snapshot: `2026-04-22`

Active PhosR parity evidence currently resolves through the following
PhosPy-owned files:

- `tests/parity/test_prediction_science_parity.py`
- `tests/parity/test_kinase_workflow_parity.py`
- `tests/parity/test_l6_prediction_parity.py`
- `tests/parity/test_adaptive_prediction_parity.py`
- `tests/parity/test_adaptive_replay_parity.py`
- `tests/parity/test_public_predmat_parity.py`
- `tests/parity/test_activity_stage_parity.py`
- `tests/parity/test_preprocessing_science_parity.py`
- `tests/parity/test_signalome_workflow_parity.py`

## PhosR Science Coverage Inventory

This table is the PhosR parity governance truth source for tracked science
areas. `Status` and `Coverage tier` are intentionally separate columns so that
`PHOSR_ALIGNED` does not imply closure-grade PhosR parity coverage.

| PhosR science area | Status | Coverage tier | Contract relation | PhosPy coverage summary |
| --- | --- | --- | --- | --- |
| profile policy behaviour | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | `strict` + `median_skipna` profile behaviour is supported and PhosR parity-tested. |
| core kinase scoring/prediction lane | CONTRACT_CHANGED | PHOSR_PARITY_ACTIVE_SCIENCE | Contract changed: supported defaults and fallback policy differ from PhosR | Numeric scoring/candidate/ranking/replay evidence is PhosR parity-gated on repaired like-for-like surfaces with explicit ranking threshold configuration and source/policy contract assertions. |
| adaptive sampling / svm_mode | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed: `adaptive_policy` replaces PhosR `svm_mode` naming | Adaptive science is implemented with PhosR parity evidence, but contract naming intentionally differs from PhosR. |
| signalome clustering/module selection | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | Clustering and module-count diagnostics are implemented and PhosR parity-backed. |
| weighted-top assignment behaviour | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | Weighted-top assignment and fractional support propagation are implemented and PhosR parity-backed. |
| network policy variants | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | `positive_only`, `absolute_threshold`, and `signed` are implemented and PhosR parity-tested. |
| expanded signalome outputs | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | `expanded_signalome` is materialised in the supported workflow path and PhosR parity-tested. |
| activity comparison lock | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | Activity thresholded-substrate-mean science is implemented and guarded by CI gates. |
| preprocessing transformation establishment | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed: narrow builder establishment policy | Supported builder lane establishes only `linear` pass-through transformation state. |
| total/protein correction | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | `total_protein_correction.policy="subtract_log_total"` is PhosR parity-gated in PhosPy-owned fixture tests, with explicit log-scale formula `log2(phospho + pseudocount) - log2(total + pseudocount)` and strict phospho/total alignment behaviour. |
| site-matrix construction | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | Supported site-matrix construction (`build_from_metadata`) is PhosR parity-gated with PhosPy-owned fixture expectations for row retention, site identity, and output matrix values. |
| comparison-building | PHOSR_ALIGNED | PHOSR_PARITY_ACTIVE_SCIENCE | PhosR-aligned in supported lane | Sample-metadata comparison construction is PhosR parity-gated for explicit and inferred pair policies, including pair identity/order and expected output values. |
| site-to-protein resolution fallback behaviour | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed: no PhosR-style site-id-prefix fallback | Signalome requires explicit `site_metadata.protein_id`. |
| signalome input route contraction | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed: workflow entrypoint narrowed | Supported signalome entrypoint is `SignalomeWorkflowRequest(kinase_result=...)`. |
| dataset-vs-reference sequence authority decisions | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed: reference bundle is sequence authority | Motif sequence authority in supported kinase lane is `references.site_sequences`. |

Open PhosR-science areas in this inventory snapshot:

- none currently classified `OPEN_GAP` in the audited list.

Visibility check:

- inventory rows and coverage-tier claims are validated directly in
  `docs/parity.md` and the science audit.

## Fixture Locations

### PhosR parity inputs and expectations

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

These files are the normal source for active PhosR parity tests in
`tests/parity/` and helpers in `tests/support/rewrite_fixture_data.py`.

### PhosPy workflow regression expectations

- public `predMat` benchmark inputs and committed PhosPy outputs:
  `tests/fixtures/public_workflow_reference/predmat_input_*.{csv,json}`,
  `tests/fixtures/public_workflow_reference/predmat_rewrite_*.csv`,
  `tests/fixtures/public_workflow_reference/predmat_rewrite_contract.json`
- provenance and promotion history:
  `tests/fixtures/public_workflow_reference/PROVENANCE.md`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*.csv`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json`

### Historical Reference Archive

- Historical PhosR parity and debug artifacts live in repository history.
- Active PhosR parity tests resolve only from committed fixture roots under
  `tests/fixtures/`.

## Run the PhosR Parity Suite

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

## Where Next

- Inventory and evidence detail:
  [Scientific coverage audit](architecture/science_gap_audit.md)
- High-level project status:
  [Roadmap](roadmap.md)
- Contract behaviour context:
  [API Guide](api.md), [Validation Guide](validation.md)

