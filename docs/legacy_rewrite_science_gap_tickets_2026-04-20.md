# Legacy Science Gap Tickets for PhosPy Rewrite

Date: 2026-04-20  
Source review: `legacy_rewrite_science_gap_review_2026-04-20.md`

These tickets turn the confirmed review findings into implementation-ready work items. They are written from the perspective of **matching the legacy scientific support story where that story is still intended**, and making the public contract truthful where support is intentionally narrowed.

---

## Ticket 1 — Replace the narrow donor inventory with a truthful legacy-science coverage inventory

**Priority:** P1  
**Type:** Governance / Scientific Coverage / Audit Truthfulness  
**Area:**

- `docs/architecture/legacy_science_gap_audit.md`
- `docs/parity.md`
- `tests/support/legacy_donor_inventory.py`
- parity governance docs and any release-facing parity statements

### Summary

Replace the current scoped donor inventory and “no open gaps” wording with a legacy-science inventory that truthfully represents what has and has not been ported.

### Problem

The current audit/governance story presents the rewrite as having no open legacy science gaps, but the tracked donor inventory only covers a limited set of kinase/activity/signalome surfaces. It omits major legacy scientific domains including preprocessing, total/protein correction, site-matrix construction, comparison-building, site-to-protein fallback behavior, and dataset-vs-reference sequence authority.

That makes the audit look broader than the supporting evidence actually is.

### Why this matters

This is the core trust problem. As long as the inventory is selective, downstream “parity complete” language is misleading even if the currently tracked tickets are genuinely closed.

### Goals

- Make the inventory describe the **whole legacy science surface**, not just the currently promoted donor set.
- Distinguish clearly between:
  - fully ported legacy scientific areas,
  - intentionally retired areas,
  - open parity gaps,
  - contract changes that narrow supported behavior.
- Prevent future “no open gaps” claims unless the inventory actually supports them.

### Scope

**In scope**

- donor inventory model and docs;
- classification of all currently confirmed omissions from the review;
- explicit status for preprocessing, site-to-protein resolution, and sequence-authority decisions;
- tests or lint-style checks that keep inventory/docs consistent.

**Out of scope**

- porting the missing science itself, unless separately ticketed.

### Required changes

- Expand the donor inventory so it includes the omitted legacy science domains found in the review.
- Reclassify the current status language so “scoped pass” and “full legacy science parity” are not conflated.
- Add an explicit status field per donor area, for example:
  - `PORTED`
  - `INTENTIONALLY_RETIRED`
  - `OPEN_GAP`
  - `CONTRACT_CHANGED`
- Ensure docs and inventory are derived from the same source of truth, or at minimum checked for consistency in tests.
- Remove or narrow any blanket “no open gaps” wording until the open items are resolved.

### Acceptance criteria

- The inventory includes preprocessing, total/protein correction, site-matrix construction, comparison-building, site-to-protein resolution, signalome input route contraction, and sequence-authority decisions.
- The docs no longer imply full legacy-science parity when only a subset is covered.
- There is a failing test or check if docs and inventory drift on open/closed status.
- A maintainer can read the parity docs and understand exactly which legacy scientific areas are still open.

---

## Ticket 2 — Restore or explicitly retire the legacy preprocessing science lane in the public dataset-builder story

**Priority:** P1  
**Type:** Scientific Correctness / Public Contract / Preprocessing  
**Area:**

- `src/phospy/api/requests.py`
- `src/phospy/api/configs.py`
- `src/phospy/datasets/builders/*`
- `src/phospy/validation/datasets/*`
- `docs/api.md`
- `docs/validation.md`
- README dataset construction guidance
- parity and builder tests

### Summary

Decide whether PhosPy still intends to support the legacy preprocessing science lane for turning raw phosphoproteomics inputs into analysis-ready datasets. If yes, port it as a supported builder story with user-visible configuration. If no, explicitly retire it and narrow the public contract truthfully.

### Problem

The legacy system had a real preprocessing layer, including localization filtering, coverage filtering, sentinel replacement, minimum-observed filtering, duplicate-gene handling, phospho-to-protein correction, pairwise augmentation, and site-matrix construction controls. The rewrite builder currently accepts already-separated tables and exposes only a narrow missing-data rule plus minimum-observed filtering.

That is not a like-for-like scientific port of the legacy input-conditioning story.

### Why this matters

This is a product-boundary science decision, not just an implementation detail. If PhosPy claims to be a maintainable port of PhosR, users need a truthful answer to whether the package supports getting from raw-ish quantitative inputs to analysis-ready matrices, or whether that work is now out of scope.

### Goals

- Choose the supported preprocessing contract deliberately.
- Keep preprocessing choices user-visible rather than hidden inside implementation defaults.
- Match legacy scientific behavior where that behavior remains in scope.
- Keep the dataset boundary and docs truthful if the supported contract is intentionally narrower.

### Scope

**In scope**

- request/config surface for preprocessing;
- builder orchestration;
- docs and examples;
- tests for supported preprocessing behaviors.

**Out of scope**

- specific sub-areas that warrant their own tickets, such as total/protein correction, site-matrix policy details, and comparison-building, except where needed to define the boundary.

### Required changes

- Decide between two explicit product directions:
  - **Option A — Supported preprocessing lane:** port the intended legacy preprocessing stages and expose them through a user-facing request/config surface.
  - **Option B — Explicitly analysis-ready-only boundary:** remove any implication that the builder is a legacy preprocessing equivalent and document exactly what the user must provide before entering PhosPy.
- If Option A is chosen:
  - introduce a preprocessing request/config model that can carry the supported scientific controls;
  - port only the intended legacy stages and document unsupported ones clearly;
  - keep configuration user-provided rather than silently inferred.
- If Option B is chosen:
  - simplify the builder story around normalization/validation only;
  - remove wording that implies legacy preprocessing equivalence;
  - state clearly that raw-table-to-analysis-ready conditioning is outside the supported product.

### Acceptance criteria

- The public dataset-builder docs clearly state whether PhosPy supports legacy-style preprocessing or only accepts already-analysis-ready inputs.
- The request/config surface matches that decision and does not hide scientific policy behind internal defaults.
- Tests lock the chosen contract so new documentation drift cannot reintroduce ambiguity.
- A user can tell, before running the package, whether PhosPy will perform legacy-style scientific preprocessing for them.

---

## Ticket 3 — Restore phospho-to-protein correction as an active supported scientific stage, or remove `total` from the active workflow story

**Priority:** P1  
**Type:** Scientific Correctness / Preprocessing / Workflow Contract  
**Area:**

- `src/phospy/datasets/models.py`
- `src/phospy/datasets/builders/*`
- `src/phospy/workflows/*`
- `src/phospy/api/requests.py`
- `src/phospy/api/configs.py`
- `docs/api.md`
- `docs/validation.md`
- README and examples
- preprocessing/workflow parity tests

### Summary

Decide whether the supported product treats the total/protein matrix as a scientifically active input. If yes, port the legacy phospho-to-protein correction path and ensure downstream workflows operate on the corrected phospho values. If no, stop presenting `total` as part of the active scientific workflow story.

### Problem

The rewrite accepts and validates `total`, stores it in the dataset, and carries it through bundle/publishing paths, but there is no evidence that workflow science consumes `dataset.total`. In legacy, total/protein values were not passive attachments; they participated in phospho-to-protein correction and related preprocessing outputs.

### Why this matters

This is a material legacy science loss. If users provide total/protein abundance expecting the legacy scientific effect, the current rewrite gives them a misleadingly rich dataset object with no comparable downstream behavior.

### Goals

- Make the role of `total` scientifically real or explicitly non-scientific.
- Match legacy correction behavior if the feature remains supported.
- Prevent passive metadata from being mistaken for active workflow input.

### Scope

**In scope**

- correction-stage design;
- relationship between preprocessing output and workflow input;
- docs and examples;
- parity tests against legacy correction behavior.

**Out of scope**

- unrelated changes to kinase or signalome algorithms except where corrected phospho values affect their input matrices.

### Required changes

- Choose one supported policy:
  - **Option A — Active total-correction lane:** port legacy phospho-to-protein correction and make it part of the supported builder/preprocessing story.
  - **Option B — Passive total metadata only:** keep `total` attachable but document clearly that it is stored, validated, and published only, not used in workflow science.
- If Option A is chosen:
  - identify the authoritative legacy correction behavior;
  - implement it in the rewrite preprocessing/builder lane;
  - ensure downstream workflow input uses corrected phospho values where appropriate;
  - add fixture-backed tests that prove corrected outputs match legacy behavior.
- If Option B is chosen:
  - remove or narrow wording that suggests `total` participates in analysis;
  - consider whether the dataset model should still accept `total` in the main supported lane.

### Acceptance criteria

- A user can tell from docs and API shape whether supplying `total` changes scientific results.
- If correction is supported, there are parity tests showing the corrected phospho matrix matches legacy outputs for representative fixtures.
- If correction is not supported, docs and examples no longer imply otherwise.
- There is no ambiguity between “stored on the dataset” and “used by workflow science”.

---

## Ticket 4 — Restore legacy site-matrix construction and duplicate-site policy controls, or explicitly keep them outside the product boundary

**Priority:** P1  
**Type:** Scientific Correctness / Preprocessing / Public Contract  
**Area:**

- `src/phospy/datasets/builders/*`
- `src/phospy/api/requests.py`
- `src/phospy/api/configs.py`
- `src/phospy/validation/datasets/*`
- `docs/api.md`
- `docs/validation.md`
- matrix/preprocessing tests

### Summary

Decide whether PhosPy supports the legacy site-matrix construction lane, including duplicate-site strategy and missing-data policy choices. If yes, port it as a user-visible builder capability. If no, document that canonical site matrices must be produced outside PhosPy.

### Problem

Legacy preprocessing could build canonical site matrices from raw identifiers and apply explicit policy controls for duplicate handling and missing-data treatment. The rewrite expects canonical aligned matrices to already exist and offers no equivalent site-matrix construction route.

### Why this matters

This is one of the largest dropped scientific input-conditioning capabilities. It moves a major part of the legacy scientific workflow outside the package without a clearly declared product decision.

### Goals

- Decide whether site-matrix construction is in or out of scope.
- If in scope, expose duplicate-site and missing-data policy as user-visible choices.
- If out of scope, make the analysis-ready-only boundary explicit and unambiguous.

### Scope

**In scope**

- canonical site ID construction policy;
- duplicate-site handling;
- missing-data policy at site-matrix construction time;
- docs and tests.

**Out of scope**

- downstream kinase/signalome science except as consumers of the resulting site matrix.

### Required changes

- Choose one supported direction:
  - **Option A — Port the site-matrix builder:** support the intended legacy matrix-construction behaviors and policy controls.
  - **Option B — Explicitly externalize site-matrix construction:** require users to provide canonical site matrices and remove all implication that PhosPy performs this step.
- If Option A is chosen:
  - define the supported duplicate-site strategies;
  - define the supported missing-data policies at this stage;
  - implement them through a stable request/config surface;
  - add regression tests based on legacy matrix fixtures.
- If Option B is chosen:
  - tighten docs to state that duplicate resolution and raw identifier collapse are upstream responsibilities.

### Acceptance criteria

- The package either supports site-matrix construction with documented policy controls or clearly does not.
- Duplicate-site behavior is not left to ad hoc caller preprocessing.
- Tests lock the chosen contract.
- The docs no longer blur “analysis-ready dataset” with “raw site matrix construction”.

---

## Ticket 5 — Restore legacy comparison-building support, or explicitly retire grouped/comparison derivation from the builder contract

**Priority:** P2  
**Type:** Scientific Capability / Builder Contract / Preprocessing  
**Area:**

- `src/phospy/api/requests.py`
- `src/phospy/api/configs.py`
- `src/phospy/datasets/builders/*`
- `docs/api.md`
- README examples
- comparison/preprocessing tests

### Summary

Decide whether the supported builder story includes legacy comparison-building from grouped/raw quantitative inputs. If yes, port the relevant `schema`/`comparisons` capabilities. If no, retire that route explicitly in the public API and docs.

### Problem

Legacy analysis-ready construction supported `schema` and `comparisons` inputs and could produce pairwise comparison-oriented outputs. The rewrite builder request does not include this surface and instead assumes already-shaped inputs.

### Why this matters

This is part of the missing raw-to-analysis-ready scientific path. It narrows the package from a workflow-building system into a consumer of already-prepared matrices, which may be acceptable, but only if declared plainly.

### Goals

- Decide whether comparison-building remains a supported scientific capability.
- Preserve truthful API/docs around grouped input support.
- Match legacy comparison derivation where the feature remains in scope.

### Scope

**In scope**

- request/config surface;
- builder support for grouped/comparison derivation if retained;
- docs and regression tests.

**Out of scope**

- unrelated workflow algorithm changes.

### Required changes

- Choose one direction:
  - **Option A — Support comparison-building:** port the relevant legacy grouped/comparison derivation behavior.
  - **Option B — Retire comparison-building from the product boundary:** keep builder inputs analysis-ready only and document that comparison derivation is out of scope.
- If Option A is chosen:
  - define the supported schema/comparison request model;
  - implement the intended derivation logic;
  - add tests covering representative grouped-input cases.
- If Option B is chosen:
  - remove any lingering language that suggests the builder can derive comparison matrices from grouped/raw data.

### Acceptance criteria

- The public builder contract clearly states whether grouped/comparison derivation is supported.
- The API shape matches that statement.
- Tests cover the supported direction and prevent ambiguous future expansion.

---

## Ticket 6 — Restore a truthful and scientifically supported site-to-protein resolution contract for SignalomeWorkflow

**Priority:** P1  
**Type:** Scientific Correctness / Validation / Public Contract  
**Area:**

- `src/phospy/workflows/signalome/validator.py`
- `src/phospy/workflows/signalome/interpreter.py`
- `src/phospy/datasets/models.py`
- `src/phospy/api/requests.py`
- `docs/api.md`
- `docs/validation.md`
- signalome parity tests

### Summary

Decide and enforce the supported site-to-protein resolution contract for signalome execution. Either keep the current strict explicit `protein_id` boundary and document the legacy support contraction honestly, or restore a supported user-controlled resolution system with diagnostics equivalent to the legacy behavior.

### Problem

Legacy signalome execution supported strict and metadata-fallback site-to-protein resolution, optional gene-symbol fallback, ambiguity handling, and structured diagnostics. The rewrite removes that machinery and now requires `dataset.site_metadata.protein_id` to already be explicit and non-blank.

### Why this matters

Protein identity is not a cosmetic label in signalome grouping. This is a real scientific contract change. It may be the right contract, but it is not a parity-complete one.

### Goals

- Make the supported signalome protein-identity contract explicit and scientifically honest.
- Preserve user control over interpretation policy if legacy-style fallback remains in scope.
- Avoid silent or lossy protein-identity inference.

### Scope

**In scope**

- signalome input validation;
- dataset/model support for resolution helpers if retained;
- diagnostics and error reporting;
- docs and parity tests.

**Out of scope**

- unrelated signalome downstream clustering/network logic.

### Required changes

- Choose one supported policy:
  - **Option A — Strict explicit protein identity only:** require trusted `protein_id` and document that legacy fallback routes are intentionally not supported.
  - **Option B — Supported controlled fallback:** restore site-to-protein resolution modes, candidate-column selection, optional gene-symbol fallback, ambiguity handling, and structured diagnostics.
- If Option B is chosen:
  - ensure fallback behavior is user-controlled rather than silently inferred;
  - preserve diagnostics so ambiguous or lossy mappings are inspectable;
  - add regression tests matching representative legacy resolution cases.
- If Option A is chosen:
  - tighten docs so users know fallback behavior is gone;
  - ensure errors clearly tell users what metadata must be provided.

### Acceptance criteria

- The signalome API/docs make the protein-identity contract unambiguous.
- If fallback is supported, diagnostics and policy controls are exposed and tested.
- If fallback is not supported, the docs explicitly state the legacy contraction and the validator enforces it consistently.
- There is no silent conversion of gene-like labels into claimed protein identity.

---

## Ticket 7 — Restore explicit signalome input routes for externally provided site-to-protein interpretation, or document their retirement

**Priority:** P2  
**Type:** Scientific Workflow Contract / Public API / Signalome  
**Area:**

- `src/phospy/api/requests.py`
- `src/phospy/workflows/signalome/interpreter.py`
- `docs/api.md`
- signalome API and integration tests

### Summary

Decide whether SignalomeWorkflow should support the legacy-style input routes that allowed explicit external site-to-protein interpretation choices, rather than deriving protein grouping only from embedded dataset metadata.

### Problem

Legacy signalome execution supported running from analysis-ready datasets with explicit site-to-protein mapping choices. The rewrite public request accepts a `KinaseWorkflowResult` plus config and derives site-to-protein grouping only from `dataset.site_metadata.protein_id`.

### Why this matters

This removes a user-controlled interpretation lane that existed in the legacy scientific surface. Even if strict embedded metadata remains the default story, the loss of explicit external mapping input is still a support contraction.

### Goals

- Decide whether explicit external mapping is part of the supported signalome surface.
- Keep the public API aligned with the intended scientific flexibility.
- Avoid a hidden assumption that all protein interpretation must already be baked into the dataset.

### Scope

**In scope**

- signalome request shape;
- interpreter input handling;
- docs and integration tests.

**Out of scope**

- site-to-protein fallback mechanics themselves, except where needed to support explicit supplied mappings.

### Required changes

- Choose one direction:
  - **Option A — Restore explicit mapping input support:** allow the signalome request surface to accept supported site-to-protein interpretation input directly.
  - **Option B — Keep metadata-only grouping:** document that external mapping control is no longer supported and all grouping inputs must already exist in dataset metadata.
- If Option A is chosen:
  - define the supported request model for explicit mapping input;
  - ensure precedence between explicit mappings and dataset metadata is well defined;
  - add integration tests for both embedded and explicit-mapping routes.
- If Option B is chosen:
  - state the retirement clearly in docs and migration guidance.

### Acceptance criteria

- The signalome public API truthfully states whether externally supplied mapping interpretation is supported.
- The request/interpreter surface matches that statement.
- Integration tests cover the chosen supported lane.

---

## Ticket 8 — Decide and enforce the authoritative sequence source for kinase motif scoring

**Priority:** P1  
**Type:** Scientific Correctness / Public Contract / Kinase Workflow  
**Area:**

- `src/phospy/workflows/kinase/interpreter.py`
- `src/phospy/datasets/models.py`
- `src/phospy/datasets/builders/*`
- reference bundle loading code
- `docs/api.md`
- `docs/validation.md`
- README and kinase parity tests

### Summary

Decide which sequence source is authoritative for motif scoring in the supported kinase workflow: dataset-provided `site_sequence`, reference-bundle `site_sequences`, or an explicitly validated precedence rule. Then enforce that policy consistently in code, validation, docs, and tests.

### Problem

Legacy kinase execution consumed sequence information owned by the analysis-ready dataset. The rewrite kinase workflow uses `references.site_sequences` as the supported scoring input. Dataset-level `site_sequence` is still present in the model and builder story, but it is not authoritative for motif scoring.

### Why this matters

This is a scientific authority change, not a cosmetic refactor. It affects whether user-provided site sequences can drive motif scoring, whether custom datasets can be scored without matching reference bundle sequence entries, and what the dataset contract actually means.

### Goals

- Choose a single truthful authority model for motif-scoring sequences.
- Avoid carrying scientifically important sequence data on the dataset while ignoring it in the workflow.
- Match legacy behavior if dataset-owned sequence remains the intended scientific source.

### Scope

**In scope**

- kinase interpreter sequence sourcing;
- dataset/reference validation around sequence availability;
- docs and tests.

**Out of scope**

- unrelated changes to downstream kinase ranking logic.

### Required changes

- Choose one supported policy:
  - **Option A — Dataset sequence authoritative:** motif scoring uses `dataset.site_metadata.site_sequence` and validates it as required scientific input.
  - **Option B — Reference bundle authoritative:** motif scoring uses reference bundle sequences and dataset `site_sequence` is non-authoritative metadata.
  - **Option C — Explicit precedence contract:** define when dataset sequence overrides reference sequence, when mismatch is an error, and when one may be absent.
- Implement the chosen policy consistently in the interpreter and validators.
- Update docs so the authoritative source is obvious.
- Add regression tests for the chosen authority rule, including mismatch and missing-data cases where relevant.

### Acceptance criteria

- The supported kinase workflow has a single documented authoritative sequence source or precedence contract.
- Code, validation, and docs all reflect the same rule.
- Tests cover success and failure cases for the chosen policy.
- A user can tell whether supplying `dataset.site_metadata.site_sequence` will affect motif scoring.

---

## Suggested execution order

1. Ticket 1 — fix the governance/trust story first.  
2. Ticket 2 — decide the preprocessing product boundary.  
3. Ticket 3 — resolve the role of `total`.  
4. Ticket 4 — resolve site-matrix construction support.  
5. Ticket 6 — resolve signalome protein-identity contract.  
6. Ticket 8 — resolve motif-sequence authority.  
7. Ticket 5 — comparison-building, depending on how broad the preprocessing lane remains.  
8. Ticket 7 — explicit signalome input route support, depending on the outcome of Ticket 6.
