# PhosPy Tickets: Preprocessing Structure First, Then Missing Legacy Science

Date: 2026-04-20
Repository snapshot reviewed: `app-src77.zip`

This ticket set assumes the current rewrite already has the right **outer builder shell**:

- `src/phospy/datasets/builders/public.py`
- `src/phospy/datasets/builders/validator.py`
- `src/phospy/datasets/builders/interpreter.py`
- `src/phospy/datasets/builders/executor.py`
- `src/phospy/datasets/builders/normalizer.py`
- `src/phospy/datasets/builders/reader.py`
- `src/phospy/datasets/builders/sequence_derivation.py`

The missing work is primarily inside the dataset-builder preprocessing boundary.

The current confirmed open legacy-science areas are:

- total/protein correction
- site-matrix construction
- comparison-building

The current governance docs also need to be corrected because they still claim there are no open scientific gaps in scope.

---

## Ticket 1 — Re-open the preprocessing science inventory and make the governance story truthful

**Priority:** P1  
**Type:** Governance / Scientific Correctness / Public Contract  
**Area:**

- `docs/architecture/legacy_science_gap_audit.md`
- `tests/support/legacy_donor_inventory.py`
- `docs/parity.md`
- `docs/roadmap.md`
- any linked architecture notes that currently imply preprocessing parity is complete

### Summary

Update the project’s parity and governance truth sources so they explicitly show that preprocessing science is still open in the rewrite lane. The audit should no longer imply that legacy-science coverage is complete when total/protein correction, site-matrix construction, and comparison-building are still missing from the supported builder path.

### Problem

The current audit states that there are no confirmed open scientific gaps in scope. That is not consistent with the current codebase or with the legacy donor surface still present in `legacy_archive/phospy_legacy/preprocessing/*`.

This creates two risks:

1. maintainers can make planning decisions from an inaccurate parity picture;
2. users and reviewers can over-trust the rewrite builder story.

### Why this matters

If the governance story is ahead of the implementation, all downstream ticketing and release confidence become unreliable. Before adding new science, the project needs a truthful statement of what remains open.

### Goals

- Mark preprocessing science as explicitly open.
- Add tracked open-ticket identifiers for the three confirmed preprocessing gaps.
- Extend the donor inventory so preprocessing donor areas are represented alongside the currently tracked kinase/signalome donors.
- Make docs and parity notes say clearly that the open work is in the dataset-builder preprocessing lane, not in workflow interpreters.

### In scope

- audit and donor inventory updates
- roadmap/parity wording updates
- creation of new preprocessing gap ticket identifiers

### Out of scope

- implementing preprocessing science itself
- changing the public builder API in this ticket

### Required changes

- Replace the “no open scientific gaps” claim with an explicit list of open preprocessing gaps.
- Add donor inventory entries for:
  - total/protein correction
  - site-matrix construction
  - comparison-building
- Link each donor area to one concrete open ticket.
- Update wording so the project distinguishes:
  - landed workflow science
  - still-open dataset-builder preprocessing science

### Acceptance criteria

- `legacy_science_gap_audit.md` no longer claims the open-gap count is zero.
- `tests/support/legacy_donor_inventory.py` includes the three preprocessing donor areas.
- parity/governance docs describe the open work as builder-boundary science.
- ticket references are stable enough to be cited from future reviews.

### Suggested ticket IDs

- `SCI-GAP-13` — preprocessing architecture
- `SCI-GAP-14` — total/protein correction
- `SCI-GAP-15` — site-matrix construction
- `SCI-GAP-16` — comparison-building

---

## Ticket 2 — Establish a real internal preprocessing subsystem under the existing dataset-builder shell

**Priority:** P1  
**Type:** Architecture / Scientific Correctness / Maintainability  
**Area:**

- `src/phospy/datasets/builders/executor.py`
- `src/phospy/datasets/builders/contracts.py`
- `src/phospy/datasets/builders/preprocessing.py`
- `src/phospy/api/configs.py`
- `src/phospy/validation/datasets/preprocessing.py`
- new internal preprocessing modules under `src/phospy/datasets/preprocessing/` or an equally clear dataset-owned package
- builder unit and integration tests

### Summary

Replace the current thin `DatasetPreprocessor` missing-value helper with a staged internal preprocessing subsystem that sits inside the existing builder executor path. Keep the public builder story unchanged, but give preprocessing its own internal state, plan, and stage boundaries so the missing legacy science can land without turning `builders/preprocessing.py` into a god module.

### Problem

The current builder architecture is good at the outer orchestration level, but the actual preprocessing lane is still a single narrow helper centered on `missing_data_policy`. That is too small for the remaining science and creates a predictable failure mode: all future preprocessing logic gets shoved into one builder-side file.

### Why this matters

The project already has the right shell. The risk now is internal collapse of responsibilities:

- executor becomes orchestration plus science plus policy routing;
- one preprocessing file becomes a dumping ground;
- restoration of missing science becomes harder to test and reason about.

This ticket exists to establish the structure before science restoration begins.

### Goals

- Keep one public builder entry point.
- Keep workflow interpreters downstream of the dataset boundary.
- Introduce a dataset-owned preprocessing subsystem with narrow responsibilities.
- Make the next three science tickets land as small stage additions rather than ad hoc executor growth.

### In scope

- internal preprocessing package/module layout
- preprocessing stage protocol and orchestration
- internal preprocessing state and plan objects
- executor integration with the new subsystem
- migration of existing missing-data behavior into the staged path

### Out of scope

- full implementation of total/protein correction
- full implementation of site-matrix construction
- full implementation of comparison-building
- widening the public builder return type beyond `AnalysisReadyPhosphoDataset`

### Required design direction

Keep the current public shell:

- `AnalysisReadyDatasetBuilder.run(request)`
- validator -> interpreter -> executor

Introduce internal preprocessing responsibilities such as:

- `PreprocessingPlan`
- `PreprocessingState`
- `PreprocessingPipeline`
- `PreprocessingStage`

Suggested package direction:

```text
src/phospy/datasets/preprocessing/
    models.py
    pipeline.py
    stages/
        missing_data.py
        total_protein_correction.py
        site_matrix.py
        comparisons.py
```

The exact filenames may be adjusted to fit the current codebase, but the responsibility split should remain.

### Required changes

- Refactor `builders/preprocessing.py` so it stops being the long-term home of all preprocessing science.
- Introduce an internal preprocessing state object that carries:
  - phospho matrix
  - site metadata
  - optional total matrix
  - optional sample metadata
  - preprocessing config / interpreted plan
- Introduce a preprocessing pipeline that applies ordered stages.
- Move current missing-data handling into its own stage.
- Keep `DatasetBuildExecutor` as the orchestration boundary rather than the owner of all preprocessing details.
- Keep preprocessing below the strict dataset boundary and above final dataset construction.

### Acceptance criteria

- the builder still exposes one public entry point and returns `AnalysisReadyPhosphoDataset`;
- the executor delegates preprocessing to a dedicated internal subsystem;
- missing-data behavior still works through the new staged path;
- there is a clear place to add the next three missing-science stages without bloating executor or one helper file;
- tests prove that workflow interpreters do not absorb preprocessing responsibilities.

### Test expectations

Add or update tests to cover:

- preprocessing pipeline ordering;
- stage isolation and state passing;
- executor delegation to the preprocessing subsystem;
- regression coverage for current `forbid` and `impute_row_median` behavior.

---

## Ticket 3 — Expand the public dataset preprocessing config into a builder-owned scientific policy surface

**Priority:** P1  
**Type:** Public Contract / Architecture / Validation  
**Area:**

- `src/phospy/api/configs.py`
- `src/phospy/api/requests.py`
- `src/phospy/validation/datasets/preprocessing.py`
- `docs/api.md`
- `docs/validation.md`
- dataset-builder tests

### Summary

Replace the current single-purpose `DatasetPreprocessingConfig` with a grouped builder-owned preprocessing config that can truthfully express the supported dataset-conditioning policies. Keep the configuration public and user-visible, rather than burying major scientific choices inside the builder interpreter or executor.

### Problem

The current config only models missing-data handling plus `min_observed_values`. That is not enough surface area for the preprocessing science the package still intends to support.

If total/protein correction, site-matrix construction, and comparison-building are going to be supported, those choices need a public home. Otherwise the package either hard-codes important science silently or pushes policy selection into internal components where the user cannot control it.

### Why this matters

The builder is the product boundary for preprocessing science. If a scientific policy is supported there, it should be visible in the builder request/config shape.

### Goals

- keep preprocessing configuration at the dataset-builder boundary;
- make major supported preprocessing policies user-visible;
- validate those policies centrally;
- avoid a giant flat config full of unrelated booleans.

### In scope

- redesign of `DatasetPreprocessingConfig`
- grouped config sub-objects or equivalent structured fields
- validation updates
- docs and examples

### Out of scope

- full implementation of all policies in this ticket
- changes to workflow request objects

### Required design direction

The config should stay builder-owned and should not move into kinase or signalome requests.

A healthy direction is grouped options such as:

- missing-data options
- total/protein correction options
- site-matrix options
- comparison-building options

The exact field names can be refined during implementation, but the public shape should support the missing science without becoming a bag of hidden defaults.

### Required changes

- Extend `DatasetPreprocessingConfig` beyond missing-data handling.
- Add validation rules for any new grouped options.
- Ensure builder interpretation resolves config once and passes an internal plan downstream.
- Update docs so users can see which preprocessing policies are supported, required, or unsupported.

### Acceptance criteria

- supported preprocessing policies have a clear public config home;
- validation errors are raised for unsupported or contradictory config combinations;
- workflow configs remain free of builder preprocessing options;
- docs/examples show preprocessing configuration at the dataset-builder boundary.

### Notes

This ticket is intentionally structural. It should land before the individual missing-science tickets so those tickets can attach to a stable config surface.

---

## Ticket 4 — Restore legacy total/protein correction in the supported dataset-builder preprocessing lane

**Priority:** P1  
**Type:** Scientific Correctness / Preprocessing / Parity  
**Area:**

- new preprocessing stage/module for total/protein correction
- `src/phospy/datasets/builders/executor.py`
- `src/phospy/api/configs.py`
- `src/phospy/validation/datasets/preprocessing.py`
- `src/phospy/datasets/models.py`
- `docs/api.md`
- `docs/validation.md`
- unit, integration, and parity-style donor tests
- legacy donor references in `legacy_archive/phospy_legacy/preprocessing/protein_correction.py`, `services.py`, and `steps.py`

### Summary

Port the legacy phospho-to-protein correction science into the supported rewrite builder lane as an explicit preprocessing stage. The dataset builder should be able to use total/protein measurements to generate corrected phosphosite values when that policy is requested and the required inputs are present.

### Problem

The current rewrite builder accepts an optional `total` matrix and carries it into the final dataset, but it does not perform the legacy correction step that used total/protein measurements to adjust phosphosite intensities before downstream analysis.

That means the rewrite currently supports transporting total data, but not using it for the legacy scientific purpose that justified it in preprocessing.

### Why this matters

This is one of the clearest remaining scientific holes in the port. If users provide total/protein data, the package needs a truthful answer to whether that data participates in supported preprocessing science or is merely preserved as auxiliary input.

### Goals

- restore total/protein correction as builder-boundary science;
- make correction opt-in or otherwise explicitly policy-driven;
- validate required inputs and failure modes clearly;
- align rewrite behavior with the legacy donor lane closely enough to support parity locking.

### In scope

- correction-stage implementation
- public config surface for enabling/configuring correction
- builder validation and error messages
- dataset assembly changes required to carry corrected outputs truthfully
- donor/regression tests

### Out of scope

- inventing a new correction algorithm unrelated to the legacy donor
- moving correction into workflow execution

### Required changes

- Introduce a dedicated preprocessing stage for total/protein correction.
- Define the supported public policy for when correction runs.
- Validate that required total/protein inputs are present before correction is attempted.
- Decide and document how corrected quantitative values are represented in the supported dataset boundary.
- Keep correction inside the dataset-builder preprocessing path, before workflow execution.

### Acceptance criteria

- builder preprocessing can execute a supported total/protein correction path;
- clear errors are raised when correction is requested but required inputs are absent or incompatible;
- docs explain whether correction is disabled by default, enabled by request, or required for some supported input story;
- tests lock the rewrite stage against legacy donor behavior for representative fixtures.

### Test expectations

Add tests covering:

- successful correction with supported phospho + total inputs;
- mismatch and unmatched-protein failure paths;
- behavior when correction is disabled;
- regression/parity comparisons against promoted donor fixtures or donor-derived expectations.

---

## Ticket 5 — Restore legacy site-matrix construction as a dedicated preprocessing stage below the analysis-ready dataset boundary

**Priority:** P1  
**Type:** Scientific Correctness / Preprocessing / Architecture  
**Area:**

- new preprocessing stage/module for site-matrix construction
- `src/phospy/api/configs.py`
- `src/phospy/validation/datasets/preprocessing.py`
- `src/phospy/datasets/builders/executor.py`
- `src/phospy/datasets/models.py`
- docs/tests around dataset-build behavior
- legacy donor references in `legacy_archive/phospy_legacy/preprocessing/site_matrix.py` and `core.py`

### Summary

Restore the legacy site-matrix construction behavior as a dedicated preprocessing stage in the rewrite builder path, with explicit policy control and a truthful relationship to the final `AnalysisReadyPhosphoDataset` boundary.

### Problem

The legacy preprocessing lane treated site-matrix construction as part of the path from prepared phospho/protein tables to an analysis-ready dataset. The rewrite currently constructs datasets directly from already-shaped phospho + site metadata tables and does not expose an equivalent supported site-matrix construction stage.

This leaves a scientific and structural hole: part of the legacy preparation logic still exists only in the archive, while the rewrite builder jumps straight to the final dataset contract.

### Why this matters

If site-matrix construction remains part of the intended ported preprocessing story, it should be owned by the builder preprocessing subsystem, not by workflows and not by ad hoc helpers.

### Goals

- restore site-matrix construction in the rewrite preprocessing path;
- make the policy explicit and testable;
- keep the final dataset boundary clean while still supporting the needed upstream construction logic.

### In scope

- site-matrix stage design and implementation
- public config surface for site-matrix policy
- integration into preprocessing pipeline
- dataset-boundary documentation updates

### Out of scope

- moving site-matrix logic into `AnalysisReadyPhosphoDataset` itself
- placing site-matrix science in kinase or signalome workflow code

### Required changes

- Introduce a dedicated site-matrix preprocessing stage.
- Decide what the supported public policy is for building or requiring site-matrix-ready inputs.
- Wire the stage into preprocessing ordering appropriately relative to missing-data handling and correction.
- Document how the stage relates to the final dataset boundary and which intermediate concerns remain private.

### Acceptance criteria

- the rewrite builder has an explicit supported path for site-matrix construction or enforcement;
- site-matrix logic lives in the preprocessing subsystem, not in workflows;
- tests cover both the supported success path and the failure path for unsupported/incomplete inputs;
- docs explain the supported site-matrix policy clearly.

### Test expectations

Add tests covering:

- stage execution on supported input shapes;
- ordering relative to other preprocessing stages;
- unsupported configuration or missing-input errors;
- donor-derived regression checks where practical.

---

## Ticket 6 — Restore legacy comparison-building in the supported dataset-builder preprocessing lane

**Priority:** P1  
**Type:** Scientific Correctness / Preprocessing / Public Contract  
**Area:**

- new preprocessing stage/module for comparison-building
- `src/phospy/api/configs.py`
- `src/phospy/validation/datasets/preprocessing.py`
- `src/phospy/datasets/builders/executor.py`
- `src/phospy/datasets/models.py`
- docs/tests around sample metadata and comparisons
- legacy donor references in `legacy_archive/phospy_legacy/preprocessing/steps.py`, `services.py`, and `dataset.py`

### Summary

Port the legacy comparison-building behavior into the supported rewrite builder lane so the package can construct dataset-level comparisons during preprocessing when that policy is requested and the required metadata is available.

### Problem

Legacy preprocessing could add pairwise comparisons as part of the dataset-conditioning path. The rewrite builder currently accepts sample metadata but does not expose a supported comparison-building stage that turns that metadata into dataset-level comparisons in the preprocessing lane.

That leaves another break in the raw/semi-processed-to-analysis-ready story.

### Why this matters

Comparison-building is part of the scientific conditioning story, not a cosmetic convenience. If the rewrite supports it, it needs a truthful builder-boundary contract. If it does not, the contract needs to say so clearly. This ticket assumes support is being restored.

### Goals

- restore comparison-building below the dataset boundary;
- make policy and prerequisites explicit;
- keep comparison logic out of workflow interpreters.

### In scope

- public config surface for comparison-building
- preprocessing stage implementation
- validation of required metadata/prerequisites
- docs and regression tests

### Out of scope

- moving comparison inference into workflow configuration
- hidden automatic comparison creation without policy control

### Required changes

- Introduce a dedicated comparison-building preprocessing stage.
- Define what sample metadata shape or explicit inputs are required.
- Decide whether comparisons are inferred, passed through, or both, and make that policy explicit in builder config.
- Wire the stage into the preprocessing pipeline at the correct point.
- Document the supported comparison-building story and failure modes.

### Acceptance criteria

- builder preprocessing can build supported comparisons from supported inputs;
- failure messages are clear when required sample metadata or grouping inputs are absent;
- workflow requests remain free of comparison-building policy;
- tests lock the restored behavior with donor-derived expectations where practical.

### Test expectations

Add tests covering:

- successful comparison construction from supported sample metadata;
- pass-through versus inferred behavior if both are supported;
- invalid metadata / missing-group errors;
- regression coverage against legacy donor behavior.

---

## Recommended implementation order

1. **Ticket 1** — fix the governance story first so planning is truthful.
2. **Ticket 2** — land the internal preprocessing subsystem.
3. **Ticket 3** — land the stable public config surface.
4. **Ticket 4** — total/protein correction.
5. **Ticket 5** — site-matrix construction.
6. **Ticket 6** — comparison-building.

This order keeps structure ahead of science and avoids pushing three new scientific features into the current one-file preprocessor.

## Suggested commit sequence

- `docs: reopen preprocessing science gaps in audit and donor inventory`
- `refactor: introduce staged dataset preprocessing subsystem`
- `feat: expand builder preprocessing config surface`
- `feat: restore total protein correction in dataset preprocessing`
- `feat: restore site matrix construction in dataset preprocessing`
- `feat: restore comparison building in dataset preprocessing`
