# Parity Model

This document explains what parity means in PhosPy, where parity evidence currently exists, and where it does not.

PhosPy is an unofficial Python port of selected PhosR workflow components. Some parts of the project are tested directly
against R-generated fixtures, while newer native workflow pieces are better described as PhosR-style implementations
with targeted fixture-backed seams rather than blanket claims of numerical equivalence.

## What Parity Means Here

In this repository, parity means that a specific Python path has been compared against outputs generated from the R
package and shown to agree within the limits of the fixture-backed tests for that path.

Parity in this project is therefore:

- explicit
- test-backed
- limited to named workflow seams
- narrower than full package equivalence

It does **not** mean that the repository as a whole is a complete behavioural, numerical, or feature-level replacement
for PhosR.

## What Parity Does Not Mean

Parity should not be read as a claim of:

- full feature parity with the original R package
- full behavioural parity across all PhosR methods
- full output parity for workflows that have not been fixture-backed
- broad numerical equivalence for the newer native kinase workflow outside the seams documented here
- a native Python replacement for `Signalomes()`

Where parity has not been established, PhosPy should be described as an evolving Python port inspired by and translating
selected parts of the PhosR workflow.

## Current Parity-Backed Areas

The strongest current parity evidence in this repository is for:

- deterministic preprocessing and matrix-building seams backed by small synthetic fixtures
- downstream kinase-analysis summaries backed by R-generated fixtures
- selected native kinase workflow seams where Python outputs are checked against reference expectations captured in
  fixtures and tests

These checks provide confidence in the implemented paths, but they are still narrower than a claim of package-wide
equivalence.

## Fixture Paths

This repository currently uses two main R fixture paths.

### Small Synthetic Fixture Path

Use this path for deterministic preprocessing and core matrix-building parity.

Generate fixtures with:

```bash
Rscript scripts/generate_r_fixtures.R
```

This writes CSV fixtures into `tests/fixtures/r_reference/` for:

- corrected phosphosite values
- PhosR input rows and site matrix
- `predMat`
- weighted kinase activity
- KSEA scores
- substrate counts
- `sessionInfo()` for provenance

This path is useful for logic-level parity and regression protection in the core preprocessing and downstream summary
flow. It should not be treated as strong evidence for broader downstream equivalence beyond the implemented and tested
wrapper path.

### Bundled PhosR L6 Fixture Path

Use this path for a more realistic downstream kinase-analysis parity check based on PhosR’s bundled rat L6 myotube
example dataset, which is used throughout the original package examples and vignette.

Generate fixtures with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

This writes CSV fixtures into `tests/fixtures/r_reference_l6/` for:

- the filtered standardised L6 phosphosite matrix used for kinase analysis
- `predMat`
- weighted kinase activity
- KSEA scores
- kinase target counts
- `sessionInfo()` for provenance

That directory also includes committed native-workflow reference tables used in parity tests for Python-native scoring
and prediction seams.

This is the better current evidence for parity of the implemented downstream kinase-analysis methods.

## Native Kinase Workflow

PhosPy now includes a native end-to-end kinase workflow covering:

- substrate-profile construction
- motif scoring
- profile scoring
- weighted score combination
- candidate-substrate selection
- adaptive SVM prediction

This workflow is implemented natively in Python and is intended to provide a coherent PhosR-style kinase scoring path.

At present, it should be described carefully:

- it is a live native workflow, not just a thin wrapper
- parts of it are tested and fixture-backed
- it is **not** yet a blanket claim of numerical equivalence to the R package
- parity claims for this path should stay limited to the specific seams that are explicitly covered by fixtures and
  tests

That distinction matters. The project can be both useful and scientifically careful at the same time.

## Parity-Sensitive Prediction Settings

The prediction layer also exposes `svm_mode="r_parity"`.

That setting is useful when you want learner behaviour that more closely tracks the R seam used by PhosR. It helps keep
parity-sensitive prediction settings explicit and testable rather than hidden behind defaults. It should still be
treated as one seam within a larger workflow, not as a standalone claim of full workflow equivalence.

## Prediction Trace Directories

The repository also contains prediction-trace exports for seam-level debugging.

### Committed R Trace Fixtures

The directory `tests/fixtures/r_reference_l6/prediction_trace/` contains committed R-side prediction traces for selected
kinases on the bundled L6 path.

These files are useful for inspecting candidate selection, negative sampling, iteration-level probabilities, and final
ensemble predictions at a more granular level than end-result tables alone.

### Committed Python Trace Fixtures

The directory `tests/fixtures/python_reference_l6/prediction_trace/` contains committed Python-side prediction traces
exported from the native predictor.

These files exist so that Python learner-stage behaviour can be compared directly against the committed R traces. They
are part of the repository’s seam-level debugging and comparison story, and they help make prediction-stage differences
easier to inspect and discuss.

### Temporary Trace Output

The directory `tmp_trace_out/` exists as a documented scratch area for ad hoc Python prediction-trace exports during
debugging.

Its purpose is different from the committed fixture paths:

- it is for short-lived or investigation-specific trace runs
- its contents may be regenerated or replaced freely
- it should not be treated as stable reference data
- files should only be considered part of the parity contract when they are intentionally promoted into a committed
  fixture path under `tests/fixtures/`

This distinction helps keep the repository honest about which artefacts are stable evidence and which are working-output
diagnostics.

## Running Parity Tests

If the fixtures are present, parity tests can be run with:

```bash
pytest -m parity
```

These tests should be treated as the executable definition of the repository’s current parity contract.

## Maintenance Rule

When a parity-backed workflow changes, at least one of the following should also change in the same line of work:

- the fixtures
- the tests
- the documented scope of the parity claim

Do not silently broaden parity claims in the README or other project documentation beyond the fixture-backed seams
described here without adding corresponding evidence.

When temporary debugging artefacts in `tmp_trace_out/` reveal an important new seam or a stable comparison case, promote
the relevant outputs into a committed fixture directory and update the associated tests or documentation accordingly.

## Provenance

Generated fixtures should keep provenance information, including `sessionInfo()` where available, so that reference
outputs can be tied back to the R environment used to create them.