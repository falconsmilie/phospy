# Validation and PhosR Parity

This is the short guide to:

- the validation rules you are most likely to hit first
- what PhosPy means by parity to the R `PhosR` package

For the fuller parity guide, see [`docs/parity.md`](parity.md).

## Quick Checks

From the repository root:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

Matching Make targets:

```bash
make pre-commit
make test-unit
make test-parity
```

## Validation Rules You Will Notice First

### Input Tables

- total input must include `genes` plus `group1` to `group6`
- phospho input must include `uid`, `gene_names`, `gene_p_site`, `localization_prob`, `centralized_sequence`, and `p_group1` to `p_group6`
- total and phospho files are read as TSV
- `predMat` is read as CSV with the first column used as the phosphosite index
- file-loaded total and phospho headers are cleaned to lowercase snake case before validation
- duplicate cleaned headers are rejected
- `gene_p_site` must split cleanly into gene and site parts, such as `BTK_Y551`
- `localization_prob` and `predMat` scores must stay in `[0, 1]`

### Compatibility Checks

- total, phospho, and corrected value columns must align in count
- by default, protein correction allows no silent phosphosite row loss
- downstream kinase analysis requires overlap between `predMat` and the phosphosite matrix
- native workflow inputs must share phosphosite IDs across the matrix, substrate map, and sequence inputs
- motif-aware native workflow runs require both `motif_sequences` and matching `site_sequences`

### Request Checks

- file paths must exist and point to files
- comparison pairs must use known schema groups and must not be duplicated
- `substrate_map` must not be empty
- `site_sequences` must be a phosphosite-keyed mapping or a pandas `Series`
- profile-only native prediction must opt in with `allow_profile_only_fallback=True`

## What Parity Means Here

In this repository, parity means parity to the R `PhosR` package for a specific tested seam.

That means all three of these are true:

- a committed fixture exists for the seam being discussed
- an automated test compares the Python output against that fixture
- the claim stays limited to that seam

It does **not** mean:

- the whole package is numerically identical to `PhosR`
- every `PhosR` workflow branch or option is implemented
- every native Python path is expected to match `PhosR`

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is part of the supported public API, but it remains a native Python workflow.

- `svm_mode="default"` is the normal native mode
- `svm_mode="r_parity"` is a narrower learner-seam mode used for selected parity checks against `PhosR`

Example:

```python
from phospy import KinaseWorkflow

native = KinaseWorkflow(svm_mode="default")
comparison = KinaseWorkflow(svm_mode="r_parity")
```

## Read Next

- [`docs/parity.md`](parity.md) for the detailed parity guide
- [`docs/api.md`](api.md) for method signatures and validation by API entry point
- [`docs/fixtures.md`](fixtures.md) for the fixture and trace layout
