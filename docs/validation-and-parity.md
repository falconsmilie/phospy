# Validation and Parity

This is the short guide to two things:

- what PhosPy validates
- what its current parity claim does and does not cover

If you want the detailed parity guide, including metrics flags, cross-platform commands, and sample output, see
[`docs/parity.md`](parity.md).

## At a Glance

PhosPy uses three practical checks:

1. **Non-parity tests** check package behaviour without relying on R fixtures.
2. **Parity tests** compare selected Python seams against committed R/PhosR-backed reference outputs.
3. **Pre-commit checks** keep the repository consistent.

Run them from the repository root:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

If you prefer Make targets, the matching shortcuts are:

```bash
make pre-commit
make test-unit
make test-parity
```

## Validation Rules You Will Notice First

### Table validation

- total input requires `genes` plus `group1` to `group6`
- phospho input requires `uid`, `gene_names`, `gene_p_site`, `localization_prob`, `centralized_sequence`, and
  `p_group1` to `p_group6`
- file-loaded total and phospho headers are normalised to lowercase snake case before validation
- duplicate raw headers that collapse to the same cleaned name are rejected
- `gene_p_site` must split cleanly into gene and site parts, such as `BTK_Y551`
- `localization_prob` and `predMat` scores must stay in `[0, 1]`
- `predMat` must use a unique, non-null phosphosite index
- site matrices must use a unique, non-null phosphosite index

### Compatibility validation

- total, phospho, and corrected value columns must align in count
- protein correction refuses mismatched joins by default through `max_unmatched_fraction=0.0`
- protein correction also refuses duplicate protein identifiers in the filtered total table
- downstream kinase analysis requires `predMat` and the phosphosite matrix to overlap by at least one row and by at
  least 10% of the phosphosite matrix
- native workflow inputs must share phosphosite IDs across the matrix, substrate map, and sequence inputs
- motif-aware native workflow runs require both `motif_sequences` and matching `site_sequences`

### Request validation

Request models also check practical edges such as:

- input paths must exist and be files
- comparison pairs must use known groups and must not be duplicated
- native workflow requests must not use an empty `substrate_map`
- `site_sequences` must be a phosphosite-keyed mapping or a pandas Series with a phosphosite index
- profile-only native prediction must opt in with `allow_profile_only_fallback=True`

## What Parity Means in PhosPy

In this repository, parity means all three of these are true:

- a committed fixture exists for the seam being discussed
- an automated parity test checks Python output against that fixture
- the claim stays limited to that seam

Parity does **not** mean:

- the whole package is numerically identical to PhosR
- every PhosR option or workflow branch is implemented
- every Python-native path should match the R implementation exactly

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is part of the supported public API, but it is still a native Python workflow.

- `svm_mode="default"` is the preferred Python-native mode
- `svm_mode="r_parity"` is a narrower learner-seam comparison mode backed by committed parity fixtures

That option narrows one comparison seam. It does not turn the whole workflow into a package-wide PhosR-equivalence
claim.

Configuration example:

```python
from phospy import KinaseWorkflow

native = KinaseWorkflow(svm_mode="default")
comparison = KinaseWorkflow(svm_mode="r_parity")
```

## Where to Go Next

- Read [`docs/parity.md`](parity.md) for the full parity guide.
- Read [`docs/fixtures.md`](fixtures.md) for the fixture and trace directory layout.
- Read [`CONTRIBUTING.md`](../CONTRIBUTING.md) for the local development and release workflow.
