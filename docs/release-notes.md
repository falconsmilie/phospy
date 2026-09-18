# PhosPy Release Notes

## Version 1.7.5

Release date: 2026-09-18.

These notes describe only the changes since Version 1.7.4.

## Release Overview

PhosPy 1.7.5 completes the package's SPS discovery and RUV-III path:

```text
reference phosphoproteomics datasets
    -> SPS discovery
    -> ControlSiteSet
    -> ruv_iii_style correction
    -> corrected dataset
```

SPS discovery is a separate, explicit operation over suitable independent
reference phosphoproteomics evidence. It is not mandatory for correction:
existing governed caller-supplied controls remain valid. The new
`ruv_iii_style` method is also explicitly opt-in, so existing analyses are not
automatically migrated.

The established `method="sps_ruv_style"` method predates this release. It
remains the PhosPy-native SPS/control-based unwanted-variation estimator and is
not retrospectively described as RUV-III. Version 1.7.5 extends the existing
infrastructure with discovery and a separate replicate-aware RUV-III route; it
does not replace `sps_ruv_style`.

## SPS Discovery

`SpsDiscoveryWorkflow` discovers stable phosphosite candidates across multiple
prepared reference datasets. References share governed, protein-scoped
`site_key` identity and must declare coherent organism, biological baseline,
reference context, source identity, and established condition-relative log2
semantics. The caller remains responsible for choosing scientifically
appropriate reference evidence and its biological control baseline.

Discovery provides explicit configuration, per-reference contribution and
attrition records, deterministic consensus ranking, serializable provenance,
and a `ControlSiteSet` suitable for correction. Partial-reference candidates
are ordered conservatively by descending contributing-reference count, then by
the Fisher-style consensus score and governed `site_key`.

The discovery identity is immutable and evidence-sensitive. It passes through
the generated `ControlSiteSet` into correction provenance, so identical selected
keys derived from different reference evidence retain distinct lineage.

## Native SPS/RUV-Style Replicate-Aware RUV-III Compatibility

`method="ruv_iii_style"` requires replicate-set metadata and uses replicate
membership directly in estimation. Protected biological-condition terms remain
explicit. Singleton replicate sets are rejected, as is a requested `k` that is
not estimable from the validated controls and replicate structure; PhosPy does
not silently cap it. These are deliberate, stricter input and error contracts,
not failures of the RUV-III mathematics.

Governed row-median completion can be used internally by `ruv_iii_style` for
estimation, after which original missing positions are restored. Normal
analysis-ready construction still requires complete numeric output, so upstream
missing-data processing remains necessary. This behavior is a PhosPy contract
and is not an external PhosR missing-data parity claim.

The pinned external fixtures establish only bounded evidence:

- SPS ranking and selection agree with PhosR `getSPS` in the tested
  complete-reference supported domain.
- RUV-III corrected matrices agree with pinned `ruv::RUVIII` for the tested
  complete-data finite-`k` cases.
- Partial-reference SPS ordering, missing-data execution, `sps_ruv_style`, and
  the complete PhosR `RUVphospho` workflow are outside those parity claims.

Scientific validation covers planted nuisance-factor removal, protected-signal
retention, downstream differential-analysis compatibility, lineage and
serialization, and SPS-specific performance contracts. Normal execution and
tests consume the pinned fixture outputs without requiring R or network access.

See [SPS Discovery and Native SPS/RUV-Style Correction](api/sps-ruv.md),
[Parity](parity.md), and [Scientific Coverage](scientific-coverage.md) for the
complete contracts and evidence boundaries.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
