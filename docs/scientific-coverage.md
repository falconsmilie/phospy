# Scientific Coverage

PhosPy implements selected PhosR-style phosphoproteomics workflows. It does not
claim full package equivalence with PhosR.

## Current Supported Lane

The supported public lane is:

1. build an `AnalysisReadyPhosphoDataset`
2. run kinase scoring and prediction
3. optionally run signalome analysis from the kinase result

Bundled runtime references in `1.5.0` are rat-only. Human and mouse analysis can
be run by passing an explicit `ReferenceBundle` in Python.

## Scientific Confidence Labels

Use these labels when discussing coverage:

| Label | Meaning |
| --- | --- |
| `PARITY_GATED_ACTIVE_SCIENCE` | protected by active fixture-backed parity tests |
| `PHOSPY_VALIDATED_SCIENCE` | validated by PhosPy contract, unit, and integration tests |
| `SUPPORTED_CONTRACT_CHANGED` | intentionally supported with a changed public contract |
| `OPEN_GAP` | not yet covered or not yet claimed |

## Active Coverage

| Area | Current status |
| --- | --- |
| Dataset boundary | strict PhosPy dataset contract |
| Kinase scoring and prediction | active fixture-backed parity and workflow tests |
| Activity output | supported as a thresholded substrate-mean and weighted activity output; not full KSEA enrichment |
| Signalome workflow | supported from kinase result with explicit `protein_id` |
| Output publishing | supported simple publishers and reloadable bundle services |
| Human/mouse bundled references | open gap for bundled runtime data in this release |

## Interpretation Limits

- Activity output named `thresholded_substrate_mean_activity` is a simple summary
  over predicted substrates above threshold.
- Rank-weighted fusion scores combine profile-correlation and motif-frequency
  evidence using rank-derived weights.
- Signalome module/network scores are derived summaries, not probabilities,
  calibrated confidence values, or causal proof.
- Missing kinase correlations stay missing. `0.0` means a finite near-zero
  correlation was estimated.

## Where Details Live

- [Parity](parity.md) tracks PhosR comparison evidence and fixture locations.
- [Performance Contracts](performance.md) covers scale limits.
- [ADR Index](adr/index.md) stores maintainer decision records.
