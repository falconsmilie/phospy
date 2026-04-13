# Package layout

PhosPy is organised by domain capability first.

This is the contributor-facing package layout to use when adding or moving code.
The older migration map remains in this directory as refactor history, not as the
main guide for where new code should live.

## Domain packages

| Package | Owns |
| --- | --- |
| `phospy.api` | Supported public workflows and other thin public orchestration |
| `phospy.datasets` | Dataset models, schema objects, loaders, and dataset-bound builders |
| `phospy.preprocessing` | Raw-to-analysis-ready preprocessing, including filtering, correction, site-matrix preparation, and preprocessing result models |

Preprocessing is intentionally organised around a small number of seams:

1. `DatasetPreprocessing` binds one dataset workspace to the preprocessing path
2. `CoreProcessor` orchestrates full or phospho-only preprocessing
3. step services and `SiteMatrixBuilder` perform the concrete transforms
4. `AnalysisReadyDatasetBuilder` adapts user-shaped inputs into the analysis-ready boundary
| `phospy.prediction` | Prediction engines, scoring, `predMat` execution, and prediction result models |
| `phospy.activities` | Kinase activity analysis and activity result models |
| `phospy.signalomes` | Signalome construction plus map and network outputs |
| `phospy.references` | Bundled biological reference assets, species and reference resolution, and reference-facing models |
| `phospy.io` | Shared table reading, mapping-file loading, publishing, and writer helpers |
| `phospy.validation` | Validation grouped by validation type |
| `phospy.errors` | Shared application error classes |
| `phospy.internal` | Narrow internal-only constants and type aliases |

## Root package policy

`phospy.__init__` is intentionally small.

The root package keeps a supported convenience surface for a limited set of
high-level workflow and result types. New code should prefer importing from the
owning domain package directly.

Use root imports only when a short, user-facing entry point is genuinely helpful,
for example in simple examples or interactive sessions. Prefer domain-package
imports for implementation code, tests that exercise internals, and new
contributor-facing documentation.

## Placement rules

- Put implementation code in the domain package that owns the behaviour.
- Keep `phospy.api` thin. It coordinates work; it does not own preprocessing,
  prediction, reference resolution, activity analysis, or signalome logic.
- Keep package `__init__` surfaces narrow. Do not re-export trusted input bundles,
  low-level validators, or orchestration-only helpers from convenience barrels.
- Keep `phospy.internal` narrow. Do not turn it into a general-purpose helper
  bucket.
- Put shared file and table access in `phospy.io` only when it is genuinely
  cross-domain and not part of a scientific domain.
- Prefer adding a small focused module to an existing domain package over
  creating a new cross-cutting bucket.

## Related documents

- [`docs/api.md`](../api.md) for the supported Python surface
- [`docs/adr/0004-reorganise-by-domain.md`](../adr/0004-reorganise-by-domain.md)
  for the architecture decision
- [`docs/architecture/root-package-migration-map.md`](./root-package-migration-map.md)
  for the historical migration record


## Scientific policy placement

Keep explicit scientific policies with the domain code that owns the behaviour.
Do not hide them in generic utility modules.

- duplicate-site policy belongs with preprocessing and site-matrix construction
- profile missing-value policy belongs with kinase profile aggregation
- signalome module-selection policy belongs with signalome clustering

These policies should stay visible in the public workflow or builder boundary that
owns them so PhosR parity and intentional divergence remain reviewable.
