# Package Layout

PhosPy is organised by domain capability first.

Use this as the contributor-facing guide for where new code should live. The migration map in this directory is historical context, not the main guide.

## Domain Packages

| Package | Owns |
| --- | --- |
| `phospy.api` | Supported public workflows and thin public orchestration |
| `phospy.datasets` | Dataset models, schema objects, loaders, and dataset-bound builders |
| `phospy.preprocessing` | Raw-to-analysis-ready preprocessing, including filtering, correction, site-matrix preparation, and preprocessing result models |
| `phospy.prediction` | Prediction engines, scoring, `predMat` execution, and prediction result models |
| `phospy.activities` | Kinase activity analysis and activity result models |
| `phospy.signalomes` | Signalome construction plus map and network outputs |
| `phospy.references` | Bundled biological reference assets, species and reference resolution, and reference-facing models |
| `phospy.io` | Shared table reading, mapping-file loading, publishing, and writer helpers |
| `phospy.validation` | Validation grouped by validation type |
| `phospy.errors` | Shared application error classes |
| `phospy.internal` | Narrow internal-only constants and type aliases |

Preprocessing is intentionally split across a small number of seams:

1. `DatasetPreprocessing` binds one dataset workspace to the preprocessing path.
2. `CoreProcessor` orchestrates full or phospho-only preprocessing.
3. Step services and `SiteMatrixBuilder` perform the concrete transforms.
4. `AnalysisReadyDatasetBuilder` adapts user-shaped inputs into the analysis-ready boundary.

## Root Package Policy

`phospy.__init__` is intentionally small.

The root package keeps a supported convenience surface for a limited set of high-level workflow and result types. New code should prefer the owning domain package directly.

Use root imports only when a short user-facing entry point is genuinely helpful, such as examples or interactive sessions.

## Placement Rules

- Put implementation code in the domain package that owns the behaviour.
- Keep `phospy.api` thin. It coordinates work; it does not own preprocessing, prediction, references, activity analysis, or signalome logic.
- Keep package `__init__` surfaces narrow.
- Keep `phospy.internal` narrow. Do not turn it into a helper bucket.
- Put shared file and table access in `phospy.io` only when it is genuinely cross-domain and not scientific logic.
- Prefer adding a small focused module to an existing domain package over creating a new catch-all bucket.

## Scientific Policy Placement

Keep explicit scientific policies with the domain code that owns the behaviour.

- duplicate-site policy belongs with preprocessing and site-matrix construction
- profile missing-value policy belongs with kinase profile aggregation
- signalome module-selection policy belongs with signalome clustering

These policies should stay visible in the public workflow or builder boundary that owns them so PhosR parity and intentional divergence remain reviewable.

## DataFrame Ownership Rule

Preprocessing and validation follow one explicit ownership rule:

1. External DataFrame boundaries may copy once when taking ownership.
2. Internal preprocessing services then work on owned tables without another whole-frame defensive copy.
3. Snapshot-style outputs must say clearly whether they reuse owned internal tables or detach copies.

In practice this means:

- schema validators keep their default defensive copy at public input boundaries
- internal fast paths use explicit owned constructors such as `prepare_owned()`, `correct_owned()`, `build_owned()`, and `AnalysisReadyPhosphoDataset.from_owned()`
- `CoreProcessor.process()` and `process_phospho_only()` copy caller-managed input tables once at entry, then stay on the owned path

This keeps mutation boundaries reviewable and reduces repeated full-frame copy churn through preprocessing.

## Related Documents

- [`docs/api.md`](../api.md) for the supported Python surface
- [`docs/adr/0004-reorganise-by-domain.md`](../adr/0004-reorganise-by-domain.md) for the architecture decision
- [`docs/architecture/root-package-migration-map.md`](./root-package-migration-map.md) for the historical migration record
