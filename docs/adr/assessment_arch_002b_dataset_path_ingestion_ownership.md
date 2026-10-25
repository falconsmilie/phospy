# ARCH-002B: Dataset Path Ingestion Ownership Assessment

## Record

- **Kind:** bounded architecture assessment; not an ADR
- **Assessment date:** 2026-10-25
- **Source snapshot:** `fcdb7a46692e17c1d7add296730aa05211c23223`
- **Recommendation:** retain the current narrowly classified Datasets-to-IO
  dependency
- **Authority effect:** none; this note neither accepts a new policy nor
  authorises a reader migration

## Current status and governing decision

The dependency remains present. `phospy.datasets.builder` imports and constructs
`phospy.io.readers.dataset_inputs.DatasetPathTableReader`. Static ownership
reporting therefore includes `datasets -> io`; IO also imports Dataset contracts
and models for importer handoff, bundle reconstruction, and publishing. The raw
owner graph has a classified cycle, while the policy-filtered graph used by the
DAG gate removes the exact classified import.

The 6 October 2026 Wave 5B authorisation is historical rationale. Its current
policy force comes from accepted and amended ADR-0053, the
`dataset-builder-path-reader` record in
`architecture/permanent-condensation-classifications.json`, and the executable
architecture guards. They permit exactly:

```text
phospy.datasets.builder
  -> phospy.io.readers.dataset_inputs.DatasetPathTableReader
```

ADR-0053 also requires the eager import closure to remain acyclic and states
that IO retains filesystem ownership. ADR-0054 separately keeps
`AnalysisReadyDatasetBuilder` as a Datasets-owned public object and makes direct
implementation-package imports unsupported external routes. The current
architecture documentation reports the raw and policy-filtered graphs
separately. A classified owner cycle is therefore an accepted tradeoff, not
evidence of a Python import failure.

## Responsibility and consumer findings

The concrete reader is a small boundary adapter with mixed context but a clear
primary responsibility:

- IO owns file-format detection, filesystem access, CSV/TSV/Parquet parsing,
  numeric and metadata parsing policy, and low-level read errors in
  `phospy.io.readers.tables`.
- The adapter selects those IO readers using Dataset field names (`phospho`,
  `site_metadata`, `sample_metadata`, and `total`), falls back to the generic
  table reader for other Dataset inputs, and adds Dataset-build error context.
- Datasets owns the `DatasetPathTableReader` protocol, DataFrame-or-path source
  validation, path normalization, optional-input handling, input-resolution
  order, and the rest of Dataset construction.

The builder is the concrete reader's only production consumer. The
`phospy.io.readers` initializer re-exports the class, but that implementation
route is not in the governed public API. Dataset importer implementations do
not consume this reader: importers own their source ingestion and hand the
builder DataFrames through `DatasetBuildRequest`. IO's reverse dependencies on
Datasets are broader and independent of this adapter, so changing this one
class would not make IO independent of Dataset contracts or models.

## Options

| Concern | Retain the classified IO adapter | Move the reader into Datasets |
| --- | --- | --- |
| Responsibility cohesion | Keeps filesystem and format behavior with IO. The small field-dispatch and error-context layer is appropriate boundary knowledge for a Dataset-specific IO adapter. | Improves proximity to Dataset field semantics, but places a concrete filesystem adapter under Datasets against current policy. |
| Change coupling | Dataset source semantics can change behind its protocol; IO table parsing remains reusable by other consumers. The exact classification exposes the coupling. | Moving only the class still imports `phospy.io.readers.tables`, so the owner edge and cycle remain under a different source module. Removing the edge requires moving or duplicating the underlying table readers and supported-format policy as well. |
| Testing | Existing source-boundary tests cover Dataset validation and injection; IO reader tests cover parsing plus end-to-end builder path reads; architecture tests cover the exact edge and import orders. | Tests would need ownership relocation and new protection against divergence from the general table readers. A full move expands the behavior under migration well beyond this adapter. |
| Import cost | Importing the public builder loads the narrow IO reader/table path. The eager closure has no path back to Datasets, and both import orders pass. No import-time regression is currently demonstrated. | A class-only move retains the IO table imports and gives no material import reduction. Moving all parsing could avoid loading the IO namespace but would load equivalent pandas/file-format code under Datasets; no measured benefit supports that cost. |
| Public object identity | No change. The public builder continues to resolve to `phospy.datasets.builder:AnalysisReadyDatasetBuilder`. The concrete reader is absent from the public API baseline. | Builder identity can be preserved, and reader identity has no supported public guarantee. Retaining an old reader alias would add maintenance without a governed compatibility need. |
| Maintenance cost | Retains one exact classification and its narrow architecture/import-order guards. | Requires source, owner inventory, classification, architecture tests, reader tests, and documentation changes. A full edge removal also requires an explicit policy change for filesystem ownership. |

## Behaviour that ownership changes must preserve

The current contract is distributed across request validation, input
resolution, and the concrete adapter:

- `DatasetBuildRequest` accepts DataFrames, strings, `pathlib.Path`, and generic
  `os.PathLike[str]` values. A generic PathLike becomes a `Path`. String paths
  retain identity during validation and are stripped only when converted to a
  `Path` for reading; blank strings fail before any reader call with
  `dataset build request <field> path cannot be empty`.
- Required `None` values fail with
  `dataset build request <field> is required`. Optional values use
  `allow_none=True` during request validation, and the source resolver returns
  `None` without calling the reader.
- Request validation precedes interpretation and reading. In the site-level
  lane it validates `phospho`, then `site_metadata`; in the peptide-evidence
  lane it validates `peptide_evidence`, then optional `peptide_site_mapping`.
  Both lanes then validate optional `sample_metadata` and `total` before later
  request policy checks. That ordering determines which error is observed.
- The concrete adapter uses schema-aware readers for `phospho`,
  `site_metadata`, `sample_metadata`, and `total`; other field names use
  `read_table`. `UnsupportedInputFormatError` remains that exact exception type
  and gains the Dataset field, path, and supported formats. `PhosPyInputError`
  remains that exact type and gains Dataset field/path context.
- `DatasetInputReader` propagates exceptions from injected source validators and
  path readers. It requires an injected path reader for path inputs but returns
  an injected DataFrame directly.
- Current fallback semantics use truthiness in selected constructors. A
  false-valued source validator passed to `DatasetInputReader` is replaced by
  its default, while its `path_reader` is stored directly. A false-valued
  `path_reader` passed through
  `AnalysisReadyDatasetBuilder._with_components(...)` is replaced by the
  concrete default because builder composition uses `path_reader or
  DatasetPathTableReader()`. Any future refactor must preserve these observable
  choices unless it deliberately changes and tests the injection contract.

## Recommendation and consequences

Retain the current ownership and classification. It preserves the useful
separation between Dataset-owned input semantics and IO-owned filesystem/table
parsing, and it avoids a broad migration that would be necessary to remove the
owner edge in fact. Current tests demonstrate an acyclic eager module graph and
stable opposing import orders, so there is no runtime-cycle defect to repair.

The consequence is that the complete static owner graph continues to report the
cycle and the policy-filtered DAG continues to exclude the one exact D05B-1
edge. That cost is explicit, narrow, and enforced. Reconsider ownership only if
there is a measured import/runtime problem, the adapter gains substantial
Dataset-only behavior, or policy deliberately moves schema-aware filesystem
parsing out of IO. No reader-migration follow-up is recommended from current
evidence.

## Evidence and uncertainty

Current source/import searches found no additional production consumer of the
concrete reader. Focused source-boundary, table-reader, exact-edge,
classification, eager-graph, and opposing-import-order checks passed (51 tests,
Python 3.11). The assessment did not benchmark import time, and direct external
use of the unsupported `phospy.io.readers` route cannot be ruled out; neither is
a current public-contract constraint.
