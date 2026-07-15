# ADR-0035: Provenance Immutability and Stable Serialization

## Status

- **ADR ID:** ADR-0035
- **Title:** Provenance Immutability and Stable Serialization
- **Status:** Accepted
- **Date:** 2026-07-15
- **Decision Type:** Result and Dataset Provenance Contract

## Context

Run, dataset, reference, preprocessing, and policy provenance are used for
reproducibility, diagnostics, saved-bundle round trips, and stable hashing.
Many provenance fields carry JSON-like dictionaries and arrays supplied by
callers or assembled by workflow components.

Frozen dataclasses and type annotations do not protect nested containers. A
caller can mutate a nested dictionary after construction, mutate a list stored
inside provenance, or mutate a serialized payload if serialization exposes an
internal alias. Any of those changes can invalidate reproducibility claims and
make stable hashes dependent on object ownership accidents.

## Decision

All JSON-like provenance values must be normalized through the
`phospy.provenance` recursive freezer at construction boundaries.

The freezer:

- defensively copies constructor input;
- stores dictionaries as recursively immutable mappings;
- stores arrays as immutable sequences;
- accepts only JSON-compatible scalar values, mappings, and non-string
  sequences;
- rejects unsupported objects explicitly; and
- rejects non-finite floats explicitly.

Serialization must thaw provenance into fresh JSON payloads every time. Public
payloads remain ordinary `dict` and `list` containers, but those containers are
copies and are not aliases of internal provenance state.

Stable hashing and round-trip serialization must operate on the same validated
JSON-compatible value space. Hashing therefore validates and serializes fresh
payloads with sorted keys and no non-finite number support.

## Consequences

Code that reads stored provenance must treat arrays as immutable sequences, not
as concrete `list` instances. Code that needs a mutable JSON payload must call a
serialization function and mutate that fresh payload, not the provenance object.

Provenance constructors are now the enforcement boundary for JSON-like
immutability. Workflow assemblers may still normalize domain objects into JSON
facts before constructing provenance, but unsupported values that reach a
provenance dataclass fail at construction rather than being stringified.

This decision complements ADR-0016's dataframe ownership boundary: dataframe
inputs are copied at public boundaries, while JSON-like provenance inputs are
recursively frozen and serialized from fresh payloads.

## Implementation Notes

The recursive freezer lives in `src/phospy/provenance/immutability.py`.
Provenance dataclasses with JSON-like mapping fields route those fields through
the freezer in `__post_init__`. Serialization helpers route frozen provenance
through the thawing path before returning payloads.

The current audited fields include:

- `ReproducibilityCaveat.details`;
- `TableFingerprint.index_structure` and `column_index_structure`;
- environment dependency, platform, BLAS/LAPACK, thread, locale, and
  constraints mappings;
- preprocessing stage `parameters` and `diagnostics`;
- batch-correction provenance mappings and rejected-entity `details`;
- reference and Kinase Library resource `sequence_window`, `source_files`, and
  `manifest` mappings;
- run-level `workflow_parameters`;
- scientific-policy `parameters`; and
- derived-quantitative missingness, matrix-transformation, and parameter
  mappings.
