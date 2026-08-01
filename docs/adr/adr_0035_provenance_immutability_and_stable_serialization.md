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

Update note (2026-07-17, recursive immutability hardening): provenance
immutability cannot be implemented by subclassing `dict` or `list`, wrapping
caller-owned storage with `MappingProxyType`, relying on frozen dataclasses, or
stringifying keys during construction. Those approaches either leave
base-class mutation paths open, retain aliases to mutable caller state, or
collapse distinct key identities before stable serialization.

Update note (2026-07-17, trusted construction assertion linkage): Trusted
dataset construction assertions are immutable provenance values. Their nested
evidence details are recursively frozen, serialized from fresh payloads, and
linked into trusted-table reconstruction provenance by a stable assertion fingerprint.
`from_trusted_tables(...)` also verifies supplied table fingerprints against
the actual constructed tables; table hashes, shape, axis alignment, organism
coherence, and `site_sequence` presence are mechanically checked facts rather
than waivable assertion text.

Update note (2026-07-20, scale and kinase attrition evidence closure):
`IntensityScaleEstablishmentProvenance.parameters` and
`KinaseWorkflowAttritionProvenance` JSON evidence (`metrics`, `policy`, and each
`policy_violations` entry) now use the shared recursive freezer directly at
their model boundaries. The transformation and result contracts retain typed
top-level fields while rejecting invalid JSON keys, unsupported nested objects,
and non-finite floats in evidence mappings.

Update note (2026-07-21, quantitative-meaning transition provenance):
`QuantitativeMeaningTransitionProvenance` is a schema-versioned immutable domain
object. Its `parameters`, `input_table_fingerprints`, and diagnostic caveat
codes are frozen at construction and serialized into fresh stable payloads.
New `IntensityScaleState` payloads include `quantitative_meaning_provenance`.
Historical payloads that lack the field may be loaded only through the explicit
legacy migration path, which labels the evidence as `legacy_unverified` and must
not claim that a derived operation occurred.

Update note (2026-08-01, trusted construction assertion schema v4):
`TrustedDatasetConstructionAssertions` schema version 4 adds the optional
`numeric_semantic_domain` evidence/waiver dimension. The seven existing
trusted-construction dimensions remain required for complete trusted table
construction. A numeric-semantic-domain waiver is serialized, included in
`waived_assertions`, and contributes to the assertion fingerprint so trusted
construction cannot hide a scale/meaning/value-domain conflict.

Update note (2026-07-30, normalized fingerprint axis semantics): provenance
fingerprints that intentionally normalize row and column order use one central
axis-label policy in `phospy.provenance.hashing`. Labels are sorted by the
collision-safe typed sort key `typed-axis-label-sort-v1`. Supported labels are
non-missing strings, integers, and tuple/MultiIndex labels composed only of
supported labels. Integer labels are ordered numerically and remain distinct from
string labels with the same display form, such as `1` and `"1"`. Duplicate typed
axis labels and unsupported labels fail with `ProvenanceFingerprintError` before
hashing; normalized provenance never falls back to caller-supplied order.

## Decision

All JSON-like provenance values must be normalized through the
`phospy.provenance` recursive freezer at construction boundaries.

The freezer:

- defensively copies constructor input;
- stores JSON objects as a private-storage immutable mapping that does not
  inherit from `dict`;
- stores arrays as tuples;
- accepts only JSON-compatible scalar values, mappings, and non-string
  sequences;
- accepts only `str` JSON object keys;
- rejects duplicate keys yielded by a mapping before construction;
- rejects unsupported objects explicitly; and
- rejects non-finite floats explicitly.

JSON object keys are never coerced with `str()`. A non-string key is invalid,
including keys that would otherwise collide after stringification, such as
`1` and `"1"`. A custom mapping that yields the same string key more than once
is also invalid. This is the only supported key policy for constructor input,
deserialized payloads, and stable JSON hashing.

Serialization must thaw provenance into fresh JSON payloads every time. Public
payloads remain ordinary `dict` and `list` containers, but those containers are
copies and are not aliases of internal provenance state. Documented public
provenance JSON accessors such as run `workflow_parameters`, preprocessing
stage `parameters` and `diagnostics`, and batch-correction provenance mapping
fields follow the same rule: each read returns a fresh JSON-shaped payload while
the dataclass storage remains frozen. Serialized JSON objects always have
string keys because provenance objects cannot contain any other key type.

Stable hashing and round-trip serialization must operate on the same validated
JSON-compatible value space. Hashing therefore validates and serializes fresh
payloads with sorted string keys and no non-finite number support. It rejects
unsupported objects and invalid JSON object keys instead of normalizing them by
representation or stringification.

Order-normalized table fingerprints must normalize both axes before computing
the exact and tolerance hashes, and both hashes must use that same normalized
view. Axis normalization belongs to `phospy.provenance.hashing`; workflows such
as kinase and signalome provenance assembly must call the shared normalized-axis
fingerprinting helper rather than sorting tables locally or catching pandas sort
failures. If the provenance layer cannot prove deterministic axis semantics, it
must fail closed with an actionable provenance exception instead of producing an
apparently normalized, order-sensitive fingerprint.

## Consequences

Code that reads internal stored provenance must treat arrays as immutable
sequences, not as concrete `list` instances. Code that reads documented public
JSON provenance accessors receives mutable-looking fresh `dict` and `list`
payloads, but mutating those copies must not affect the provenance object.

Provenance constructors are now the enforcement boundary for JSON-like
immutability. Workflow assemblers may still normalize domain objects into JSON
facts before constructing provenance, but unsupported values that reach a
provenance dataclass fail at construction rather than being stringified.

Code must not test internal provenance JSON objects with `isinstance(value,
dict)`. Internal provenance mappings satisfy `collections.abc.Mapping`; mutable
`dict` objects are produced only by serialization/thawing helpers or by
documented public JSON accessors. Code that needs an array payload must use one
of those public/thawing paths, because internal provenance arrays are tuples.

This decision complements ADR-0016's dataframe ownership boundary: dataframe
inputs are copied at public boundaries, while JSON-like provenance inputs are
recursively frozen and serialized from fresh payloads.

## Implementation Notes

The recursive freezer lives in `src/phospy/provenance/immutability.py`.
Provenance dataclasses with JSON-like mapping fields route those fields through
the freezer in `__post_init__`. Serialization helpers and selected public JSON
accessors route frozen provenance through the thawing path before returning
payloads.

The same primitive is also the repository-wide policy for exported JSON-like
scientific and result state. Domain owners may adapt freezer failures to their
public boundary exception type, but they must not introduce a second immutable
mapping container or stringify unsupported keys and values.

The concrete container policy is:

- `FrozenJsonMapping` stores a tuple of validated `(str, frozen_value)` pairs
  in private storage and implements `collections.abc.Mapping`;
- `FrozenJsonMapping` has no mutable mapping base class, so `dict.__setitem__`,
  `dict.update`, and related base-class mutation paths cannot operate on it;
- JSON arrays are stored as tuples, so list base-class mutation cannot operate
  on them;
- every nested mapping and sequence is recursively frozen before storage;
- `thaw_json_mapping()` and `thaw_json_value()` always allocate fresh `dict`
  and `list` containers for serialization output; and
- `FrozenJsonMapping.copy()` and dataclass deep-copy helpers thaw into fresh
  ordinary JSON containers; direct internal mapping access remains immutable.

The current audited fields include:

- `ReproducibilityCaveat.details`;
- `TableFingerprint.index_structure` and `column_index_structure`;
- environment dependency, platform, BLAS/LAPACK, thread, locale, and
  constraints mappings;
- preprocessing stage `parameters` and `diagnostics`;
- dataset processing-state missing-data and total-protein correction
  diagnostic JSON mappings;
- batch-correction provenance mappings and rejected-entity `details`;
- reference and Kinase Library resource `sequence_window`, `source_files`, and
  `manifest` mappings;
- run-level `workflow_parameters`;
- scientific-policy `parameters`;
- derived-quantitative missingness, matrix-transformation, and parameter
  mappings;
- trusted dataset construction evidence `details` and construction workflow
  parameters;
- intensity-scale establishment `parameters`;
- quantitative-meaning transition `parameters`, input/output table fingerprint
  payloads, and diagnostic caveat codes;
- public `ResultCaveat.details`;
- importer `ImporterQualityReport.format_specific`,
  `ImporterMissingIntensitySummary` count mappings, and
  `PhosphositeImportResult.diagnostics`;
- enrichment result `diagnostics`, `method_metadata`, `background_summary`, and
  `set_collection_summary`; and
- kinase workflow attrition `metrics`, `policy`, and `policy_violations`
  entries.

`tests/architecture/test_exported_json_state_immutability.py` owns the exported
result-state inventory. New public result dataclass fields typed as mappings
must either enter the immutable JSON policy registry, be typed to a reviewed
immutable domain object, or be added to a narrow follow-up allowlist with a
ticket and reason.
