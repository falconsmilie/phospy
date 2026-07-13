# ADR: Reference and Fixture Data Policy for PhosPy

## Document Control

- **ADR ID:** ADR-0015
- **Title:** Reference and Fixture Data Policy for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how reference data and test fixture data should be managed in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, the package must distinguish clearly between bundled runtime reference data, internal development fixtures, and parity reference outputs used for scientific comparison.

The decision is to adopt an explicit data policy with three distinct categories:

- bundled runtime reference data
- non-runtime test fixtures
- parity reference outputs

Each category should have a clear purpose, location, update rule, and authority level. Runtime reference data exists to support the package's actual scientific workflows. Test fixtures exist to support validation of behaviour. Parity reference outputs exist to support explicit scientific comparison, primarily against PhosR outputs. These categories must not be mixed casually.

## Status

Accepted.

This ADR defines the data policy that supports the reference-resolution
architecture in ADR-0004, the parity strategy in ADR-0013, the test structure
in ADR-0014, and the broader fresh-start rewrite defined in ADR-0012.

## Context and Problem Statement

Earlier ADRs already depend on several kinds of data artefacts:

- bundled reference resources for supported organisms
- builder and workflow test data
- parity fixtures and expected scientific outputs

Without a clear policy for these artefacts, the package risks several kinds of structural and scientific confusion:

- runtime reference data becoming mixed with test-only data
- parity outputs being treated as normal runtime resources
- unclear authority about which data source defines expected behaviour
- stale or duplicated fixture sets
- accidental shipping of unnecessary data in the package
- unclear update rules when bundled references or parity expectations change

PhosPy therefore needs an explicit decision about how reference and fixture data are classified, stored, updated, and treated as authoritative.

## Decision Drivers

The decision is driven by the following considerations:

1. **Scientific trust.** The package must make it clear what runtime reference data it depends on and why.
2. **Boundary clarity.** Runtime data and test data serve different purposes and should be treated differently.
3. **Maintainability.** Data artefacts need clear ownership and update expectations.
4. **Packaging discipline.** Only data required for supported runtime behaviour should ship as bundled package data.
5. **Parity integrity.** Scientific comparison outputs should remain explicit and reviewable.
6. **Fresh-start discipline.** Data policy should support the new architecture rather than preserve accidental patterns from the previous application.

## Decision

PhosPy will classify reference and fixture data into three explicit categories:

1. **Bundled runtime reference data**
2. **Test fixture data**
3. **Parity reference outputs**

Each category has a different role:

- bundled runtime reference data supports real package behaviour
- test fixture data supports ordinary testing
- parity reference outputs support explicit scientific comparison

These categories should be physically and conceptually separated.

## Core Design Principle

Data artefacts should be organised by **purpose** and **authority**, not simply by file type or historical origin.

The same CSV-like file format may appear in multiple places, but the package must still distinguish whether that file is:

- a bundled runtime resource
- a normal test fixture
- a scientific reference expectation for parity

## Category 1: Bundled Runtime Reference Data

Bundled runtime reference data is the data that PhosPy ships with the package to support real workflow execution.

Examples include:

- organism-specific bundled reference resources
- supported kinase-substrate mapping resources
- bundled sequence-oriented resources required for supported reference resolution

### Policy

Bundled runtime reference data should:

- live with the runtime reference system, not in test directories
- be treated as package resources
- be versioned as part of the package source
- be used only for supported runtime scenarios
- remain as small and explicit as practical
- include a machine-readable manifest (reference data card) per bundled runtime lane

### Authority

Bundled runtime reference data is authoritative for PhosPy's bundled reference behaviour.

It is not automatically authoritative for parity expectations unless explicitly used as part of a parity fixture definition.

Each runtime manifest should declare organism scope, identifier namespace, source metadata, license or redistribution status, supported workflows, sequence-window definition, and known limitations.

### Runtime Reference Data Card Requirements

Every bundled runtime lane must include a machine-readable manifest (reference
data card) that is treated as runtime governance metadata, not optional notes.
Every bundled reference manifest must contain the schema fields required by
`src/phospy/science/references/validation.py`:

- identity and scope: `reference_id`, `display_name`, `organism`,
  `taxonomy_id`, `protein_namespace`, and `reference_version`
- source and build provenance: `source_name`, `source_url`, `source_version`,
  `retrieved_at`, `derived_from`, `generated_by`, `generated_at_utc`, and
  `manifest_schema_version`
- license and redistribution metadata: `license_name`, `license_url`,
  `redistribution_status`, `redistribution_allowed`,
  `redistribution_notes`, and `redistribution_evidence`
- package integrity metadata: `table_sha256` and `files`

Every `files` entry must contain `relative_path`, `role`, `format`, `sha256`,
`row_count`, and `column_names`. Sequence-aware bundled references should also
record the sequence-context policy, sequence-window length, center index,
allowed alphabet, supported uses, and known limitations when applicable.

`redistribution_status` is the governing redistribution field. The derived
`redistribution_allowed` compatibility value must not be used as the review
authority. If JSON includes `redistribution_allowed`, it must be a boolean and
must mirror `redistribution_status`: `true` for `approved`, `false` for
`external_only` or `unresolved`.

### Redistribution Status Policy

The allowed `redistribution_status` values are:

- `approved`: release eligible only when the manifest contains typed
  upstream-package license evidence for the exact snapshot and exact files being
  packaged.
- `external_only`: reference source is known, but the package must not ship the
  data. Users may supply local files under their own source terms.
- `unresolved`: redistribution evidence is missing, incomplete, or not yet
  reviewed.

Unresolved bundled references block release. External-only references must not
be shipped as bundled data. Codex agents and human developers must not mark
references approved without verified evidence. Source lineage, file hashes,
and upstream package names are provenance evidence, not redistribution approval
by themselves.

For approved bundled manifests, `redistribution_evidence` is release-significant
typed state rather than approval prose. It records:

- `evidence_type="upstream_package_license"`
- upstream package name, package version, machine-readable license name, and
  optional license URL
- scope reference ID/version, exact packaged-file coverage, duplicate-free
  relative POSIX packaged-file paths, and no future-bundle coverage
- repository notice and bundle-local attribution paths
- `independent_database_permission_claimed=false`

The release gate checks that evidence license metadata agrees with the manifest
license, scope ID/version match the manifest, scope files exactly equal
`files[].relative_path`, attribution files exist, and contradictory legal
approval text or affirmative independent database permission claims are absent.

## Category 2: Test Fixture Data

Test fixture data is data used to exercise behaviour in unit and integration tests.

Examples include:

- small input tables for builder tests
- small reference-like tables for validator tests
- small transformed datasets for dataset model tests
- workflow test inputs for integration scenarios

### Policy

Test fixture data should:

- live under the test tree
- be small and focused where possible
- serve behavioural verification rather than shipping behaviour
- be shareable across integration and parity suites where this improves consistency

Fixture sharing between integration and parity suites is acceptable where it keeps both test classes grounded on the same meaningful input data.

### Authority

Test fixture data is authoritative only for the tests that intentionally use it.

It is not runtime reference data and should not be used as though it defines supported package behaviour outside the test suite.

## Category 3: Parity Reference Outputs

Parity reference outputs are expected scientific outputs used to validate scientific parity, primarily against PhosR-derived expectations.

Examples include:

- expected scoring outputs
- expected prediction outputs
- expected activity outputs
- expected signalome outputs
- analysis-ready dataset expectations where parity is defined at that boundary

### Policy

Parity reference outputs should:

- live under the parity test area or a parity-specific fixture area
- be treated as high-value scientific assets
- remain clearly separate from normal runtime reference data
- have explicit comparison rules associated with them
- be kept current and correct

Parity fixtures do not require extra repository-level versioning or annotation machinery beyond staying current and scientifically appropriate.

Stale parity fixtures should not exist. Removed or superseded parity fixtures should be deleted rather than retained in a way that creates confusion.

### Authority

The primary authority for parity reference outputs should be PhosR outputs where available and appropriate.

Selected outputs from the old application may still be used as secondary reference material where scientifically helpful, but they are not the primary authority.

## Non-Interchangeability Rules

The three categories are intentionally non-interchangeable:

1. Runtime reference loaders must not read test fixture or parity-output trees.
2. Test fixtures must not be promoted to runtime references without data-card
   review.
3. Parity outputs must never be treated as runtime lookup resources.

## Physical Layout Direction

A likely healthy direction is:

```text
phospy/
  references/
    resources/

  data/

tests/
  support/
  unit/
  integration/
  parity/
    fixtures/
```

This captures the intended separation while keeping future runtime data needs open-ended:

- bundled runtime reference data may live under the runtime reference system and a broader runtime data area where justified later
- ordinary fixtures live with tests
- parity outputs live with the parity suite

Runtime resolution code must not treat files under `tests/fixtures/` as bundled runtime references.

A broader runtime `data/` area is acceptable if future package needs go beyond the current reference-resource model.

## Shipping Policy

Only bundled runtime reference data should ship as part of the runtime package by default.

Ordinary test fixtures and parity reference outputs are not runtime resources and should not be packaged as though they are part of the runtime reference system.

This prevents the runtime package from silently accumulating unnecessary scientific comparison artefacts or bulky test data.

## Update Policy

### Bundled Runtime Reference Data

Update when supported runtime reference behaviour changes intentionally.

### Test Fixture Data

Update when test scenarios change or when fixture quality needs improvement.

### Parity Reference Outputs

Update only with deliberate review, because changing these artefacts changes what the project considers scientifically equivalent or acceptable.

Parity outputs must not drift casually.

Realism is the higher priority for parity assets. Large parity fixtures should not be minimised so aggressively that they stop representing meaningful real scientific scenarios.

## Source-Of-Truth Policy

The source of truth depends on the category.

### For Runtime Reference Behaviour

The source of truth is the bundled runtime reference data shipped by the package for supported built-in reference resolution.

### For Ordinary Behavioural Tests

The source of truth is the explicitly defined fixture or scenario in that test.

### For Scientific Parity

The source of truth is the parity reference output, whose primary source should be PhosR where available and appropriate.

The old application may inform parity, but it must not silently become the default authority.

## Review Policy for Data Changes

Changes to bundled runtime reference data or parity reference outputs should be treated as meaningful review events.

A reviewer should be able to answer:

- what data changed?
- why did it change?
- what category of data is this?
- does this alter runtime behaviour, scientific expectation, or only a local test scenario?

This does not require a heavy process, but it does require deliberate review.

Bundled runtime reference updates and parity-output updates should be called out explicitly in changelog or release notes when they materially affect scientific behaviour.

## Relationship to Reference Resolution

ADR-0004 defines how references are resolved at runtime.

This ADR defines how the data underpinning that system is classified and managed.

In particular:

- bundled runtime reference data supports `ReferencePreset` resolution
- user-supplied `ReferenceBundle` remains outside bundled data policy as caller-supplied input
- parity outputs must not be mistaken for runtime bundled references

## Relationship to Parity Policy

ADR-0013 defines what parity means.

This ADR defines how the data used for parity should be handled.

In particular:

- parity outputs are explicit scientific comparison assets
- they are not runtime reference resources
- they should primarily be derived from PhosR outputs

## Loader Policy

Common test data loaders may live wherever they are easiest to use and where the ownership remains clear.

A shared location such as `tests/support/` is appropriate when loaders are genuinely reusable. Parity-specific loaders may stay close to parity fixtures when that keeps usage and maintenance simpler.

The governing rule is practicality and clarity, not forced uniformity.

## Consequences

### Positive Consequences

- Runtime data and test data become easier to reason about.
- The package avoids mixing reference resources with parity artefacts.
- Scientific comparison assets remain explicit and reviewable.
- Shipping policy becomes clearer and tighter.
- Data changes become easier to classify during review.

### Negative Consequences

- The project must stay disciplined about where new data artefacts are placed.
- Some existing fixtures or resources may need to be moved to fit the policy.
- Parity artefacts require deliberate maintenance rather than casual accumulation.

### Neutral Consequences

- Some test fixtures may still resemble runtime reference data structurally, even though they serve a different purpose.
- Internal helper utilities may still load data across categories when tests intentionally require it, as long as the category boundary remains clear.

## Rejected Alternatives

### Alternative 1: Keep All Scientific Data Files Together Regardless of Purpose

This option was rejected because it weakens the distinction between runtime reference data, ordinary fixtures, and parity expectations.

### Alternative 2: Treat Old-Application Outputs as the Default Authority for Parity Fixtures

This option was rejected because the project is a PhosR port and the primary parity target should remain PhosR outputs.

### Alternative 3: Ship Parity Outputs with the Runtime Package by Default

This option was rejected because parity outputs are test assets, not runtime reference resources.

### Alternative 4: Create Heavy Repository-Level Metadata/Versioning Rules for Parity Fixtures Immediately

This option was rejected because the current direction is to keep the system simple and rely on keeping parity fixtures up to date rather than adding extra management machinery too early.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. A broader runtime `data/` area is acceptable for future package needs beyond the current `references/resources/` model.
2. Realism is the higher priority for large parity fixtures.
3. Common test data loaders may live wherever they are easiest to use and maintain, provided ownership remains clear.
4. Removed or superseded parity fixtures should be deleted so stale fixtures do not create confusion.
5. Bundled runtime reference updates and parity-output updates should be called out explicitly in changelog or release notes when they materially affect scientific behaviour.

## Implementation Guidance

A likely healthy direction is:

- keep bundled runtime reference data with the runtime reference system and broader runtime data area where justified
- keep ordinary test fixtures in the test tree
- keep parity outputs in the parity suite
- keep PhosR outputs as the primary reference source for parity where possible
- review data changes according to purpose, not just file format
- remove stale parity artefacts rather than carrying them forward

Reviewers should reject data placement or reuse patterns that blur the distinction between runtime references, ordinary fixtures, and parity expectations.

## Scope Boundaries

This ADR defines reference and fixture data policy only.

It does not define:

- the scientific algorithms themselves
- the full test execution policy beyond what other ADRs already establish
- release-note formatting rules
- external data acquisition tooling
- future benchmark dataset policy

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check future changes against the following questions:

1. What category of data is this artefact?
2. Is it stored in the right place for that category?
3. Is its authority level clear?
4. Should it ship with the runtime package, or stay test-only?
5. Does this change strengthen or weaken scientific trust and maintainability?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-0004 defines the reference resolution strategy and `ReferenceBundle`
  contract.
- ADR-0008 defines the builder architecture below the dataset boundary.
- ADR-0013 defines scientific parity strategy and parity-testing policy.
- ADR-0014 defines the overall test suite structure and policy.
- ADR-0015 defines how reference data and fixture data should be classified and
  managed.
- ADR-0016, ADR-0017, ADR-0018, and ADR-0019 introduce feature lanes that
  require explicit fixture-vs-runtime governance and category boundaries.

Together, these ADRs establish:

- one explicit runtime reference system
- one clear distinction between runtime data and test data
- one parity policy driven primarily by PhosR outputs
- one disciplined policy for where scientific data artefacts belong

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub. https://github.com/PYangLab/PhosR
