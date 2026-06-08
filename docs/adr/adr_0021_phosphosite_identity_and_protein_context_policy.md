# ADR: Phosphosite Identity and Protein Context Policy

## Document Control

- **ADR ID:** ADR-0021
- **Title:** Phosphosite Identity and Protein Context Policy
- **Status:** Superseded
- **Date:** 2026-05-12
- **Decision Type:** Architecture Decision Record

Superseded in part by ADR-0023 and then amended by ADR-0024. ADR-0024 is the
controlling decision for supported analysis-ready row identity.

## Context

PhosPy historically centered site identity on display-oriented `GENE;SITE;`
tokens. That format is convenient and remains useful as a display label, but it
is not sufficient scientific identity for accession-aware phosphoproteomics.

The same display token can correspond to different proteins, isoforms, source
records, or namespaces. If those contexts are silently collapsed, downstream
kinase/signalome/sequence-aware interpretation can become irreproducible or
scientifically ambiguous.

## Decision

PhosPy introduces an explicit phosphosite identity model that separates:

- **Display identity:** standardized `GENE;SITE;` token, now represented as
  `display_id`
- **Scientific identity context:** protein/accession/provenance fields

Identity is represented in a dedicated site-domain module and validated at
dataset/preprocessing boundaries before workflow execution.

ADR-0024 adds the analysis-ready row key: `site_key` is the unique
protein-scoped row identity for `AnalysisReadyPhosphoDataset.phospho.index` and
`AnalysisReadyPhosphoDataset.site_metadata.index`. `display_id` is metadata and
may repeat when distinct `site_key` values preserve the protein context.
ADR-0024 is the current authority for final required identity metadata:
`site_key`, `display_id`, `organism`, `protein_namespace`,
`protein_identifier`, `gene_symbol`, `site`, and `site_sequence`.

### Identity Model

The identity model includes:

- `display_id`
- `gene_symbol`
- `residue`
- `position`
- optional context fields:
  - `organism`
  - `protein_namespace`
  - `protein_identifier`
  - `isoform_id`
  - `source_namespace`
  - `source_site_id`

### Default Site-Token Strictness

By default, phosphosite `site` values must be strict phosphoproteomics tokens:
`S`, `T`, or `Y` followed by a positive integer (for example `S123`, `T45`,
`Y999`).

Malformed or non-phosphorylatable tokens (for example `FOO`, `A123`, `S0`,
blank, or null values) are rejected before workflow execution.

Opaque site values are only valid behind an explicit waiver
(`allow_opaque_site_values=True`), and that waiver must be provenance-visible
in any workflow path that enables it.

### Sequence Context Distinction

`site_sequence` presence is a dataset boundary requirement, but sequence-aware
centred context is a separate execution requirement. For sequence-aware
workflows, centred context requires:

- odd-length sequence windows
- centre residue in `S/T/Y`
- centre residue agreement with the site token residue
- no silent bypass of centre checks via gapped/underscore characters under
  strict mode

If relaxed/gapped/unknown context is allowed, that allowance must be explicit
in the validator policy.

Base dataset/site-metadata validation guarantees only sequence plausibility, not
full centred context semantics. The base `site_sequence` policy is explicit:

- non-empty, non-whitespace strings are required
- allowed amino-acid residues: `A C D E F G H I K L M N P Q R S T V W Y`
- allowed unknown residue: `X`
- allowed gap/placeholders: `_` and `-`
- unsupported letters such as `B`, `J`, `O`, `U`, `Z` are rejected

Base validation deliberately does **not** guarantee odd-length windows,
phosphorylatable centre residue, or centre/site token agreement. Those stricter
rules remain owned by sequence-aware workflow identity contracts.

### Compatibility

Standardized `GENE;SITE;` parsing remains supported for display labels.
Display identity is a reporting surface, not complete scientific identity and
not the analysis-ready row key.

### Ambiguity and Collision Policy

When multiple builder-input rows map to the same display label, PhosPy rejects
the input only when protein-scoped identity cannot be derived safely or when the
rows collapse to the same `site_key` without an explicit duplicate policy that
resolves them before final dataset construction.

Repeated `display_id` values are valid in analysis-ready datasets when each row
has a unique `site_key`.

## Workflow Requirements

Workflow validators must compose shared identity validation through
`src/phospy/validation/workflows/identity.py` and declare one explicit contract.

| Workflow | Contract ID | Identity minimum | Required additions | Explicitly not required |
| --- | --- | --- | --- | --- |
| Differential | `site_key_identity_minimum` | Analysis-ready `site_key` row identity plus `display_id` metadata | Collision checks for unsafe builder-input duplicates when duplicate display labels are present | A one-to-one display-label model |
| Kinase | `sty_site_identity_plus_sequence_context` | Differential minimum | Strict site-token parsing (`S/T/Y<position>` unless opaque waiver), centred sequence context | Treating reference display IDs as analysis-ready row identity |
| Signalome | `protein_scoped_site_identity` | Kinase minimum | Explicit non-empty signalome protein grouping metadata in `protein_id` per retained site | Inference of protein identity or grouping metadata from display IDs |

Reference-organism compatibility remains enforced at the workflow runtime
reference-resolution boundary (`ReferenceCompatibilityValidator` via kinase
workflow interpreter), while workflow identity contracts remain focused on
phosphosite scientific identity semantics.

## Responsibility Audit

Ownership boundaries are:

- **Site identity domain module:** parsing + identity construction + collision
  validation.
- **Dataset/table schemas:** structural shape and column contracts.
- **Workflow identity contract composer:**
  `enforce_workflow_site_identity_contract` in
  `src/phospy/validation/workflows/identity.py`.
- **Workflow validators:** select and apply exactly one workflow identity
  contract; they do not parse identity rows ad hoc.
- **Executors/interpreters:** consume validated identity; they do not infer or
  synthesize identity from fallback heuristics.

## Consequences

### Positive

- Scientific ambiguity in same-display-ID rows is explicit and testable.
- Identity expectations become workflow-aware and auditable.
- Compatibility with legacy display IDs is preserved where acceptable.

### Tradeoffs

- Some duplicate-site datasets that previously collapsed silently now fail fast
  until protein-scoped identity ambiguity is resolved upstream.
- Users may need to provide richer identity metadata for strict workflows.
