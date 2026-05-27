# ADR: Phosphosite Identity and Protein Context Policy

## Document Control

- **ADR ID:** ADR-0021
- **Title:** Phosphosite Identity and Protein Context Policy
- **Status:** Superseded
- **Date:** 2026-05-12
- **Decision Type:** Architecture Decision Record

Superseded in part by ADR-0023 and ADR-0024 for supported analysis-ready row
identity scope and duplicate display-site handling.

## Context

PhosPy historically centered site identity on display-oriented `GENE;SITE;`
tokens. That format is convenient and remains useful for matrix indexing, but
it is not sufficient scientific identity for accession-aware phosphoproteomics.

The same display token can correspond to different proteins, isoforms, source
records, or namespaces. If those contexts are silently collapsed, downstream
kinase/signalome/sequence-aware interpretation can become irreproducible or
scientifically ambiguous.

## Decision

PhosPy introduces an explicit phosphosite identity model that separates:

- **Display identity:** standardized `GENE;SITE;` token
- **Scientific identity context:** protein/accession/provenance fields

Identity is represented in a dedicated site-domain module and validated at
dataset/preprocessing boundaries before workflow execution.

### Identity Model

The identity model includes:

- `display_id`
- `gene_symbol`
- `residue`
- `position`
- optional context fields:
  - `organism`
  - `protein_id`
  - `protein_accession`
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

Standardized `GENE;SITE;` parsing remains supported and is still required for
display IDs at strict dataset boundaries. Display identity remains an indexing
surface, not complete scientific identity.

### Ambiguity and Collision Policy

When multiple rows map to the same display site ID, PhosPy now rejects rows
that carry conflicting scientific identity contexts (for example conflicting
`protein_id`/`protein_accession` values) before duplicate-site collapse.

Semantically identical repeated identities remain valid.

## Workflow Requirements

Workflow validators must compose shared identity validation through
`src/phospy/validation/workflows/identity.py` and declare one explicit contract.

| Workflow | Contract ID | Identity minimum | Required additions | Explicitly not required |
| --- | --- | --- | --- | --- |
| Differential | `display_site_identity_minimum` | Standardized display IDs and coherent `gene_symbol/site` rows | Collision checks for conflicting scientific context when duplicate display IDs are present | Protein/accession fields |
| Kinase | `sty_site_identity_plus_sequence_context` | Differential minimum | Strict site-token parsing (`S/T/Y<position>` unless opaque waiver), centred sequence context | Mandatory protein/accession on every row |
| Signalome | `protein_scoped_site_identity` | Kinase minimum | Explicit non-empty `protein_id` per retained site | Inference of protein identity from display IDs |

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
  until identity ambiguity is resolved upstream.
- Users may need to provide richer identity metadata for strict workflows.
