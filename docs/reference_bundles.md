# Reference Bundles

`ReferenceBundle` is the workflow-ready reference contract for kinase analysis.
Human and mouse are valid organisms, but this release does not ship bundled
human or mouse runtime references. Build those bundles from local source files
with `ReferenceBundleBuilder`.

## Packaged Bundles

Packaged runtime reference data lives only under
`src/phospy/data/reference_bundles/<organism>/<bundle_id>/`.

Current packaged lanes:

| Organism | Bundle ID | Status |
| --- | --- | --- |
| rat | `l6_native` | Bundled in this release |
| human | N/A | Not bundled; no approved redistributable lane is committed |
| mouse | N/A | Not bundled; no approved redistributable lane is committed |

Human or mouse lanes may be added only when the committed manifest documents:

- source name, version, URL, and retrieval method
- license name and license URL
- explicit redistribution status saying redistribution is approved or allowed
- redistribution basis explaining why the license permits packaging
- limitations and supported uses

Do not commit restricted scientific reference datasets. In particular, do not
copy PhosphoSitePlus or Kinase Library data into packaged runtime data unless
their license or a written permission record explicitly allows redistribution.

## Local Builder

```python
from phospy.api import (
    Organism,
    ReferenceBundleBuildRequest,
    ReferenceBundleBuilder,
)

references = ReferenceBundleBuilder().run(
    ReferenceBundleBuildRequest(
        organism=Organism.MOUSE,
        kinase_substrate_path="mouse_kinase_substrates.csv",
        site_sequence_path="mouse_site_sequences.csv",
        source_name="local curated kinase reference",
        source_version="2026-06-11",
        retrieved_at="2026-06-11",
        license="record the source license here",
        redistribution_status="record redistribution status here",
        identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
    )
)
```

The result is a normal validated `ReferenceBundle` and can be passed directly
to `KinaseWorkflowRequest(references=references)`.

## Expected Files

The kinase-substrate file must contain:

- `kinase`
- `substrate_site` or `site_id`

The site-sequence file must contain:

- `site_id` or `substrate_site`
- `site_sequence`

Optional columns are preserved and normalized where present:

- `display_id`
- `gene_symbol` or `gene`
- `protein_accession`
- `protein_id`
- `organism` or `species`

If an `organism` or `species` column is present, every value must match the
requested `Organism`.

## Provenance

Builder-created bundles include local-source provenance:

- source name and version
- retrieval date
- license and redistribution status in the manifest
- identifier namespace
- sequence-window metadata
- source-file paths, SHA-256 digests, and byte counts
- table fingerprints
- identifier-normalisation diagnostics

If `sequence_window` is not supplied, the builder infers it only from uniform
odd-length centered `site_sequence` windows whose central residues are all
`S`, `T`, or `Y`. Otherwise pass an explicit `SequenceWindowDefinition`.

## Failure Policy

The builder fails instead of repairing scientific inputs:

- missing kinase or site columns produce diagnostics listing accepted aliases
  and available columns
- malformed site identifiers are rejected
- duplicate kinase-substrate edges are rejected
- every kinase-substrate site must have a sequence row
- remote URLs are rejected; only local filesystem paths are supported

Do not commit restricted source reference datasets. Keep source files local or
manage them according to the original provider license.
