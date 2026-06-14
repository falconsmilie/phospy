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
| rat | `l6_native` | Bundled in this release with explicit provenance and redistribution caveats |
| human | N/A | Not bundled; no approved redistributable lane is committed |
| mouse | N/A | Not bundled; no approved redistributable lane is committed |

### Rat `l6_native` Provenance

The packaged rat lane is a PhosR-derived snapshot packaged on 2026-04-16.
Its manifest records the generation lineage as PhosR data objects
`phospho.L6.ratio.pe`, `PhosphoSite.mouse`, and `motif.mouse.list`, generated
through `scripts/active/generate_r_l6_fixtures.R` and redistributed as CSVs
under `src/phospy/data/reference_bundles/rat/l6_native/`.

The upstream PhosR package metadata declares GPL-3 + file `LICENSE`. That is
not the same as independent approval for every underlying scientific source:
PhosR documentation identifies `PhosphoSite.mouse` as extracted from
PhosphoSitePlus and identifies the L6 phosphoproteome object with PRIDE
accession notes. The rat manifest therefore records redistribution for this
exact derived CSV snapshot as not independently verified. Do not use this
bundle as approval precedent for human, mouse, PhosphoSitePlus, Kinase Library,
PRIDE, or other third-party reference data.

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

## Kinase Library-Style Matrix Resources

PhosPy also has schema and loader support for local Kinase Library-style motif
matrix resources. This is reference-resource support only: these matrices are
not wired into `KinaseWorkflow`, and existing kinase scoring output is
unchanged.

Use `KinaseLibraryResourceLoader` for local files:

```python
from phospy.api import KinaseLibraryResourceLoader

resource = KinaseLibraryResourceLoader().run("kinase_library_matrices.csv")
```

The loader accepts local CSV, TSV, or parquet tables. A long table must contain:

- `kinase`
- `residue_class` with values `ser_thr` or `tyr`
- `position` as integer positions relative to the phospho-acceptor residue
- `amino_acid`
- `score`

It also requires resource metadata, either as repeated table columns or as a
`KinaseLibraryResourceLoadRequest`:

- `source_name`
- `source_version`
- `license`
- `score_scale`
- `organisms`
- `upstream_residues`
- `downstream_residues`
- `central_residue_required`

Optional kinase metadata columns are preserved:

- `kinase_family`
- `kinase_group`

Scores are preserved as numeric provider-scale values. PhosPy does not treat
Kinase Library scores as probabilities, does not normalize them into the
internal motif-frequency scorer, and does not force them into the
kinase-substrate map shape.

Loaded resources include explicit provenance:

- source name and version
- license
- score scale
- organism applicability
- sequence-window definition
- source-file path, SHA-256 digest, and byte count
- matrix table fingerprints
