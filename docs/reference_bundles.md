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
| rat | `l6_native` | Bundled locally with a hash-verifiable manifest and explicit provenance caveats |
| human | N/A | Not bundled; no approved redistributable lane is committed |
| mouse | N/A | Not bundled; no approved redistributable lane is committed |

Each packaged bundle must have one manifest for the logical reference bundle,
not one manifest per CSV.

### Bundled Manifest Schema

Every bundled reference manifest must contain these top-level fields:

- identity and scope: `reference_id`, `display_name`, `organism`,
  `taxonomy_id`, `protein_namespace`, and `reference_version`
- source and build provenance: `source_name`, `source_url`, `source_version`,
  `retrieved_at`, `derived_from`, `generated_by`, `generated_at_utc`, and
  `manifest_schema_version`
- license and redistribution metadata: `license_name`, `license_url`,
  `redistribution_status`, and `redistribution_notes`
- package integrity metadata: `table_sha256` and `files`

Every `files` item must contain:

- `relative_path`
- `role`
- `format`
- `sha256`
- `row_count`
- `column_names`

Sequence-aware manifests may also carry `source_publication`,
`sequence_context_policy`, `sequence_window_length`, `sequence_center_index`,
`allowed_sequence_alphabet`, `organism_common_name`, `supports`, and
`limitations`. When one of `sequence_window_length` or
`sequence_center_index` is present, both must be present and the center index
must be inside the declared window.

`redistribution_allowed` is a derived compatibility boolean, not the
authoritative redistribution field for review or release. The authoritative
field is `redistribution_status`.

### Redistribution Status

`redistribution_status` has three allowed values:

- `approved`: the only release-eligible bundled status. Use it only when the
  manifest records verified license or permission evidence for the exact files
  being packaged. `license_name` and `license_url` are required, and
  `redistribution_notes` must state the evidence basis.
- `external_only`: the reference source is known, but users must obtain it
  outside the package under the source provider's terms. External-only
  references must not be shipped as bundled data.
- `unresolved`: redistribution review is incomplete, missing, or not strong
  enough to support packaging. Unresolved bundled references block release.

Codex agents and human developers must not mark references approved without
verified evidence in the manifest. A lineage note, hash, upstream package name,
or optimistic interpretation of a third-party license is not enough.

Runtime bundled-reference loading validates the manifest and file hashes before
the tables are exposed to workflows. The release gate enforces stricter
publication rules: every packaged file must be listed, every declared file hash
must match, required organism/source/license metadata must be present, and each
bundled manifest must declare `redistribution_status="approved"`. Packaged
manifests that declare `external_only` or `unresolved` fail release validation.

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
accession notes. The rat manifest records explicit caveats for this exact
derived CSV snapshot. Do not use this bundle as approval precedent for human,
mouse, PhosphoSitePlus, Kinase Library, PRIDE, or other third-party reference
data.

Human or mouse lanes may be added only when the committed manifest documents:

- source name, version if known, URL, and publication where applicable
- license name/text and license URL
- file-level SHA-256 hashes for every packaged reference file
- `redistribution_status="approved"`
- redistribution notes explaining why the license or permission permits
  packaging
- sequence-context policy when sequence windows are included
- limitations and supported uses in the manifest payload

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

## Validation Report

Every constructed `ReferenceBundle` has a structured validation report:

```python
report = references.validation_report

print(report.bundle_name)
print(report.kinase_substrate_record_count)
print(report.compatibility_warnings)
```

The report is informational. It does not repair reference data and does not make
invalid bundles usable. Missing required tables, missing required columns,
duplicate kinase-substrate records, malformed identifiers, and missing sequence
rows still fail validation. Reference validators and manifest-validation
internals are implementation routes, not `phospy.api` public exports.

Report fields include:

- bundle name and version, when manifest or provenance metadata provides them
- organism and identifier namespace metadata, when available
- required table status and required columns
- required source-file metadata from the manifest, when available
- kinase-substrate record count
- duplicate-record count for accepted tables
- missing-value counts for important fields
- available provenance fields such as source name, source version, license, and
  redistribution status
- compatibility warnings for limited metadata

Kinase workflows keep the validated `ReferenceBundle` on the result, so you can
inspect the same report before or after a run:

```python
result = KinaseWorkflow().run(request)
report = result.references.validation_report
```

Use the report to check whether the reference source, organism, namespace, and
provenance are suitable for your dataset before interpreting kinase scores. It
does not change the scientific meaning of the workflow outputs.

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
- source-file paths, SHA-256 digests, row counts, and column names when
  available
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
matrix resources. These resources can be used by opt-in `KinaseWorkflow`
scoring modes through `KinaseWorkflowRequest.kinase_library_resource`:

- `KinaseScoringConfig(scoring_mode="kinase_library_motif")`
- `KinaseScoringConfig(scoring_mode="combined_profile_motif")`

This is support for caller-supplied local resources with explicit provenance.
PhosPy does not bundle official Kinase Library data and does not claim official
Kinase Library parity.

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

Scores are preserved as numeric provider-scale values in the loaded resource.
Workflow motif scores are normalized support scores for within-run ranking; they
are not calibrated probabilities and do not replace kinase-substrate reference
membership or activity inference.

Loaded resources include explicit provenance:

- source name and version
- license
- score scale
- organism applicability
- sequence-window definition
- source-file path, SHA-256 digest, and byte count
- matrix table fingerprints
