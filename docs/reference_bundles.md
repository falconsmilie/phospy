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
| rat | `l6_native` | Approved bundled snapshot with hash-verifiable files and structured exact-file upstream-package license evidence |
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
  `redistribution_status`, `redistribution_allowed`, `redistribution_notes`, and
  `redistribution_evidence`
- package integrity metadata: `table_sha256` and `files`

Every `files` item must contain:

- `relative_path`
- `role`
- `format`
- `sha256`
- `row_count`
- `column_names`

Manifest JSON uses a strict schema. Unsupported extension fields are rejected
rather than ignored at the top level and in every `files[]` entry. Compatibility
aliases emitted by runtime payloads, such as `bundle_id`, `source_files`, or
`sequence_window`, are derived outputs and are not accepted manifest-input
fields.

Sequence-aware manifests may also carry `source_publication`,
`sequence_context_policy`, `sequence_window_length`, `sequence_center_index`,
`allowed_sequence_alphabet`, `organism_common_name`, `supports`, and
`limitations`. When one of `sequence_window_length` or
`sequence_center_index` is present, both must be present and the center index
must be inside the declared window.

`redistribution_status` is the governing redistribution field.
`redistribution_allowed` compatibility value must not be used as the review
authority. If present in JSON, it must be a boolean and must mirror
`redistribution_status`: `true` only for `approved`, `false` for
`external_only` and `unresolved`.
When the raw JSON key `redistribution_allowed` is present, its value must be a
JSON Boolean. JSON `null` is invalid. Omitting the compatibility key is
distinct from supplying `null`.

### Redistribution Status

`redistribution_status` has three allowed values:

- `approved`: the only release-eligible bundled status. Use it only when the
  manifest records verified structured upstream-package license evidence for the
  exact snapshot and exact files being packaged. `source_version`,
  `license_name`, `license_url`, and `redistribution_evidence` are required,
  and `redistribution_notes` must state the evidence basis without contradicting
  the approved state.
- `external_only`: the reference source is known, but users must obtain it
  outside the package under the source provider's terms. External-only
  references must not be shipped as bundled data.
- `unresolved`: redistribution review is incomplete, missing, or not strong
  enough to support packaging. Unresolved bundled references block release.

Codex agents and human developers must not mark references approved without
verified evidence in the manifest. A lineage note, hash, upstream package name,
or optimistic interpretation of a third-party license is not enough.
File hashes establish source and package byte identity for validation; hashes
alone do not establish redistribution approval.

Runtime bundled-reference loading validates the manifest and file hashes before
the tables are exposed to workflows. The release checks enforce stricter
publication rules: every packaged file must be listed, every declared file hash
must match, required organism/source/license metadata must be present,
approved bundled references require structured exact-snapshot redistribution
evidence, and each bundled manifest must declare
`redistribution_status="approved"`. Packaged manifests that declare
`external_only` or `unresolved` fail release validation.

Release-validation diagnostics use stable labels so maintainers and CI logs can
identify the failed manifest and field without reverse-engineering prose. Every
semantic release error includes `reference_id=`, `display_name=`, `organism=`,
`namespace=`, `field=`, `redistribution_status=`, `actual_value=`, and
`reason=`. Source-tree file integrity diagnostics also include `file=`,
`expected_digest=`, and `actual_digest=`; a missing declared file reports
`actual_digest=missing`. These messages are diagnostic output only. External
tools may assert the labels, but must not reconstruct redistribution approval
state from error text.

For `redistribution_status="approved"`, `redistribution_evidence` must be an
object with these enforced fields:

- `evidence_type`: currently `upstream_package_license`
- `upstream_package`: `package_name`, `package_version`, `license_name`, and
  optional `license_url`; `license_name` must be machine-readable and must
  agree with manifest-level `license_name`
- `scope`: `reference_id`, `reference_version`,
  `applies_to_exact_packaged_files`, `packaged_files`, and
  `applies_to_future_bundles`; the scope ID/version must match the manifest,
  exact-file scope must be `true`, future-bundle scope must be `false`, and
  `packaged_files` must be duplicate-free relative POSIX paths that exactly
  equal the manifest `files[].relative_path` set
- `attribution`: `repository_notice_path` and `bundle_attribution_path`; the
  repository notice file must exist during source release validation, and the
  bundle attribution path must be listed in `files` and exist in the bundle
- `independent_database_permission_claimed`: must be `false`
- `verified_at`: optional in the general typed model, but mandatory for
  approved bundled release evidence; the verification date must be supplied
  explicitly as an ISO `YYYY-MM-DD` calendar date and is never inferred from
  retrieval dates, generation dates, filesystem metadata, Git history, notes,
  or wheel metadata
- optional `evidence_url` and narrow-scope `notes`

Unrecognized fields inside `redistribution_evidence` or its nested objects fail
validation. Unsupported extension fields are rejected rather than ignored for
the evidence object, its `upstream_package`, `scope`, and `attribution`
objects.

### Rat `l6_native` Provenance

The packaged rat lane is an exact PhosPy-packaged snapshot derived from PhosR
1.20.0 package data and packaged on 2026-04-16. Its manifest records the
generation lineage as PhosR data objects `phospho.L6.ratio.pe`,
`PhosphoSite.mouse`, and `motif.mouse.list`, generated through
`scripts/active/generate_r_l6_fixtures.R` and redistributed as CSVs under
`src/phospy/data/reference_bundles/rat/l6_native/`.

The upstream PhosR 1.20.0 package metadata declares `GPL-3 + file LICENSE`, and
the rat manifest records typed exact-file license evidence for this committed
PhosR 1.20.0-derived snapshot. That approval applies only to the exact files in
`src/phospy/data/reference_bundles/rat/l6_native/`. It does not claim
independent direct permission from PhosphoSitePlus, PRIDE, Kinase Library, or
other upstream databases, and it must not be generalized to human, mouse,
future PhosR, future PhosPy, or other third-party reference data.

Human or mouse lanes may be added only when the committed manifest documents:

- source name, version if known, URL, and publication where applicable
- license name/text and license URL
- file-level SHA-256 hashes for every packaged reference file
- `redistribution_status="approved"`
- structured exact-snapshot `redistribution_evidence` with upstream package
  license metadata, exact packaged-file scope, attribution paths, no
  future-bundle scope, and no independent database permission claim
- redistribution notes explaining why the license or permission permits
  packaging without generalizing approval to other bundles or external datasets
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
        source_version="upstream-package-7.4",
        retrieved_at="2026-07-14",
        license="record the source license here",
        redistribution_status="record redistribution status here",
        identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        reference_version="local-snapshot-2026-07-14",
    )
)
```

The result is a normal validated `ReferenceBundle` and can be passed directly
to `KinaseWorkflowRequest(references=references)`.

`ReferenceBundleBuildRequest.source_version` identifies the upstream package,
database, or caller-source version. `reference_version` identifies the local
PhosPy snapshot being built. For example, a local bundle can use
`source_version="upstream-package-7.4"` and
`reference_version="local-snapshot-2026-07-14"` at the same time.

If `reference_version` is omitted, the builder derives a deterministic local
snapshot version from the two source-file SHA-256 fingerprints:

```text
kinase_substrate:<kinase SHA-256>\nsite_sequences:<sequence SHA-256>\n
```

The builder SHA-256 hashes that specified ASCII byte string and emits
`local-snapshot-sha256-<64-character digest>`. Rebuilding identical files
produces the same local snapshot version; changing either source file changes
it.

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

`source_version` is the upstream package, database, or caller-supplied source
identity version. `reference_version` is the local PhosPy reference snapshot
version for the packaged or constructed bundle. These values are intentionally
separate: a missing upstream source identity is unknown, not implicitly equal to
the local snapshot version, and provenance must not copy `reference_version`
into `source_version`. For builder-created local bundles, omit
`reference_version` only when the content-derived
`local-snapshot-sha256-...` identity is the intended local snapshot version.

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

- `KinaseScoringConfig(reliability_profile=KinaseReliabilityProfile.CUSTOM, scoring_mode="kinase_library_motif")`
- `KinaseScoringConfig(reliability_profile=KinaseReliabilityProfile.CUSTOM, scoring_mode="combined_profile_motif")`

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
