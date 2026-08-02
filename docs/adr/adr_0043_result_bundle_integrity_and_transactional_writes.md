# ADR-0043: Result Bundle Integrity and Transactional Writes

## Status

- **ADR ID:** ADR-0043
- **Title:** Result Bundle Integrity and Transactional Writes
- **Status:** Accepted
- **Date:** 2026-08-02
- **Decision Type:** Artifact Integrity and I/O Lifecycle Contract

## Context

Reloadable kinase and signalome result bundles persist scientific result tables,
reference tables, provenance, and config snapshots to a directory tree. The
previous bundle manifest recorded logical paths, but the loader trusted those
serialized files after parsing the manifest. A partially written directory or a
later table/config edit could therefore be mistaken for a complete bundle until
table parsing or model validation happened to fail.

Result-bundle integrity must be owned by `io.bundles`. Scientific result models
should continue to represent workflow outputs and provenance, not filesystem
lifecycle, temporary directories, promotion semantics, or byte-level artifact
verification.

## Decision

Kinase and signalome result-bundle manifests use manifest version 2.

Version-2 manifests are content-addressed at the file-entry level:

- each manifest-owned payload file except `manifest.json` records `path`,
  `sha256`, `byte_size`, and `logical_type`;
- table file entries additionally record `shape.rows` and `shape.columns`;
- `manifest.json` remains the trust root and is written last.

Writers serialize bundles into a sibling temporary directory. They write all
tables and JSON sidecars first, hash those exact bytes, build the manifest from
the observed file records, write the manifest last, and then promote the staged
directory to the requested bundle root.

Overwrite is explicit. Existing output directories are rejected unless callers
pass `overwrite=True`. Replacement writes are still staged in a separate
directory; the existing directory is not modified in place.

Loaders parse the manifest, verify all declared file sizes and SHA-256 digests,
and reject undeclared bundle-local files before reconstructing result models.
Digest mismatches, missing files, and stale extra files are fatal.

Existing version-1 result bundles are not compatibility-read. They are rejected
with the existing migration message directing users to regenerate bundles with
the current PhosPy version.

## Consequences

Interrupted new writes cannot leave a target directory with a completed
manifest. If staged serialization fails, the staged directory is cleaned up and
the requested bundle root is not published.

Failed overwrite attempts before promotion leave the previous target directory
unchanged. Successful overwrite attempts replace the previous directory with a
fresh staged tree, preventing stale files from surviving by accident.

The manifest cannot record its own digest without an external trust root or
sidecar signature. This ADR scopes integrity to files addressed by the manifest;
tamper-evident or signed manifest roots remain a separate future decision.

Filesystem lifecycle remains isolated to `phospy.io.bundles._shared` and bundle
writer/loader orchestration. Science and result model constructors receive
already-loaded tables and metadata only.
