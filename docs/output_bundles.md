# Output Bundles

PhosPy can write workflow tables to an ordinary directory for inspection,
sharing, or use by another tool. These public publisher helpers are the normal
choice for user code.

PhosPy also maintains reloadable bundle formats for durable-state and release
verification. Those codecs are internal. They are documented later on this
page so that the persistence guarantees and limitations remain explicit.

## Choose What to Save

Use a publisher when you need a supported directory of result tables:

| Result in Memory | Publisher |
| --- | --- |
| `AnalysisReadyPhosphoDataset` | `publish_dataset(...)` |
| `KinaseWorkflowResult` | `publish_kinase_workflow(...)` |
| `SignalomeWorkflowResult` | `publish_signalome_workflow(...)` |

Publisher directories are table publications. They are not reloadable Python
result objects.

## Publish Workflow Tables

Publisher helpers are part of the advanced supported API. Import the helper
for the result you want to write:

```python
from pathlib import Path

from phospy.advanced import publish_kinase_workflow

output_directory = Path("kinase-output")

written_files = publish_kinase_workflow(
    kinase_result,
    output_directory,
    output_format="csv",
)
```

Supported table formats are `csv`, `tsv`, and `parquet`. Install the optional
Parquet dependency before selecting `parquet`.

The returned value identifies the files written by the publisher. Keep the
directory together because its manifest describes the publication.

## Protect Existing Output

Publisher helpers reject an existing output root by default. Existing files,
directories, and symlinks at that root are treated as conflicts before any
output file is written:

```python
publish_kinase_workflow(kinase_result, output_directory)
```

Pass `overwrite=True` only when you want to replace the entire publication
root:

```python
publish_kinase_workflow(
    kinase_result,
    output_directory,
    overwrite=True,
)
```

Each publisher stages the complete output tree in a hidden sibling directory on
the same filesystem as the final path. The final publication root is created or
replaced only after all tables and manifests are written in staging. On
successful overwrite, the previous root is replaced as a whole, so stale files
from an older publication cannot survive.

When `overwrite=True`, PhosPy moves the existing root to a sibling recovery
backup, promotes the staged replacement, then removes the backup. If promotion
fails, PhosPy restores the original root when the filesystem allows it. If a
cleanup step fails after a safe promotion, the raised error reports the retained
backup path instead of silently deleting the only recoverable copy.

Filesystem atomicity is limited by the local filesystem and operating system.
The transaction uses sibling directory moves to avoid cross-device promotion,
but it is not a lock service for concurrent writers and is not designed for
object stores or network filesystems with non-standard rename semantics.

## Understand the Publisher Layout

`publish_dataset(...)` writes:

```text
dataset/
  phospho.*
  site_metadata.*
  sample_metadata.*   # optional
  total.*             # optional
  manifest.json
```

`publish_kinase_workflow(...)` also writes:

```text
kinase/
  scoring/
    profile_scores.*
    motif_scores.*                    # optional
    rank_weighted_fusion_scores.*     # optional
    score_fusion_weights.*            # optional
    substrate_contributions.*         # optional
  prediction/
    pred_mat.*
    substrate_list.*                  # optional
  activity/                           # optional
  references/
    kinase_substrate_map.*
    site_sequences.*
  manifest.json
```

`publish_signalome_workflow(...)` also writes:

```text
signalome/
  module_assignments.*
  signalome_modules.*
  kinase_network_edges.*
  kinase_network_nodes.*                  # optional
  kinase_network_candidate_correlations.* # optional
  expanded_signalome.*                    # optional
  manifest.json
```

Site-level outputs preserve the enforced identity contract. Dataset
`site_metadata` is `site_key` indexed and includes both `site_key` and
`display_id`. Kinase and signalome site-level tables that materialize row
identity also include both columns where applicable; internal workflow alignment
remains `site_key` based.

Publisher outputs and reloadable bundles are generated artefacts. They are useful
analysis outputs, but they are not source archives and should not be packaged
into clean source/release archives. Build source archives from tracked source
state instead, excluding generated bundle directories, `site/`, `build/`,
`dist/`, cache directories, and prior archive files.

## Understand Reloadable Bundles

Kinase and signalome owners maintain reloadable bundle codecs for durable-state
and release verification. Those codecs and their config-snapshot types are
owner-internal implementation surfaces, not supported user import routes. The
supported aggregate output API is the `publish_*` family shown above.

A saved bundle contains:

- data tables
- reference tables
- workflow output tables
- `manifest.json`
- `config/snapshot.json`

Loaders return `LoadedKinaseWorkflowBundle` or `LoadedSignalomeWorkflowBundle`.
Each object contains the reconstructed result, the config snapshot, and the
manifest version.

## Bundle Integrity and Limits

### Verify Bundle Integrity

Result bundle manifests are versioned. Current kinase bundle writers emit
manifest version 3. Current signalome bundle writers emit manifest version 2.

Bundle-relative manifest paths must stay inside the bundle root. Absolute paths
and paths that escape the root are rejected.

The manifest is the bundle trust root and is written last. It records each
bundle-owned payload file's relative path, SHA-256 digest, byte size, and logical
type. Table entries also record row and column counts. Loaders verify declared
file sizes and digests, reject missing files, and reject undeclared stale files
before reconstructing workflow result models.

Writers stage bundles in a sibling temporary directory and promote the completed
directory only after all tables, JSON sidecars, and the manifest have been
written. Existing output directories are rejected by default. Pass
`overwrite=True` to replace an existing bundle with a freshly staged bundle.

Optional tables remain optional. Missing optional activity, motif, node,
site-membership, or context tables are represented as `None` when loaded.

### Preserve Activity Semantics

The manifest records output format, bundle kind, manifest version, table file
entries, reference organism, and provenance where available. Dataset metadata in the
manifest includes `processing_state`, so explicit quantitative meaning
(`phosphosite_log_abundance`, `phospho_total_log_ratio`, or mixed state) is
preserved in published outputs.

For kinase activity score outputs, manifest metadata includes explicit activity
method identity (`activity_method_id`, family, and non-KSEA/non-PhosR-equivalence
flags) when activity is enabled. Kinase manifest version 3 also persists the
typed `ActivityInputSemantics` and `ActivityProfileMetadata` payloads from
`KinaseActivityResult`. These payloads preserve the declared profile axis,
quantitative semantics, profile identifiers, sample/condition/contrast
identifiers, condition-summary aggregation metadata, and activity-matrix column
axis semantics. Loaders reconstruct current-schema activity results from these
typed payloads; they do not infer scientific activity semantics from activity
method names, table labels, config strings, or diagnostic text. When provenance
also records resolved activity profile-axis or quantitative-semantics values,
the loader checks that provenance agrees with the manifest result semantics.

Within `outputs.activity`, schema version 3 uses the fields `enabled`, `method`,
`summary`, `input_semantics`, `profile_metadata`, `membership_selection`, and
`tables`. Enabled KSEA activity persists typed membership-selection provenance
and its inferential eligibility state. Enabled activity requires
`input_semantics` and `profile_metadata` objects. Disabled activity requires
those fields and `membership_selection` to be `null`.

KSEA-style and ssGSEA-style runs also emit an activity `statistics_table` with
method statistics, p-values when available, optional q-values, substrate counts,
background counts, and computability statuses. The statistics-table row
identity is `profile_id`.
For KSEA-style runs whose substrate membership consumed the tested quantitative
matrix, `p_value_matrix` and `q_value_matrix` are absent and `statistics_table`
uses missing `p_value`/`q_value` cells with structured inferential-status
columns. Benjamini-Hochberg adjustment is run only over eligible finite p-values.
For `ActivityProfileAxis.SAMPLE`, `profile_id` contains sample IDs; for
`ActivityProfileAxis.CONDITION_SUMMARY`, it contains condition-summary profile
IDs and an optional `condition` column may also be present with the same values;
for `ActivityProfileAxis.CONTRAST`, it contains contrast IDs; and for
`ActivityProfileAxis.EFFECT`, it contains neutral effect profile IDs. Sample,
contrast, and effect statistics tables do not persist a `condition` column.
These outputs are exploratory substrate-supported scores, not direct causal
kinase activation evidence.

### Preserve Quantitative State

Bundle loaders reconstruct `dataset.processing_state` and
`dataset.intensity_scale_state` from saved payloads, including mixed
corrected/uncorrected total-protein quantitative meaning and row-level
total-protein correction diagnostics.

### Preserve Caveats and Provenance

Structured result caveats are persistence data, not display-only warnings.
Current reloadable bundle services persist caveats for kinase and signalome
results, and `DifferentialAnalysisResult.to_payload()` serializes differential
caveats, differential policy provenance, workflow provenance, and
feature-eligibility tables for handoff. Any future differential bundle writer
must preserve those payloads, including exploratory
single-biological-replicate caveats that distinguish computable output from
production-supported inferential support and imputation inference fields that
distinguish retained imputed-row analysis from ordinary fully observed
inference.

### Respect Bundle Versions

Version-1 kinase and signalome result bundles are rejected with a migration
message. Kinase version-2 bundles are also rejected as legacy because they did
not contain enough typed activity input semantics, profile identity, or
condition-summary aggregation metadata to reconstruct every valid
`KinaseActivityResult` faithfully. Regenerate old result bundles with the
current PhosPy version before loading them.

## Choose the Supported User Path

Use the `publish_*` helpers for user-facing output directories. Do not import
owner-internal bundle writers or loaders in application code.

If you need to preserve a result beyond the published tables, keep the input
data, request settings, package version, result provenance, and result caveats
together. The current public API does not promise reloadable workflow-result
objects as a general persistence interface.
