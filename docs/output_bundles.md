# Output Bundles

PhosPy has two related output paths.

1. Publisher helpers write simple workflow output directories.
2. Bundle services write reloadable output bundles with a manifest and config
   snapshot.

Use the simple publisher layout when you only need files. Use bundle services when you
want to save and reload a workflow result object.

## Publisher Output Layout

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

Supported table formats are `csv`, `tsv`, and `parquet`.

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

## Reloadable Bundle Services

Public bundle services live in:

```python
from phospy.io.bundles.kinase import load_kinase_workflow_bundle
from phospy.io.bundles.kinase import save_kinase_workflow_bundle
from phospy.io.bundles.signalome import load_signalome_workflow_bundle
from phospy.io.bundles.signalome import save_signalome_workflow_bundle
```

Bundle writers require a config snapshot:

```python
from pathlib import Path

from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot
from phospy.io.bundles.kinase import save_kinase_workflow_bundle

snapshot = KinaseWorkflowConfigSnapshot.from_request(request)
written = save_kinase_workflow_bundle(
    result,
    Path("bundle-out"),
    config_snapshot=snapshot,
    output_format="csv",
)
```

A saved bundle contains:

- data tables
- reference tables
- workflow output tables
- `manifest.json`
- `config/snapshot.json`

Loaders return `LoadedKinaseWorkflowBundle` or `LoadedSignalomeWorkflowBundle`.
Each object contains the reconstructed result, the config snapshot, and the
manifest version.

## Bundle Safety Rules

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
`summary`, `input_semantics`, `profile_metadata`, and `tables`. Enabled activity
requires `input_semantics` and `profile_metadata` objects. Disabled activity
requires both fields to be `null`.

KSEA-style and ssGSEA-style runs also emit an activity `statistics_table` with
method statistics, p-values when available, optional q-values, substrate counts,
background counts, and computability statuses. The statistics-table row
identity is `profile_id`.
For `ActivityProfileAxis.SAMPLE`, `profile_id` contains sample IDs; for
`ActivityProfileAxis.CONDITION_SUMMARY`, it contains condition-summary profile
IDs and an optional `condition` column may also be present with the same values;
for `ActivityProfileAxis.CONTRAST`, it contains contrast IDs; and for
`ActivityProfileAxis.EFFECT`, it contains neutral effect profile IDs. Sample,
contrast, and effect statistics tables do not persist a `condition` column.
These outputs are exploratory substrate-supported scores, not direct causal
kinase activation evidence.

Bundle loaders reconstruct `dataset.processing_state` and
`dataset.intensity_scale_state` from saved payloads, including mixed
corrected/uncorrected total-protein quantitative meaning and row-level
total-protein correction diagnostics.

Structured result caveats are persistence data, not display-only warnings.
Current reloadable bundle services persist caveats for kinase and signalome
results, and `DifferentialAnalysisResult.to_payload()` serializes differential
caveats for handoff. Any future differential bundle writer must preserve the
same `caveats` payload, including exploratory single-biological-replicate
caveats that distinguish computable output from production-supported
inferential support.

Version-1 kinase and signalome result bundles are rejected with a migration
message. Kinase version-2 bundles are also rejected as legacy because they did
not contain enough typed activity input semantics, profile identity, or
condition-summary aggregation metadata to reconstruct every valid
`KinaseActivityResult` faithfully. Regenerate old result bundles with the
current PhosPy version before loading them.

## Choosing the Right Output Path

| Need | Use |
| --- | --- |
| Write workflow output tables from Python | `publish_*` helpers |
| Write result tables from Python | `publish_*` helpers |
| Save and reload a result object | bundle services |
| Keep config snapshot beside outputs | bundle services |
