# SPS Discovery and RUV-Style Correction

PhosPy keeps three scientific operations separate:

1. **Reference preparation** establishes the quantitative meaning, log2 scale,
   and biologically appropriate reference/control baseline of independent SPS
   evidence.
2. **SPS discovery** ranks consensus-stable phosphosites in those prepared
   references and returns a governed `ControlSiteSet`.
3. **RUV correction** uses that negative-control evidence to correct a target
   experiment.

They are not one automatic operation. In particular, correction does not need
the original SPS reference matrices once a `ControlSiteSet` has been produced.

```text
SPS reference datasets                target experiment
----------------------                -----------------
condition-relative log2               normal governed PhosPy input
        |                                      |
        v                                      |
  SPS discovery                                |
        |                                      |
        +------------ ControlSiteSet ----------+
                               |
                               v
                         correction
```

## What SPS Means

A stable phosphosite (SPS) is a phosphosite whose condition-relative signal is
consistently small across multiple suitable reference phosphoproteomics
datasets. Such sites can provide negative-control evidence for estimating
technical or otherwise unwanted variation.

Conventional housekeeping controls are difficult in phosphoproteomics because
phosphorylation is dynamic, condition dependent, and context specific. A site
that is stable in one tissue, perturbation, time course, or acquisition setting
may be regulated in another. Multiple independent references reduce reliance
on a single experiment: PhosPy scores within each reference and then combines
the evidence into a consensus ranking.

SPS discovery identifies stable candidates from already valid reference
evidence. It does **not**:

- choose or infer the biological reference condition;
- turn absolute abundance into condition-relative measurements;
- infer quantitative meaning from values, column names, or condition labels;
- decide that low variation in the target experiment is biological evidence of
  stability; or
- correct a target matrix.

Normally use independent reference evidence appropriate to the organism,
tissue, perturbation, and experimental context. Automatically declaring sites
stable because they change little in the same target experiment being assessed
can make the controls circular and can remove real biology. Reuse of target
data for discovery requires a deliberate, documented scientific justification.

## SPS Reference Quantitative-Input Contract

Every `SpsReferenceDataset` must already have the governed quantitative state
required for condition-relative measurements:

- one supported, typed organism shared by every reference in the request;
- explicit biological baseline and reference-context descriptions;
- a reconstructable source name, version, and URI;
- `IntensityScaleKind.LOG2` scale;
- `QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE` meaning;
- zero denoting the established reference/control baseline; and
- matrix-bound establishment evidence for both the log2 scale and the
  quantitative-meaning transition.

In practical terms, the reference matrix must contain condition-relative log2
values. The biological reference must have been established before discovery.
An arbitrary absolute-abundance matrix, an unknown quantitative state, linear
abundance, or unbound/missing establishment evidence is invalid and is rejected
before ranking. SPS discovery never examines the numbers and guesses what they
mean.

```text
processed log2 abundance
    |
    v
choose/establish biologically appropriate reference condition
    |
    v
reference-centre measurements
    |
    v
condition-relative log2 values (established baseline = zero)
    |
    v
SPS discovery
```

Subtracting a convenient condition is not automatically scientifically valid.
The caller remains responsible for choosing and documenting an appropriate
biological reference and for ensuring that the preparation operation has the
claimed meaning.

This strict contract applies to the **reference datasets used for SPS
discovery**. It is not a blanket requirement that every target matrix passed to
`sps_ruv_style` or `ruv_iii_style` be an SPS-style condition-relative matrix.
The target experiment follows the normal governed PhosPy dataset input and
preprocessing contracts.

## Prepare the Reference Evidence

The study design must establish the biological baseline before SPS discovery.
The teaching references below already contain control-relative log2 values, so
zero represents the established control baseline. The code records that fact;
it does not infer or create the baseline.

Start with the shared protein-scoped `site_key` identities.

```python
import pandas as pd

site_keys = (
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=MAPK14|residue=Y|position=182",
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=AKT1|residue=T|position=308",
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=SRC|residue=Y|position=416",
)
```

Each reference may have its own samples and conditions. These first values make
MAPK14 and SRC relatively stable while AKT1 changes strongly.

```python
reference_a_values = pd.DataFrame(
    {
        "control_1": [0.0, 0.0, 0.0],
        "control_2": [0.0, 0.0, 0.0],
        "treated_1": [0.1, 2.0, 0.2],
        "treated_2": [0.1, 2.0, 0.2],
    },
    index=pd.Index(site_keys, name="site_key"),
)
reference_a_conditions = {
    "control_1": "control",
    "control_2": "control",
    "treated_1": "treated",
    "treated_2": "treated",
}
```

```python
reference_b_values = pd.DataFrame(
    {
        "control_1": [0.0, 0.0, 0.0],
        "control_2": [0.0, 0.0, 0.0],
        "treated_1": [0.2, 3.0, 0.3],
        "treated_2": [0.2, 3.0, 0.3],
    },
    index=pd.Index(site_keys, name="site_key"),
)
reference_b_conditions = {
    "control_1": "control",
    "control_2": "control",
    "treated_1": "treated",
    "treated_2": "treated",
}
```

`from_condition_relative_log2()` is an assertion boundary. It copies the
matrix and binds the supplied scale and baseline evidence to that exact matrix;
it does not transform the values.

```python
from phospy.advanced import SpsReferenceDataset

reference_a = SpsReferenceDataset.from_condition_relative_log2(
    dataset_id="reference_a",
    intensities=reference_a_values,
    condition_by_sample=reference_a_conditions,
    log2_scale_established_by="reference_a log2 preprocessing",
    baseline_centering_established_by="reference_a control subtraction",
    organism="rat",
    baseline_context="control condition",
    reference_context="synthetic rat correction tutorial",
    source_name="governed-reference-a",
    source_version="2026-09",
    source_uri="https://example.org/references/reference-a",
)
```

```python
reference_b = SpsReferenceDataset.from_condition_relative_log2(
    dataset_id="reference_b",
    intensities=reference_b_values,
    condition_by_sample=reference_b_conditions,
    log2_scale_established_by="reference_b log2 preprocessing",
    baseline_centering_established_by="reference_b control subtraction",
    organism="rat",
    baseline_context="control condition",
    reference_context="synthetic rat correction tutorial",
    source_name="governed-reference-b",
    source_version="2026-09",
    source_uri="https://example.org/references/reference-b",
)
```

```python
reference_datasets = (reference_a, reference_b)
```

## Configure SPS Discovery

The configuration selects at most two controls and requires both references to
contribute to each ranked site. The two-reference floors are the smallest
supported consensus path.

```python
from phospy.advanced import SpsDiscoveryConfig

discovery_config = SpsDiscoveryConfig(
    top_n=2,
    minimum_reference_datasets=2,
    minimum_datasets_per_site=2,
    minimum_shared_sites=2,
)
```

## Construct and Run the Discovery Request

```python
from phospy.advanced import SpsDiscoveryRequest

discovery_request = SpsDiscoveryRequest(
    reference_datasets=reference_datasets,
    config=discovery_config,
)
```

```python
from phospy.advanced import SpsDiscoveryWorkflow

discovery_workflow = SpsDiscoveryWorkflow()
discovery_result = discovery_workflow.run(discovery_request)
```

Inspect the selected control identities and the attrition boundary.

```python
selected_site_keys = discovery_result.selected_site_keys
selection_boundaries = discovery_result.provenance.selection_boundaries.to_payload()

print("Selected controls:", selected_site_keys)
print("Selection boundaries:", selection_boundaries)
```

The result selects the MAPK14 and SRC keys. The selection-boundary payload
reports three unique, valid, and ranked sites and two selected sites. Stability
here is a property of the supplied reference evidence, not proof that these
sites are universally stable.

## Prepare the Target Experiment

The correction target follows the ordinary dataset contract. It need not be
control-centred like an SPS discovery reference. This complete log2 matrix has
three sites measured in four samples.

```python
target_phospho = pd.DataFrame(
    {
        "sample_1": [10.0, 5.0, 20.0],
        "sample_2": [10.0, 9.0, 20.0],
        "sample_3": [14.0, 8.0, 28.0],
        "sample_4": [14.0, 12.0, 28.0],
    },
    index=pd.Index(
        ["MAPK14;Y182;", "AKT1;T308;", "SRC;Y416;"],
        name="site_id",
    ),
)
```

```python
target_site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "AKT1", "SRC"],
        "protein_id": ["MAPK14", "AKT1", "SRC"],
        "site": ["Y182", "T308", "Y416"],
        "site_sequence": [
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "T" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
        ],
        "localisation_confidence": [0.95, 0.92, 0.98],
    },
    index=target_phospho.index,
)
```

Batch and protected condition labels come from aligned sample metadata. The
replicate column is recorded for diagnostics by `sps_ruv_style`; it does not
activate replicate-aware RUV-III estimation.

```python
target_sample_metadata = pd.DataFrame(
    {
        "batch": ("run_1", "run_1", "run_2", "run_2"),
        "condition": ("control", "treated", "control", "treated"),
        "replicate": ("r1", "r2", "r2", "r1"),
    },
    index=target_phospho.columns,
)
```

## Configure the Correction

The target is complete, so the default missingness policy performs no temporary
completion. The correction requests one unwanted factor and consumes the
discovered `ControlSiteSet` directly.

```python
from phospy.advanced import CorrectionMissingnessPolicy

correction_missingness = CorrectionMissingnessPolicy()
```

```python
from phospy.advanced import SpsRuvBatchCorrectionConfig

correction_config = SpsRuvBatchCorrectionConfig(
    method="sps_ruv_style",
    control_site_set=discovery_result.control_site_set,
    batch_column="batch",
    condition_columns=("condition",),
    replicate_column="replicate",
    missingness_policy=correction_missingness,
    n_unwanted_factors=1,
    diagnostics_enabled=True,
    provenance_enabled=True,
)
```

Place that operation in the dataset preprocessing configuration.

```python
from phospy.api import DatasetPreprocessingConfig

preprocessing_config = DatasetPreprocessingConfig(
    batch_correction=correction_config,
)
```

## Build and Inspect the Corrected Dataset

```python
from phospy.api import DatasetBuildRequest, Organism

target_request = DatasetBuildRequest(
    phospho=target_phospho,
    site_metadata=target_site_metadata,
    sample_metadata=target_sample_metadata,
    organism=Organism.RAT,
    input_intensity_scale="log2",
    preprocessing_config=preprocessing_config,
)
```

```python
from phospy import AnalysisReadyDatasetBuilder

dataset_builder = AnalysisReadyDatasetBuilder()
corrected_dataset = dataset_builder.run(target_request)
```

Inspect the public correction report and confirm that the output matrix changed.

```python
correction_report = corrected_dataset.preprocessing_report.batch_correction
corrected_shape = corrected_dataset.phospho.shape
matrix_changed = not corrected_dataset.phospho.equals(target_phospho)

print("Corrected shape:", corrected_shape)
print("Correction method:", correction_report.method)
print("Correction status:", correction_report.status)
print("Matrix changed:", matrix_changed)
```

The expected shape is `(3, 4)`, the method is `sps_ruv_style`, the status is
`applied`, and `matrix_changed` is `True`. Correction estimates and removes a
supported unwanted-variation component. It does not prove that all technical
variation was removed or that biological effects were preserved outside the
protected design and control assumptions.

## Discovery Behaviour and Provenance

Each reference supplies its own sample-to-condition mapping. References need
not share samples or condition labels, but they do share governed `site_key`
identity. Discovery averages available biological replicates within each
condition, measures each site's largest absolute condition mean in each
reference, ranks smaller magnitudes as more stable, and combines the ranks with
a Fisher-style consensus score.

The defaults require at least two reference datasets and contribution from at
least two references per ranked site. Configure
`minimum_reference_datasets`, `minimum_datasets_per_site`, and
`minimum_shared_sites` explicitly when the study needs different valid floors.
A site may contribute partially: a completely missing site-condition block
prevents that site from contributing to that reference, but it may still rank
when it meets `minimum_datasets_per_site`. Individual missing replicate cells
are ignored only if a finite value remains for that site and condition.
Infinite values are rejected.

Cross-count ordering is an explicit evidence policy. For site \(i\), let
\(k_i\) be its contributing-reference count, \(n_r\) the candidate count in
reference \(r\), and \(\bar a_{ir}\) its tie-aware stability midrank. Discovery
computes

\[
q_{ir}=(n_r-\bar a_{ir}+1/2)/n_r,\quad
T_i=-2\sum_r\log q_{ir},\quad
C_i=\Pr\{\chi^2_{2(k_i-1)}\ge T_i\}.
\]

The final ordering is lexicographic in
`(-k_i, -C_i, site_key_i)`: more contributing references first, then larger
consensus stability score within an equal-count stratum, then ascending
`site_key`. This conservative corroboration policy keeps missing
reference evidence from behaving like favourable stability evidence. A site
with stronger per-reference stability can therefore rank below a site with
more independent reference contributions.

`top_n` is an upper bound. If fewer eligible sites remain, all eligible sites
are selected and the requested and actual counts are recorded. Equal
within-reference values receive equal ranks; exact final consensus ties are
ordered by ascending `site_key`, so selection is deterministic. Each
`SpsDatasetSiteStatistic.rank` is the reporting competition rank. Its
`rank_quantile` is derived from the tie-aware estimator midrank and is the
value supplied directly to the Fisher-style consensus, so serialized results
retain the complete per-reference consensus inputs.

The discovery record makes the following auditable without embedding complete
reference matrices:

- typed organism, biological baseline/reference context, reconstructable
  reference source identities, sample-condition mappings, and matrix
  fingerprints;
- scale, `CONTRAST_LOG2_FOLD_CHANGE` meaning, and baseline-establishment
  evidence;
- SPS algorithm ID/version and parameters;
- the explicit cross-count ordering policy;
- per-reference contribution and selection/attrition counts;
- per-reference stability scores, reporting ranks, and consensus-input rank
  quantiles; and
- ranked and selected controls.

Each result also exposes a deterministic `discovery_identity`. Its
`sha256-stable-json-v1` digest is computed from the deterministic matrix-free
discovery result: reference identities and quantitative fingerprints, governed
quantitative evidence, algorithm/configuration, selection boundaries, ranking,
and selected controls. It is not a hash of the selected `site_key` values alone;
two runs that select the same controls from different evidence remain distinct.
The validated organism is also copied to the generated `ControlSiteSet`, so a
target dataset with a conflicting organism is rejected by the existing control
compatibility boundary.

The identity is serialized with the `phospy-sps-discovery-result-v2` schema,
which requires every per-reference `rank_quantile`. Baseline
`phospy-sps-discovery-result-v1` payloads omit that field; loading v1 explicitly
validates the v1 identity, reconstructs each quantile from the
serialized stability-score tie groups, and returns a v2 result. A v1 payload
carrying v2-only quantile fields is rejected rather than interpreted
ambiguously.

Algorithm version `2.0.0` uses evidence-count-first ordering and serializes
`contributing_dataset_count_descending_then_consensus_stability_descending` in
provenance. Version `1.0.0` results deserialize with score-first ordering and
omit that field, so persisted results are not
silently reinterpreted. Complete-reference results have a constant contributor
count and retain their PhosR-compatible score and ordering.

After running the primary example, persist `discovery_result.to_payload()` as
JSON-compatible data. Later,
`SpsDiscoveryResult.from_payload(...)` revalidates it and reconstructs a usable
`ControlSiteSet`; the source matrices are not needed for correction:

```python
from phospy.advanced import SpsDiscoveryResult, SpsRuvBatchCorrectionConfig

saved_discovery_payload = discovery_result.to_payload()
restored = SpsDiscoveryResult.from_payload(saved_discovery_payload)
reused_correction = SpsRuvBatchCorrectionConfig(
    method="sps_ruv_style",
    control_site_set=restored.control_site_set,
    batch_column="batch",
    condition_columns=("condition",),
    missingness_policy=correction_missingness,
    n_unwanted_factors=1,
)
```

Transfer a discovered set only to a compatible biological and identifier
context. Persistence makes reuse technically possible; it does not establish
that a control set from another organism, tissue, perturbation, or acquisition
regime is scientifically appropriate.

## Choose the Correction Method

| Method | Scientific role | Replicate semantics |
|---|---|---|
| `sps_ruv_style` | Existing PhosPy-native SPS/control-based unwanted-variation correction after protected-design handling | `replicate_column` is optional and provenance/diagnostics-only |
| `ruv_iii_style` | Replicate-aware native method related to RUV-III and informed by PhosR's use of negative controls | `replicate_column` is required; replicate-set membership participates directly in estimation |

`sps_ruv_style` is not RUV-III. The two estimators are numerically distinct and
can legitimately return different corrected values from the same controls and
target. For `ruv_iii_style`, each sample has exactly one replicate-set
assignment, sets provide the within-set residual structure used to estimate
unwanted factors, and a set must not cross protected condition strata. Under
the supported PhosPy contract, every replicate set must contain at least two
samples; singleton replicate groups are rejected because they are not suitable
for estimating the intended RUV-III structure.

Choose `k` (`n_unwanted_factors`) only when control count and
replicate-residual rank support it. PhosPy rejects a requested `k` that exceeds
the estimable latent rank. It does not silently cap or reduce `k`, so the
requested and effective factor counts remain the same for a valid fit. The
pinned `ruv` reference is more permissive about singleton replicate sets and
non-estimable `k`; these are deliberate PhosPy input-validation restrictions,
not numerical differences in the estimator over inputs supported by both
implementations.

### Estimator Rejection Boundaries and Recovery

`ruv_iii_style` fails closed when its requested factor space cannot be
identified unambiguously. Its user-visible restrictions are:

- **Singleton replicate sets:** every replicate set must contain at least two
  samples. Merge or redesign singleton sets only when the resulting assignment
  is scientifically valid; otherwise collect another replicate or do not use
  `ruv_iii_style`.
- **Rank restrictions:** the replicate mapping, replicate-residual matrix, and
  negative-control loading matrix must have the rank required by the requested
  fit. Review replicate assignments and controls, remove redundant or
  non-informative controls, or choose a smaller `k`; do not add noise merely to
  force full rank.
- **Control-count restrictions:** at least `k` eligible negative-control sites
  must reach the estimator. Supply more defensible controls, resolve control
  attrition, or reduce `k`.
- **Non-estimable `k`:** `k` cannot exceed either the available within-set
  residual degrees of freedom or the numerical rank of the
  replicate-residual matrix. Increase supported replication or reduce `k`.
  PhosPy never silently substitutes a smaller effective value.
- **Tied singular-value cutoff:** `k` must not split a numerically tied block
  of replicate-residual singular values. A split would select an arbitrary
  orientation-dependent subset of one latent subspace. Choose `k` immediately
  below or above the complete tied block, subject to the other rank and control
  restrictions, and rerun. Reordering samples is not a recovery strategy.

The exception reports the failed boundary and relevant count, rank, or `k`.
Treat these messages as study-design and estimator diagnostics. Do not perturb
measurements, relabel replicates, or replace controls solely to bypass a
rejection.

Correction provenance records the selected method, selected control identities,
`ControlSiteSet` source metadata, replicate definition when relevant, `k`,
missing-data strategy, control eligibility/attrition, estimator diagnostics,
warnings, and matrix/mask fingerprints. For SPS-derived controls, that source
metadata carries the same immutable `sps_discovery_identity` as the discovery
result. The full discovery provenance remains with `SpsDiscoveryResult`; the
correction record stores its digest link and does not copy the complete SPS
provenance or reference matrices. Manually or externally defined control sets
may legitimately have no SPS discovery identity.

## Missingness

`sps_ruv_style` requires a complete numeric correction-stage matrix and does
not use `row_median_temporary` as permission to accept actual missing cells.

`ruv_iii_style` can temporarily fill actual missing cells with the governed
row-median method when the explicit `CorrectionMissingnessPolicy` and
`ObservationMask` permit it. That completion exists only so the numerical
estimator can run. Actual missing positions are restored to missing in the
executor output, and completed values are never marked as observed evidence.
Cells already imputed by upstream preprocessing can remain numeric while their
observation-mask state remains false. Analysis-ready construction still
requires its existing complete-output contract, so upstream missing-data
preprocessing remains authoritative.

Neither SPS reference preparation nor RUV-III's internal completion changes
PhosPy's existing missing-data and imputation contracts.

## When to Use It

This workflow is appropriate when you have multiple relevant, independently
prepared phosphoproteomics references, can defend their biological baselines,
can map the selected sites to the target, and have target metadata and controls
adequate for the chosen estimator and `k`.

Do not use it when reference meaning or baseline is unknown, the inputs are
arbitrary absolute abundances, references are biologically incompatible with
the target, control coverage is inadequate, replicate sets are missing or
invalid for RUV-III, or the experimental design confounds unwanted and
protected biological structure. It is also inappropriate to discover controls
automatically from the target merely because they appear unchanged unless that
reuse has been scientifically justified in advance.

PhosPy implements SPS discovery, `sps_ruv_style`, and `ruv_iii_style`. Their
design is informed by SPS/RUV literature and PhosR, but the package does not
claim full PhosR workflow equivalence. See [Comparison With
PhosR](../parity.md#sps-and-ruv-iii-external-reference-evidence) for the narrow
pinned external evidence and its exclusions.

## Complete Runnable Example

```python
import pandas as pd

site_keys = (
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=MAPK14|residue=Y|position=182",
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=AKT1|residue=T|position=308",
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=SRC|residue=Y|position=416",
)

reference_a_values = pd.DataFrame(
    {
        "control_1": [0.0, 0.0, 0.0],
        "control_2": [0.0, 0.0, 0.0],
        "treated_1": [0.1, 2.0, 0.2],
        "treated_2": [0.1, 2.0, 0.2],
    },
    index=pd.Index(site_keys, name="site_key"),
)
reference_a_conditions = {
    "control_1": "control",
    "control_2": "control",
    "treated_1": "treated",
    "treated_2": "treated",
}

reference_b_values = pd.DataFrame(
    {
        "control_1": [0.0, 0.0, 0.0],
        "control_2": [0.0, 0.0, 0.0],
        "treated_1": [0.2, 3.0, 0.3],
        "treated_2": [0.2, 3.0, 0.3],
    },
    index=pd.Index(site_keys, name="site_key"),
)
reference_b_conditions = {
    "control_1": "control",
    "control_2": "control",
    "treated_1": "treated",
    "treated_2": "treated",
}

from phospy.advanced import SpsReferenceDataset

reference_a = SpsReferenceDataset.from_condition_relative_log2(
    dataset_id="reference_a",
    intensities=reference_a_values,
    condition_by_sample=reference_a_conditions,
    log2_scale_established_by="reference_a log2 preprocessing",
    baseline_centering_established_by="reference_a control subtraction",
    organism="rat",
    baseline_context="control condition",
    reference_context="synthetic rat correction tutorial",
    source_name="governed-reference-a",
    source_version="2026-09",
    source_uri="https://example.org/references/reference-a",
)

reference_b = SpsReferenceDataset.from_condition_relative_log2(
    dataset_id="reference_b",
    intensities=reference_b_values,
    condition_by_sample=reference_b_conditions,
    log2_scale_established_by="reference_b log2 preprocessing",
    baseline_centering_established_by="reference_b control subtraction",
    organism="rat",
    baseline_context="control condition",
    reference_context="synthetic rat correction tutorial",
    source_name="governed-reference-b",
    source_version="2026-09",
    source_uri="https://example.org/references/reference-b",
)

reference_datasets = (reference_a, reference_b)

from phospy.advanced import SpsDiscoveryConfig

discovery_config = SpsDiscoveryConfig(
    top_n=2,
    minimum_reference_datasets=2,
    minimum_datasets_per_site=2,
    minimum_shared_sites=2,
)

from phospy.advanced import SpsDiscoveryRequest

discovery_request = SpsDiscoveryRequest(
    reference_datasets=reference_datasets,
    config=discovery_config,
)

from phospy.advanced import SpsDiscoveryWorkflow

discovery_workflow = SpsDiscoveryWorkflow()
discovery_result = discovery_workflow.run(discovery_request)

selected_site_keys = discovery_result.selected_site_keys
selection_boundaries = discovery_result.provenance.selection_boundaries.to_payload()

print("Selected controls:", selected_site_keys)
print("Selection boundaries:", selection_boundaries)

target_phospho = pd.DataFrame(
    {
        "sample_1": [10.0, 5.0, 20.0],
        "sample_2": [10.0, 9.0, 20.0],
        "sample_3": [14.0, 8.0, 28.0],
        "sample_4": [14.0, 12.0, 28.0],
    },
    index=pd.Index(
        ["MAPK14;Y182;", "AKT1;T308;", "SRC;Y416;"],
        name="site_id",
    ),
)

target_site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "AKT1", "SRC"],
        "protein_id": ["MAPK14", "AKT1", "SRC"],
        "site": ["Y182", "T308", "Y416"],
        "site_sequence": [
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "T" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
        ],
        "localisation_confidence": [0.95, 0.92, 0.98],
    },
    index=target_phospho.index,
)

target_sample_metadata = pd.DataFrame(
    {
        "batch": ("run_1", "run_1", "run_2", "run_2"),
        "condition": ("control", "treated", "control", "treated"),
        "replicate": ("r1", "r2", "r2", "r1"),
    },
    index=target_phospho.columns,
)

from phospy.advanced import CorrectionMissingnessPolicy

correction_missingness = CorrectionMissingnessPolicy()

from phospy.advanced import SpsRuvBatchCorrectionConfig

correction_config = SpsRuvBatchCorrectionConfig(
    method="sps_ruv_style",
    control_site_set=discovery_result.control_site_set,
    batch_column="batch",
    condition_columns=("condition",),
    replicate_column="replicate",
    missingness_policy=correction_missingness,
    n_unwanted_factors=1,
    diagnostics_enabled=True,
    provenance_enabled=True,
)

from phospy.api import DatasetPreprocessingConfig

preprocessing_config = DatasetPreprocessingConfig(
    batch_correction=correction_config,
)

from phospy.api import DatasetBuildRequest, Organism

target_request = DatasetBuildRequest(
    phospho=target_phospho,
    site_metadata=target_site_metadata,
    sample_metadata=target_sample_metadata,
    organism=Organism.RAT,
    input_intensity_scale="log2",
    preprocessing_config=preprocessing_config,
)

from phospy import AnalysisReadyDatasetBuilder

dataset_builder = AnalysisReadyDatasetBuilder()
corrected_dataset = dataset_builder.run(target_request)

correction_report = corrected_dataset.preprocessing_report.batch_correction
corrected_shape = corrected_dataset.phospho.shape
matrix_changed = not corrected_dataset.phospho.equals(target_phospho)

print("Corrected shape:", corrected_shape)
print("Correction method:", correction_report.method)
print("Correction status:", correction_report.status)
print("Matrix changed:", matrix_changed)
```

## Continue

- [Prepare a dataset](dataset-build-workflow.md)
- [Differential analysis](differential-analysis.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
- [Comparison with PhosR](../parity.md#sps-and-ruv-iii-external-reference-evidence)
