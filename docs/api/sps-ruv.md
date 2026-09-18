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

## What SPS means

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

## SPS reference quantitative-input contract

Every `SpsReferenceDataset` must already have the repository-authoritative
state used for condition-relative measurements:

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

## End-to-end API example

The example makes the preparation boundary visible. Here the study design has
already established `control` as the appropriate biological reference; that
scientific decision is not made by the code below.

```python
import pandas as pd

from phospy.advanced import (
    CorrectionMissingnessPolicy,
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
    SpsRuvBatchCorrectionConfig,
)
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
)


def prepare_reference(
    *,
    dataset_id: str,
    processed_log2: pd.DataFrame,
    condition_by_sample: dict[str, str],
) -> SpsReferenceDataset:
    """Assert a scientifically established control-centred log2 reference."""
    control_samples = [
        sample
        for sample, condition in condition_by_sample.items()
        if condition == "control"
    ]
    # This centring is valid here only because the study design established
    # this control condition as the biological reference beforehand.
    control_baseline = processed_log2.loc[:, control_samples].mean(axis=1)
    condition_relative_log2 = processed_log2.subtract(control_baseline, axis=0)

    return SpsReferenceDataset.from_condition_relative_log2(
        dataset_id=dataset_id,
        intensities=condition_relative_log2,
        condition_by_sample=condition_by_sample,
        log2_scale_established_by="reference-processing-pipeline-v3",
        baseline_centering_established_by=(
            f"{dataset_id}:study-design-control-mean-v1"
        ),
        organism="human",
        baseline_context="untreated control condition",
        reference_context="human cell-line treatment time course",
        source_name="independent phosphoproteomics reference",
        source_version="2026-09",
        source_uri=f"https://example.org/references/{dataset_id}",
    )


reference_a = prepare_reference(
    dataset_id="reference-a",
    processed_log2=reference_a_processed_log2,
    condition_by_sample=reference_a_conditions,
)
reference_b = prepare_reference(
    dataset_id="reference-b",
    processed_log2=reference_b_processed_log2,
    condition_by_sample=reference_b_conditions,
)

discovery = SpsDiscoveryWorkflow().run(
    SpsDiscoveryRequest(
        reference_datasets=(reference_a, reference_b),
        config=SpsDiscoveryConfig(
            top_n=100,
            minimum_reference_datasets=2,
            minimum_datasets_per_site=2,
            minimum_shared_sites=50,
        ),
    )
)

controls = discovery.control_site_set
print(discovery.selected_site_keys)
print(discovery.discovery_identity)
print(discovery.provenance.selection_boundaries.to_payload())

correction = SpsRuvBatchCorrectionConfig(
    method="ruv_iii_style",
    control_site_set=controls,
    batch_column="batch",
    condition_columns=("condition",),
    replicate_column="replicate_set",
    missingness_policy=CorrectionMissingnessPolicy(),
    n_unwanted_factors=1,
)

target_dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=target_phospho,
        site_metadata=target_site_metadata,
        sample_metadata=target_sample_metadata,
        organism="human",
        input_intensity_scale="log2",
        preprocessing_config=DatasetPreprocessingConfig(
            batch_correction=correction,
        ),
    )
)
```

The target matrix above is ordinary governed log2 target data; it need not be
centred like an SPS discovery reference. Control keys must map to target
`site_key` rows and satisfy the correction validator.

## Discovery behaviour and provenance

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
reference evidence from behaving like favorable stability evidence. A site
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
`phospy-sps-discovery-result-v1` payloads omitted that field; loading v1
explicitly validates the legacy identity, reconstructs each quantile from the
serialized stability-score tie groups, and returns a v2 result. A v1 payload
carrying v2-only quantile fields is rejected rather than interpreted
ambiguously.

Algorithm version `2.0.0` introduced the evidence-count-first ordering and
serializes
`contributing_dataset_count_descending_then_consensus_stability_descending` in
provenance. Version `1.0.0` results deserialize with their historical
score-first ordering and omit that newer field, so persisted results are not
silently reinterpreted. Complete-reference results have a constant contributor
count and retain their PhosR-compatible score and ordering.

Persist `discovery.to_payload()` as JSON-compatible data. Later,
`SpsDiscoveryResult.from_payload(...)` revalidates it and reconstructs a usable
`ControlSiteSet`; the source matrices are not needed for correction:

```python
from phospy.advanced import SpsDiscoveryResult, SpsRuvBatchCorrectionConfig

restored = SpsDiscoveryResult.from_payload(saved_discovery_payload)
reused_correction = SpsRuvBatchCorrectionConfig(
    method="sps_ruv_style",
    control_site_set=restored.control_site_set,
    batch_column="batch",
    condition_columns=("condition",),
    missingness_policy=correction_missingness_policy,
    n_unwanted_factors=1,
)
```

Transfer a discovered set only to a compatible biological and identifier
context. Persistence makes reuse technically possible; it does not establish
that a control set from another organism, tissue, perturbation, or acquisition
regime is scientifically appropriate.

## Choose the correction method

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

## When to use it

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
