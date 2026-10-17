# ADR-0054: Public API Reset and Compatibility Removal

## Status

- **ADR ID:** ADR-0054
- **Title:** Public API Reset and Compatibility Removal
- **Status:** Accepted
- **Date:** 26 September 2026
- **Accepted:** 26 September 2026
- **Decision type:** API and compatibility governance
- **Evidence base:** Source snapshot `138b7764bcb75bc3506dc58cdd2d8b044cf70d4c` (PhosPy 1.7.6)
- **Related decision:** ADR-0053, Domain-First Capability Architecture and Dependency Governance

## Context

ADR-0053 establishes domain-first implementation ownership. Public import
routes are separate: an object may be owned by a domain package while being
re-exported through a user-facing facade.

The current surface contains five root exports, 68 stable aggregate exports,
124 advanced aggregate exports, additional facade submodules, workflow and
table routes, four groups of semi-public science routes, and 923 retained
deprecation records. Of those records, 901 are import-route combinations and
22 are aliases or legacy runtime behaviours.

PhosPy is a beta package and has explicitly chosen a clean break rather than
carrying historical physical paths into the domain-first architecture. The
reset must still design the post-cutover API deliberately rather than preserve
or remove names merely because they exist today.

## Decision drivers

1. `phospy.api` remains the official ordinary import route.
2. The package root remains deliberately minimal.
3. Specialist configuration remains available without enlarging the ordinary
   surface.
4. Domain ownership and external import convenience remain separate.
5. Python importability does not create a support promise.
6. Compatibility code must not constrain physical migration.
7. Supported scientific behaviour remains unchanged during structural work.
8. Compatibility removal and physical relocation are separate reviewable
   changes even though neither provides a backwards-compatibility window.

## Decision

### Supported external surfaces

After the reset, PhosPy has three supported aggregate surfaces.

| Surface | Purpose | Support rule |
| --- | --- | --- |
| `phospy` | Minimal convenience imports | Only names in `phospy.__all__` are supported. |
| `phospy.api` | Official ordinary user contract | Only names in `phospy.api.__all__` are supported. |
| `phospy.advanced` | Documented specialist contract | Only names in `phospy.advanced.__all__` are supported. |

`phospy.api` and `phospy.advanced` are permanent composition facades, not
compatibility shims. They own no scientific behaviour and re-export the exact
objects defined by logical owners. No facade submodule or implementation
package is an additional supported external route unless a later ADR adds it.

### Root convenience surface

The target root is exactly:

```python
(
    "AnalysisReadyDatasetBuilder",
    "KinaseWorkflow",
    "SignalomeWorkflow",
)
```

This preserves the established product decision that the root offers the
dataset entry point and the two primary downstream workflows. Differential and
Enrichment are first-class architectural capabilities and remain fully
supported through `phospy.api`; architectural status does not require a root
alias. `AnalysisReadyPhosphoDataset` remains stable in `phospy.api` but is not a
root executable entry point.

### Ordinary API classification

The current 68-name stable inventory was reviewed name by name. Sixty-six names
remain in `phospy.api`. `DesignMatrix` and `ContrastMatrix` become internal
Differential computation models: the public `DifferentialAnalysisRequest`
accepts `ExperimentalDesign` and `Contrast` values, and ordinary documentation
does not require callers to construct the lower-level matrices.

| Name | Logical owner | Decision |
| --- | --- | --- |
| `AnalysisReadyDatasetBuilder` | Datasets | Retain ordinary |
| `AnalysisReadyPhosphoDataset` | Datasets | Retain ordinary |
| `ReferenceBundleBuilder` | References | Retain ordinary |
| `DatasetBuildRequest` | Datasets | Retain ordinary |
| `PhosphositeImporter` | IO/importers | Retain ordinary |
| `PhosphositeImportRequest` | IO/importers | Retain ordinary |
| `PhosphositeImportResult` | IO/importers | Retain ordinary |
| `DatasetLocalisationConfig` | Preprocessing | Retain ordinary |
| `DatasetPreprocessingConfig` | Preprocessing | Retain ordinary |
| `BatchCovariate` | Differential/design | Retain ordinary |
| `CategoricalCovariate` | Differential/design | Retain ordinary |
| `ContinuousCovariate` | Differential/design | Retain ordinary |
| `FixedEffectCovariate` | Differential/design | Retain ordinary |
| `ExperimentalDesign` | Differential/design | Retain ordinary |
| `SampleDesignRecord` | Differential/design | Retain ordinary |
| `Contrast` | Differential/design | Retain ordinary |
| `DesignMatrix` | Differential | Remove from public API; internal computation model |
| `ContrastMatrix` | Differential | Remove from public API; internal computation model |
| `DifferentialAnalysisRequest` | Differential | Retain ordinary |
| `EnrichmentDerivedQuantitativeSetProvenance` | Enrichment | Retain ordinary |
| `EnrichmentDerivedSetMissingValueRule` | Enrichment | Retain ordinary |
| `EnrichmentDerivedSetSourceResultKind` | Enrichment | Retain ordinary |
| `EnrichmentDerivedSetThresholdDirection` | Enrichment | Retain ordinary |
| `EnrichmentDerivedSetValueMeaning` | Enrichment | Retain ordinary |
| `EnrichmentDerivedSetValueScale` | Enrichment | Retain ordinary |
| `EnrichmentIdentifierSetProvenance` | Enrichment | Retain ordinary |
| `EnrichmentIdentifierSetSourceType` | Enrichment | Retain ordinary |
| `EnrichmentSet` | Enrichment | Retain ordinary |
| `EnrichmentSetCollection` | Enrichment | Retain ordinary |
| `EnrichmentWorkflowRequest` | Enrichment | Retain ordinary |
| `InputIntensityScaleEvidence` | Quantitative | Retain ordinary |
| `EnrichmentConfig` | Enrichment | Retain ordinary |
| `GeneSetCollection` | Enrichment | Retain ordinary |
| `PtmSetCollection` | Enrichment | Retain ordinary |
| `KinaseWorkflowRequest` | Kinase | Retain ordinary |
| `SignalomeWorkflowRequest` | Signalome | Retain ordinary |
| `all_pairwise_contrasts` | Differential/design | Retain ordinary |
| `contrasts_vs_control` | Differential/design | Retain ordinary |
| `DifferentialAnalysisWorkflow` | Differential | Retain ordinary |
| `EnrichmentWorkflow` | Enrichment | Retain ordinary |
| `KinaseWorkflow` | Kinase | Retain ordinary |
| `SignalomeWorkflow` | Signalome | Retain ordinary |
| `DifferentialAnalysisResult` | Differential | Retain ordinary |
| `EnrichmentResultRecord` | Enrichment | Retain ordinary |
| `EnrichmentWorkflowResult` | Enrichment | Retain ordinary |
| `KinaseActivityResult` | Kinase/activity | Retain ordinary |
| `KinasePredictionResult` | Kinase/prediction | Retain ordinary |
| `KinaseScoringResult` | Kinase/scoring | Retain ordinary |
| `KinaseWorkflowResult` | Kinase | Retain ordinary |
| `ResultCaveat` | Shared result diagnostics | Retain ordinary |
| `SignalomeWorkflowResult` | Signalome | Retain ordinary |
| `Organism` | Taxonomy | Retain ordinary |
| `ReferenceBundle` | References | Retain ordinary |
| `ReferenceBundleBuildRequest` | References | Retain ordinary |
| `ReferencePreset` | References | Retain ordinary |
| `IntensityScaleKind` | Quantitative | Retain ordinary |
| `PhosPyError` | Errors | Retain ordinary |
| `PhosPyInputError` | Errors | Retain ordinary |
| `UnsupportedInputFormatError` | Errors | Retain ordinary |
| `PhosPyValidationError` | Errors | Retain ordinary |
| `ContractValidationError` | Errors | Retain ordinary |
| `PhosPyReferenceError` | Errors | Retain ordinary |
| `ReferenceResolutionError` | Errors | Retain ordinary |
| `ReferenceCompatibilityError` | Errors | Retain ordinary |
| `PhosPyWorkflowError` | Errors | Retain ordinary |
| `WorkflowValidationError` | Errors | Retain ordinary |
| `WorkflowBoundaryError` | Errors | Retain ordinary |
| `SignalomeScaleError` | Signalome/errors | Retain ordinary |

### Advanced API classification

All 124 current advanced names remain supported after individual review. They
are required to construct documented non-default configurations, run SPS
discovery, supply governed control sites or Kinase Library resources, inspect
typed diagnostics, filter or rank Differential tables, or publish documented
outputs. Retention is based on specialist user need rather than historical
route compatibility.

| Name | Logical owner | Decision |
| --- | --- | --- |
| `ControlSiteAnnotation` | Preprocessing/control sites | Retain advanced |
| `ControlSiteSet` | Preprocessing/control sites | Retain advanced |
| `ControlSiteSourceMetadata` | Preprocessing/control sites | Retain advanced |
| `ControlSiteStatus` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryRequest` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryWorkflow` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryValidationError` | Preprocessing/control sites | Retain advanced |
| `SpsReferenceDataset` | Preprocessing/control sites | Retain advanced |
| `CorrectionMaskPolicy` | Preprocessing/batch correction | Retain advanced |
| `CorrectionMissingnessPolicy` | Preprocessing/batch correction | Retain advanced |
| `ObservationMask` | Preprocessing/batch correction | Retain advanced |
| `OriginallyMissingCellTracking` | Preprocessing/batch correction | Retain advanced |
| `TemporaryImputationMethod` | Preprocessing/batch correction | Retain advanced |
| `TemporaryImputationPolicy` | Preprocessing/batch correction | Retain advanced |
| `DatasetBatchCorrectionConfig` | Preprocessing | Retain advanced |
| `DatasetBatchCorrectionMethod` | Preprocessing | Retain advanced |
| `DatasetPreprocessingBatchCorrectionConfig` | Preprocessing | Retain advanced |
| `SpsRuvBatchCorrectionConfig` | Preprocessing/batch correction | Retain advanced for native PhosPy SPS/RUV-style correction |
| `SpsRuvBatchCorrectionMethod` | Preprocessing/batch correction | Retain advanced |
| `SpsDiscoveryConfig` | Preprocessing/control sites | Retain advanced |
| `DatasetComparisonBuildingConfig` | Preprocessing | Retain advanced |
| `DatasetComparisonBuildingPolicy` | Preprocessing | Retain advanced |
| `DatasetComparisonPair` | Preprocessing | Retain advanced |
| `DatasetGroupCoverageFilterConfig` | Preprocessing | Retain advanced |
| `DatasetIntensityTransformConfig` | Preprocessing/quantitative | Retain advanced |
| `DatasetIntensityTransformPolicy` | Preprocessing/quantitative | Retain advanced |
| `DatasetLocalisationMode` | Preprocessing/evidence | Retain advanced |
| `DatasetMissingDataConfig` | Preprocessing | Retain advanced |
| `DatasetMissingDataInputScale` | Preprocessing/quantitative | Retain advanced |
| `DatasetMissingDataKnnNoOverlapPolicy` | Preprocessing | Retain advanced |
| `DatasetMissingDataPolicy` | Preprocessing | Retain advanced |
| `DatasetNormalisationConfig` | Preprocessing | Retain advanced |
| `DatasetNormalisationPolicy` | Preprocessing | Retain advanced |
| `DatasetProteinAwarePreparationConfig` | Preprocessing | Retain advanced |
| `DatasetProteinAwarePreparationMappingPolicy` | Preprocessing | Retain advanced |
| `DatasetProteinAwarePreparationPolicy` | Preprocessing | Retain advanced |
| `DatasetRuvReadinessConfig` | Preprocessing/batch correction | Retain advanced |
| `DatasetSiteMatrixConfig` | Preprocessing | Retain advanced |
| `DatasetSiteMatrixDuplicateSitePolicy` | Preprocessing/phosphosites | Retain advanced |
| `DatasetSiteMatrixMissingDataPolicy` | Preprocessing | Retain advanced |
| `DatasetSiteMatrixPolicy` | Preprocessing | Retain advanced |
| `DatasetSiteSequenceConflictPolicy` | Preprocessing/phosphosites | Retain advanced |
| `DatasetSiteSequenceResolutionConfig` | Preprocessing | Retain advanced |
| `DatasetSiteSequenceResolutionMode` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionConfig` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionDuplicatePolicy` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionIdentityConfig` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionIdentityMatchingPolicy` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionIdentityMode` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionPolicy` | Preprocessing | Retain advanced |
| `DatasetTotalProteinCorrectionUnmatchedPolicy` | Preprocessing | Retain advanced |
| `DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1` | Differential | Retain advanced |
| `DifferentialAnalysisConfig` | Differential | Retain advanced |
| `DifferentialImputedValuePolicy` | Differential | Retain advanced |
| `DifferentialProteinAwareModelConfig` | Differential | Retain advanced |
| `DifferentialReliabilityProfile` | Differential | Retain advanced |
| `EmpiricalBayesConfig` | Differential | Retain advanced |
| `MultipleTestingConfig` | Differential/statistics | Retain advanced |
| `MultipleTestingCorrection` | Statistics | Retain advanced |
| `MultipleTestingMethod` | Differential/statistics | Retain advanced |
| `PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION` | Differential/design | Retain advanced |
| `PAIRED_DESIGN_POLICY_FIXED_BLOCK` | Differential/design | Retain advanced |
| `PAIRED_DESIGN_POLICY_REJECT` | Differential/design | Retain advanced |
| `PairedDesignPolicy` | Differential/design | Retain advanced |
| `EnrichmentIdentifierKind` | Enrichment | Retain advanced |
| `EnrichmentMethod` | Enrichment | Retain advanced |
| `EnrichmentOutsideBackgroundPolicy` | Enrichment | Retain advanced |
| `KinaseActivityConfig` | Kinase/activity | Retain advanced |
| `KinaseActivityMethod` | Kinase/activity | Retain advanced |
| `KinaseActivityPValueMethod` | Kinase/activity | Retain advanced |
| `KinaseActivitySsgseaRankingDirection` | Kinase/activity | Retain advanced |
| `KinaseAdaptivePolicy` | Kinase/prediction | Retain advanced |
| `KinaseAttritionPolicy` | Kinase | Retain advanced |
| `KinaseAttritionViolationMode` | Kinase | Retain advanced |
| `KinasePredictionConfig` | Kinase/prediction | Retain advanced |
| `KinasePredictionMode` | Kinase/prediction | Retain advanced |
| `KinaseProfileMissingValueStrategy` | Kinase/scoring | Retain advanced |
| `KinaseReferenceDisplayAmbiguityPolicy` | Kinase/references | Retain advanced |
| `KinaseReliabilityProfile` | Kinase | Retain advanced |
| `KinaseScoringConfig` | Kinase/scoring | Retain advanced |
| `KinaseScoringMode` | Kinase/scoring | Retain advanced |
| `KinaseSiteSequenceConflictPolicy` | Kinase/phosphosites | Retain advanced |
| `LocalisationPolicy` | Evidence | Retain advanced |
| `LocalisationRequirement` | Evidence | Retain advanced |
| `ProfileSelfInclusionPolicy` | Kinase/scoring | Retain advanced |
| `ReferenceContextCompatibilityPolicy` | References | Retain advanced |
| `SignalomeAssignmentPolicy` | Signalome | Retain advanced |
| `SignalomeCandidateScoringPolicy` | Signalome/clustering | Retain advanced |
| `SignalomeClusteringConfig` | Signalome/clustering | Retain advanced |
| `SignalomeConfig` | Signalome | Retain advanced |
| `SignalomeKinaseNetworkPolicy` | Signalome/network | Retain advanced |
| `SignalomeMode` | Signalome | Retain advanced |
| `SignalomeOutputConfig` | Signalome | Retain advanced |
| `SignalomePerformanceConfig` | Signalome | Retain advanced |
| `SignalomeScientificConfig` | Signalome | Retain advanced |
| `SignalomeScorePreconditioningPolicy` | Signalome | Retain advanced |
| `SignalomeValidationConfig` | Signalome | Retain advanced |
| `TechnicalReplicatePolicy` | Differential/design | Retain advanced |
| `KinaseLibraryResource` | Kinase/references | Retain advanced |
| `KinaseLibraryResourceLoadRequest` | Kinase/references | Retain advanced |
| `KinaseLibraryResourceLoader` | Kinase/references adapter | Retain advanced |
| `load_kinase_library_resource` | Kinase/references adapter | Retain advanced |
| `filter_differential_results` | Differential/tables | Retain advanced |
| `rank_differential_results` | Differential/tables | Retain advanced |
| `publish_dataset` | IO/publishers | Retain advanced |
| `publish_kinase_workflow` | IO/publishers | Retain advanced |
| `publish_signalome_workflow` | IO/publishers | Retain advanced |
| `DifferentialModelDiagnostics` | Differential/diagnostics | Retain advanced |
| `QuantificationDepthTrendDiagnostics` | Differential/diagnostics | Retain advanced |
| `KinaseEligibilityReport` | Kinase/diagnostics | Retain advanced |
| `KinaseWorkflowAttritionProvenance` | Kinase/provenance | Retain advanced |
| `KinaseWorkflowCaveat` | Kinase/diagnostics | Retain advanced |
| `KinaseWorkflowPreprocessingAttritionSummary` | Kinase/diagnostics | Retain advanced |
| `KinaseWorkflowScoringAttritionSummary` | Kinase/diagnostics | Retain advanced |
| `KinaseWorkflowSiteAttritionSummary` | Kinase/diagnostics | Retain advanced |
| `SpsDatasetSiteStatistic` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryProvenance` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryResult` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryValidationIssue` | Preprocessing/control sites | Retain advanced |
| `SpsDiscoveryValidationResult` | Preprocessing/control sites | Retain advanced |
| `SpsReferenceDatasetProvenance` | Preprocessing/control sites | Retain advanced |
| `SpsSampleConditionAssignment` | Preprocessing/control sites | Retain advanced |
| `SpsSelectionBoundaryCounts` | Preprocessing/control sites | Retain advanced |
| `SpsSiteStabilityRecord` | Preprocessing/control sites | Retain advanced |

### Unsupported route categories

The following cease to be supported at the API reset:

- `phospy.api.builders`, `phospy.api.datasets`, `phospy.api.configs` and its
  descendants, `phospy.api.enums`, `phospy.api.requests`,
  `phospy.api.results`, and `phospy.api.workflows`;
- `phospy.advanced.configs` and `phospy.advanced.results`;
- `phospy.workflows` and `phospy.workflows.<capability>`;
- compatibility exports under `phospy.tables`;
- direct external imports from owner and implementation packages; and
- every semi-public route previously defined by ADR-0028.

Internal PhosPy code imports owner packages directly. External callers use the
aggregate facades. The same class or function is not supported under multiple
public routes.

### Removal of the semi-public tier

The semi-public category is deleted.

| Former semi-public contract | Target status |
| --- | --- |
| `PreprocessingStageMetadata` and registry helpers | Internal preprocessing composition |
| `ClusterTreeEngine` and `SignalomeClusteringEngine` protocol route | Internal Signalome backend seam |
| Exact-Python clustering facade exports | Internal Signalome implementation and parity surface |
| Rank-fusion scoring helper | Internal Kinase scoring implementation |

If third-party stages or backends are supported later, they require a new
extension API with explicit lifecycle, validation, versioning, and failure
isolation rather than retention of historical test routes.

### Compatibility and deprecation policy

The reset has no runtime compatibility period and does not wait for the former
2.0.0 removal target.

- Remove dynamic `__getattr__` compatibility exports.
- Remove compatibility-only wrapper modules and import aliases.
- Remove deprecated classes, functions, methods, properties, classmethods,
  arguments, values, and constructor behaviours recorded in
  `phospy._deprecations`.
- Remove all 923 retained records and source exemptions after the corresponding
  behaviour is gone.
- Do not emit warnings or add forwarding modules for removed routes.
- Provide concise release migration documentation only.

Compatibility removal is completed and tested as a dedicated pre-migration
change set before domain files move. This is sequencing for review safety, not
a compatibility window. Formerly deprecated inputs may fail after the reset;
supported inputs retain their scientific behaviour.

### Ownership and dependency direction

Public facades depend inward on owner-defined contracts. Owners never import
from `phospy`, `phospy.api`, or `phospy.advanced`. Cross-capability imports use
owner-defined contracts directly. Signalome therefore consumes the owner-defined
`KinaseWorkflowResult`, not its facade route.

### Serialization boundary

Supported persistence is defined by versioned schemas and documented
table/bundle formats. Pickle paths, class `__module__` values, and direct owner
module imports are not supported persistence contracts.

## Superseded and retained decisions

| Earlier ADR | Effect of ADR-0054 |
| --- | --- |
| ADR-0001 | Supersedes its root-name list and supported-route rules. Retains typed request → `run()` → typed result shape, mandatory boundary validation, result objects as data contracts, and explicit public promotion. |
| ADR-0028 | Fully superseded; the semi-public tier is removed. |
| ADR-0031 | Supersedes its export inventories, submodule support, compatibility adapters, warning period, and 2.0.0 removal target. Retains explicit stability classification, exact object-identity re-exports, and the rule that importability alone is not public support. |
| ADR-0010 | Supersedes compatibility-governance requirements that require historical module routes to survive migration. Single ownership and facade-without-behaviour rules remain. |

## Enforcement

1. Exact-set tests assert root, 66-name ordinary, and 124-name advanced
   `__all__` inventories.
2. Identity tests prove facade exports are the owner-defined objects.
3. Negative tests prove removed routes fail rather than warn or redirect.
4. Static tests reject facade `__getattr__` dispatch and implementation imports
   from public facades.
5. Documentation tests allow ordinary imports only from `phospy` and
   `phospy.api`, and specialist imports only from `phospy.advanced`.
6. Behaviour, parity, provenance, dataframe-ownership, and serialization tests
   prove unchanged supported behaviour.

## Consequences

### Positive

- The public model has one ordinary and one specialist comprehensive route.
- The root remains intentionally small.
- Compatibility machinery cannot constrain domain migration.
- Semi-public ambiguity and duplicate supported paths disappear.
- The 192 current aggregate names have explicit post-reset decisions and
  logical owners.

### Negative

- Existing callers using supported submodules or owner paths must update
  imports immediately.
- Deprecated runtime behaviour fails rather than warns.
- `AnalysisReadyPhosphoDataset` and `DifferentialAnalysisWorkflow` stop being
  root imports; `DesignMatrix` and `ContrastMatrix` stop being public imports.
- Documentation, tests, examples, and typing fixtures require a coordinated
  API-reset change.

## Scope boundaries

This ADR does not change scientific algorithms, dataset invariants, domain
ownership, physical file-move order, or the release number assigned to the
cutover.

## Acceptance criteria

This ADR is accepted with the following decisions fixed:

1. only aggregate `phospy`, `phospy.api`, and `phospy.advanced` are supported;
2. root exports exactly the builder, Kinase workflow, and Signalome workflow;
3. `phospy.api` contains the 66 retained ordinary names above;
4. `phospy.advanced` contains the 124 retained specialist names above;
5. all compatibility and semi-public routes and all deprecated runtime
   behaviours are removed without shims;
6. compatibility removal precedes physical migration as a separate change set;
7. owners never import public facades; and
8. supported scientific behaviour remains unchanged.
