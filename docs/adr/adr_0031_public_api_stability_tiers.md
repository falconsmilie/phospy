# ADR: Public API Stability Tiers

## Document Control

- **ADR ID:** ADR-0031
- **Title:** Public API Stability Tiers
- **Status:** Accepted
- **Date:** 2026-06-30
- **Decision Type:** Architecture Decision Record

## Context

`phospy.api` had grown into a broad aggregate namespace with 202 exported
names. That made implementation details look stable, increased compatibility
burden, and encouraged users to import validators, diagnostic internals, and
compatibility constants from the same route as workflow entrypoints.

ADR-0001 remains the public API governance baseline. This ADR narrows the
aggregate facade and defines stability tiers for names exported from stable and
advanced package facades.

## Decision

PhosPy now uses three API tiers with separate stable and advanced namespaces:

- Stable public API (`phospy.api`): default user-facing imports for builders, request objects,
  workflow classes, primary result objects, reference-bundle entrypoints,
  example-level enums, and common exception types.
- Advanced supported API (`phospy.advanced`): explicit opt-in imports for specialized
  configuration, advanced reference-resource loading, control-site policy
  values, diagnostic result models, compatibility result aliases, and
  result-table inspection helpers.
- Internal / experimental API: implementation details that are not exported
  from `phospy.api`, including validators, workflow interpreters/executors,
  result assemblers, low-level scoring helpers, private provenance
  serialization, processing-state internals, reference manifest validation
  internals, nested diagnostic models, and compatibility constants.

`phospy.api.__all__` contains only stable supported names. `phospy.advanced.__all__`
contains advanced supported names. Advanced names may still change through
normal deprecation and documentation policy, but their import route no longer
implies stable support. Names classified as internal / experimental are not
stable or advanced exports; when they are needed for PhosPy development, import
them from their owning modules.

Public submodule wildcard surfaces under `phospy.api` are stable-only or
compatibility-only. A name classified as advanced or internal / experimental
must not appear in a stable public submodule `__all__`, and
`phospy.api.datasets` is a stable-only route that exports only
`AnalysisReadyPhosphoDataset`. Dataset preprocessing diagnostics remain
inspectable from returned objects such as `dataset.preprocessing_report` and
workflow result properties when present; the diagnostic model classes are not
stable import targets under `phospy.api`.

Historical advanced imports from `phospy.api`, `phospy.api.configs`,
`phospy.api.requests`, and `phospy.api.results` are retained during migration
as compatibility adapters. They emit `DeprecationWarning` and identify a
replacement import, usually `phospy.advanced` or a `phospy.advanced.*`
submodule. Internal compatibility imports, where retained, warn that they are
unsupported compatibility routes.

Update note (2026-07-22): `AnalysisReadyPhosphoDataset` remains a stable public
result/domain type and import target. Its ordinary public constructor is sealed
and raises immediately; supported creation is through
`AnalysisReadyDatasetBuilder` or the advanced
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` factory.

Update note (2026-07-22, config ownership): Stable public configuration import
routes are preserved through `phospy.api.configs` and `phospy.api`. Advanced
configuration import routes are owned by `phospy.advanced.configs`. Some public
config names are compatibility re-exports of science-owned policy objects; those
routes must re-export the exact object identity. Public transport DTO classes
remain owned by `phospy.contracts.configs` and are translated by workflow
interpreters into resolved execution models before numerical science code runs.

## Dataset Diagnostics Policy

ADR-0031 chooses Option A: PhosPy does not provide a public diagnostic
re-export module for dataset processing-state internals.

Internal dataset processing-state classes, validation-state classes, and
low-level diagnostic DTOs are not part of the stable public API. They must not
be re-exported from `phospy.api` or `phospy.api.datasets`.

Users access diagnostics through stable public result surfaces:

- `AnalysisReadyPhosphoDataset`
- builder reports
- workflow result objects
- provenance records
- documented serialisable report payloads

This avoids creating a second semi-public API tier for implementation details
and prevents users from depending on internal classes that are likely to change
as validation and preprocessing evolve.

### Consequences

- `phospy.api.datasets` exports only stable dataset entrypoints.
- No `phospy.api.diagnostics` or `phospy.api.advanced_datasets` module is
  introduced for processing-state internals.
- Internal diagnostic classes may remain importable from implementation modules,
  but those imports are unsupported.
- Tests must verify that public API submodules do not re-export internal
  processing-state names.

## Compatibility

This is an intentional public-surface reduction. The project has not promised
that every previous aggregate `phospy.api` re-export remains stable. Removed
aggregate names stay available only through their explicit owning modules when
those modules are themselves supported.

Submodule `__all__` lists are also curated public surfaces. Compatibility
attributes may remain on lower-level wrappers during migration, but wildcard
exports and documented examples must not promote names classified here as
internal / experimental.

Scientific table import routes under `phospy.tables` are retained as advanced
supported compatibility routes outside the aggregate `phospy.api` facade. They
must re-export the exact owned objects from `phospy.science.tables` or
`phospy.frames`; they do not own scientific logic, validation behavior, or
duplicate class definitions.

Configuration wrapper routes under `phospy.api.configs.*` are retained only for
stable names or as compatibility adapters during migration. Advanced
configuration imports must use `phospy.advanced.configs`. Compatibility
adapters must forward to the owning contract/science route and preserve
constructor behavior, enum values, defaults, and deliberate object identity for
science-owned policies.

No validators, workflow executors, workflow interpreters, private result
assemblers, or reference manifest validation internals may be promoted into
`phospy.api` without a new contract review.

## Inventory

The implementation source of truth is `src/phospy/_api_inventory.py`, surfaced
by `src/phospy/api/__init__.py` and `src/phospy/advanced/__init__.py`:

- `_STABLE_PUBLIC_API`
- `_ADVANCED_SUPPORTED_API`
- `_INTERNAL_EXPERIMENTAL_API`

The groups below classify every current stable, advanced, and retained
compatibility-audit export name. The implementation inventory is the source
of truth; this ADR records the policy classification for review.

### Stable Public API (61 names)

- AnalysisReadyDatasetBuilder
- AnalysisReadyPhosphoDataset
- ReferenceBundleBuilder
- DatasetBuildRequest
- PhosphositeImporter
- PhosphositeImportRequest
- PhosphositeImportResult
- DatasetLocalisationConfig
- DatasetPreprocessingConfig
- BatchCovariate
- CategoricalCovariate
- ContinuousCovariate
- FixedEffectCovariate
- ExperimentalDesign
- SampleDesignRecord
- Contrast
- DesignMatrix
- ContrastMatrix
- DifferentialAnalysisRequest
- EnrichmentIdentifierSetProvenance
- EnrichmentIdentifierSetSourceType
- EnrichmentSet
- EnrichmentSetCollection
- EnrichmentWorkflowRequest
- EnrichmentConfig
- GeneSetCollection
- PtmSetCollection
- KinaseWorkflowRequest
- SignalomeWorkflowRequest
- all_pairwise_contrasts
- contrasts_vs_control
- DifferentialAnalysisWorkflow
- EnrichmentWorkflow
- KinaseWorkflow
- SignalomeWorkflow
- DifferentialAnalysisResult
- EnrichmentResultRecord
- EnrichmentWorkflowResult
- KinaseActivityResult
- KinasePredictionResult
- KinaseScoringResult
- KinaseWorkflowResult
- ResultCaveat
- SignalomeWorkflowResult
- Organism
- ReferenceBundle
- ReferenceBundleBuildRequest
- ReferencePreset
- IntensityScaleKind
- PhosPyError
- PhosPyInputError
- UnsupportedInputFormatError
- PhosPyValidationError
- ContractValidationError
- PhosPyReferenceError
- ReferenceResolutionError
- ReferenceCompatibilityError
- PhosPyWorkflowError
- WorkflowValidationError
- WorkflowBoundaryError
- SignalomeScaleError

### Advanced Supported API (100 names)

- ControlSiteAnnotation
- ControlSiteSet
- ControlSiteSourceMetadata
- ControlSiteStatus
- CorrectionMaskPolicy
- CorrectionMissingnessPolicy
- ObservationMask
- OriginallyMissingCellTracking
- TemporaryImputationMethod
- TemporaryImputationPolicy
- DatasetBatchCorrectionConfig
- DatasetBatchCorrectionMethod
- DatasetPreprocessingBatchCorrectionConfig
- SpsRuvBatchCorrectionConfig (native PhosPy SPS/RUV-style; not RUV-III)
- SpsRuvBatchCorrectionMethod
- DatasetComparisonBuildingConfig
- DatasetComparisonBuildingPolicy
- DatasetComparisonPair
- DatasetGroupCoverageFilterConfig
- DatasetIntensityTransformConfig
- DatasetIntensityTransformPolicy
- DatasetLocalisationMode
- DatasetMissingDataConfig
- DatasetMissingDataInputScale
- DatasetMissingDataPolicy
- DatasetNormalisationConfig
- DatasetNormalisationPolicy
- DatasetProteinAwarePreparationConfig
- DatasetProteinAwarePreparationMappingPolicy
- DatasetProteinAwarePreparationPolicy
- DatasetRuvReadinessConfig
- DatasetSiteMatrixConfig
- DatasetSiteMatrixDuplicateSitePolicy
- DatasetSiteMatrixMissingDataPolicy
- DatasetSiteMatrixPolicy
- DatasetSiteSequenceConflictPolicy
- DatasetSiteSequenceResolutionConfig
- DatasetSiteSequenceResolutionMode
- DatasetTotalProteinCorrectionConfig
- DatasetTotalProteinCorrectionDuplicatePolicy
- DatasetTotalProteinCorrectionIdentityConfig
- DatasetTotalProteinCorrectionIdentityMatchingPolicy
- DatasetTotalProteinCorrectionIdentityMode
- DatasetTotalProteinCorrectionPolicy
- DatasetTotalProteinCorrectionUnmatchedPolicy
- DifferentialAnalysisConfig
- DifferentialImputedValuePolicy
- DifferentialReliabilityProfile
- EmpiricalBayesConfig
- MultipleTestingConfig
- MultipleTestingCorrection
- MultipleTestingMethod
- PairedDesignPolicy
- EnrichmentIdentifierKind
- EnrichmentMethod
- EnrichmentOutsideBackgroundPolicy
- KinaseActivityConfig
- KinaseActivityMethod
- KinaseActivityPValueMethod
- KinaseActivitySsgseaRankingDirection
- KinaseAdaptivePolicy
- KinaseAttritionPolicy
- KinaseAttritionViolationMode
- KinasePredictionConfig
- KinasePredictionMode
- KinaseProfileMissingValueStrategy
- KinaseReferenceDisplayAmbiguityPolicy
- KinaseReliabilityProfile
- KinaseScoringConfig
- KinaseScoringMode
- KinaseSiteSequenceConflictPolicy
- LocalisationPolicy
- LocalisationRequirement
- ProfileSelfInclusionPolicy
- ReferenceContextCompatibilityPolicy
- SignalomeAssignmentPolicy
- SignalomeCandidateScoringPolicy
- SignalomeClusteringConfig
- SignalomeConfig
- SignalomeKinaseNetworkPolicy
- SignalomeMode
- SignalomeOutputConfig
- SignalomePerformanceConfig
- SignalomeScientificConfig
- SignalomeScorePreconditioningPolicy
- SignalomeValidationConfig
- TechnicalReplicatePolicy
- KinaseLibraryResource
- KinaseLibraryResourceLoadRequest
- KinaseLibraryResourceLoader
- load_kinase_library_resource
- filter_differential_results
- rank_differential_results
- DifferentialModelDiagnostics
- KinaseEligibilityReport
- KinaseWorkflowAttritionProvenance
- KinaseWorkflowCaveat
- KinaseWorkflowPreprocessingAttritionSummary
- KinaseWorkflowScoringAttritionSummary
- KinaseWorkflowSiteAttritionSummary

### Internal / Experimental API (102 names)

- ActivityMethodDiagnostics
- BatchCorrectionDiagnostics
- BatchCorrectionPolicy
- BatchCorrectionReport
- ComparisonState
- DatasetBuildError
- DatasetPreprocessingReport
- DatasetProcessingState
- DatasetValidationError
- DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED
- DifferentialContrastDefinition
- DifferentialDesignMatrixSummary
- DifferentialEmpiricalBayesProvenance
- DifferentialFixedEffectCovariateProvenance
- DifferentialMissingValuePolicyProvenance
- DifferentialPolicyProvenance
- DifferentialReplicatePolicyProvenance
- DifferentialStatisticalTestingProvenance
- DifferentialTechnicalReplicateGroup
- DifferentialUnsupportedDesignPolicyProvenance
- ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID
- ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
- ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE
- ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID
- ENRICHMENT_IDENTIFIER_KIND_SITE_KEY
- ENRICHMENT_METHOD_OVER_REPRESENTATION
- GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS
- IMPUTED_VALUE_POLICY_REJECT
- IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES
- IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
- IMPORTER_QUALITY_STATUS_NOT_REPORTED
- IMPORTER_QUALITY_STATUS_REPORTED
- ImporterDetectedIntensityColumn
- ImporterDuplicateKeySummary
- ImporterFlaggedRowSummary
- ImporterLocalisationConfidenceSummary
- ImporterMissingIntensitySummary
- ImporterQualityCount
- ImporterQualityReport
- ImporterQualityStatus
- ImputationObservationMetadata
- IntensityScaleState
- InvalidTransformationStateError
- KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
- KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
- KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF
- KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF
- KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
- KINASE_SCORING_MODES
- KinaseLibraryMatrix
- KinaseLibraryResidueClass
- KseaZScoreActivityDiagnostics
- MatrixIntensityScaleState
- MissingDataState
- MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
- MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI
- MULTIPLE_TESTING_CORRECTION_BONFERRONI
- MULTIPLE_TESTING_CORRECTION_HOLM
- MULTIPLE_TESTING_CORRECTION_NONE
- MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG
- MULTIPLE_TESTING_METHOD_BENJAMINI_YEKUTIELI
- MULTIPLE_TESTING_METHOD_BONFERRONI
- MULTIPLE_TESTING_METHOD_HOLM
- MULTIPLE_TESTING_METHOD_NONE
- NormalisationState
- PAIRED_DESIGN_POLICY_FIXED_BLOCK
- PAIRED_DESIGN_POLICY_REJECT
- PhosPyBuildError
- PhosPyTransformationError
- ProteinAwareMappingDiagnostics
- ProteinAwarePreparationReport
- ProteinAwarePreparationResult
- ProteinAwareSiteEligibility
- PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS
- QuantitativeMeaning
- ReferenceBundleMissingValueCount
- ReferenceBundleSourceFileValidationReport
- ReferenceBundleTableValidationReport
- ReferenceBundleValidationReport
- ReferenceIdentifierNormalisationValidationError
- ReferenceValidationError
- RuvReadinessState
- SignalomeClusteringEngine
- SiteMatrixState
- SiteSequenceResolutionRowDiagnostic
- SiteSequenceResolutionState
- SsgseaSubstrateEnrichmentActivityDiagnostics
- SUPPORTED_DIFFERENTIAL_IMPUTED_VALUE_POLICIES
- SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS
- SUPPORTED_ENRICHMENT_METHODS
- SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
- SUPPORTED_MULTIPLE_TESTING_METHODS
- SUPPORTED_PAIRED_DESIGN_POLICIES
- TotalProteinCorrectionDiagnostics
- TotalProteinCorrectionDiagnosticsV1
- TotalProteinCorrectionState
- TransformationStateEstablishmentError
- TransformationValidationError
- TransformerExecutionError
- UnsupportedOrganismError
- WeightedSubstrateActivityDiagnostics
- WorkflowStageError

## Review Criteria

Future API changes must satisfy these criteria:

1. `phospy.api.__all__` remains stable-only.
2. README examples use stable public API names.
3. Advanced names are exported from `phospy.advanced` and documented as
   advanced.
4. Validators, executors, interpreters, private assemblers, internal scoring
   helpers, private provenance serialization, and reference manifest validation
   internals are not exported through `phospy.api` or `phospy.advanced`.
5. Stable public submodule `__all__` lists do not include names classified as
   advanced or internal / experimental.
6. Any promotion from internal / experimental to supported public API requires
   code, docs, and tests in the same change.
