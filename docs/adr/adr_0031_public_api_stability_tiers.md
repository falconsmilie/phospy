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
aggregate facade and defines stability tiers for names exported from
`phospy.api`.

## Decision

PhosPy now uses three API tiers:

- Stable public API: default user-facing imports for builders, request objects,
  workflow classes, primary result objects, reference-bundle entrypoints,
  example-level enums, and common exception types.
- Advanced supported API: explicit opt-in imports for specialized
  configuration, advanced reference-resource loading, control-site policy
  values, and result-table inspection helpers.
- Internal / experimental API: implementation details that are not exported
  from `phospy.api`, including validators, workflow interpreters/executors,
  result assemblers, low-level scoring helpers, private provenance
  serialization, processing-state internals, reference manifest validation
  internals, nested diagnostic models, and compatibility constants.

`phospy.api.__all__` is the aggregate public facade and contains only stable and
advanced supported names. Advanced names may still change through normal
deprecation and documentation policy, but they are intentionally public. Names
classified as internal / experimental are not aggregate exports; when they are
needed for PhosPy development, import them from their owning modules.

Public submodule wildcard surfaces under `phospy.api` follow the same
stability-tier inventory. A name classified as internal / experimental must not
appear in a public submodule `__all__`, and `phospy.api.datasets` is a
stable-only route that exports only `AnalysisReadyPhosphoDataset`. Dataset
preprocessing diagnostics remain inspectable from returned objects such as
`dataset.preprocessing_report` and workflow result properties when present; the
diagnostic model classes are not supported import targets under `phospy.api`.

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

No validators, workflow executors, workflow interpreters, private result
assemblers, or reference manifest validation internals may be promoted into
`phospy.api` without a new contract review.

## Inventory

The implementation source of truth is `src/phospy/api/__init__.py`:

- `_STABLE_PUBLIC_API`
- `_ADVANCED_SUPPORTED_API`
- `_INTERNAL_EXPERIMENTAL_API`

The three groups below classify every one of the 202 pre-reduction aggregate
exports.

### Stable Public API

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
- PhosPyReferenceError
- ReferenceResolutionError
- ReferenceCompatibilityError
- PhosPyWorkflowError
- WorkflowValidationError
- WorkflowBoundaryError
- SignalomeScaleError

### Advanced Supported API

- CorrectionMaskPolicy
- CorrectionMissingnessPolicy
- ObservationMask
- OriginallyMissingCellTracking
- TemporaryImputationMethod
- TemporaryImputationPolicy
- ControlSiteAnnotation
- ControlSiteSet
- ControlSiteSourceMetadata
- ControlSiteStatus
- DatasetBatchCorrectionConfig
- SpsRuvBatchCorrectionConfig (native PhosPy SPS/RUV-style; not RUV-III)
- DatasetComparisonBuildingConfig
- DatasetGroupCoverageFilterConfig
- DatasetIntensityTransformConfig
- DatasetMissingDataConfig
- DatasetNormalisationConfig
- DatasetProteinAwarePreparationConfig
- DatasetRuvReadinessConfig
- DatasetSiteMatrixConfig
- DatasetSiteSequenceResolutionConfig
- DatasetTotalProteinCorrectionConfig
- DatasetTotalProteinCorrectionIdentityConfig
- DifferentialAnalysisConfig
- DifferentialImputedValuePolicy
- EmpiricalBayesConfig
- MultipleTestingConfig
- MultipleTestingCorrection
- MultipleTestingMethod
- PairedDesignPolicy
- EnrichmentIdentifierKind
- EnrichmentMethod
- KinaseActivityConfig
- KinasePredictionConfig
- KinaseScoringConfig
- LocalisationRequirement
- SignalomeConfig
- SignalomeScientificConfig
- SignalomeClusteringConfig
- SignalomeValidationConfig
- SignalomeOutputConfig
- SignalomePerformanceConfig
- TechnicalReplicatePolicy
- KinaseLibraryResource
- KinaseLibraryResourceLoadRequest
- KinaseLibraryResourceLoader
- load_kinase_library_resource
- filter_differential_results
- rank_differential_results

### Internal / Experimental API

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
- SiteMatrixState
- SiteSequenceResolutionRowDiagnostic
- SiteSequenceResolutionState
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
- WorkflowStageError

## Review Criteria

Future API changes must satisfy these criteria:

1. `phospy.api.__all__` remains the union of stable and explicitly advanced
   supported names only.
2. README examples use stable public API names.
3. Advanced names are documented as advanced.
4. Validators, executors, interpreters, private assemblers, internal scoring
   helpers, private provenance serialization, and reference manifest validation
   internals are not exported through `phospy.api`.
5. Public submodule `__all__` lists do not include names classified as
   internal / experimental.
6. Any promotion from internal / experimental to supported public API requires
   code, docs, and tests in the same change.
