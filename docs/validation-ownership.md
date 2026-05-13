# Validation Ownership Map

This map defines primary ownership for high-impact validation invariants.
Higher-level layers may compose owners, but must not duplicate ownership logic.

| Invariant | Owner | Enforcement point | Should not live in | Tests |
| --- | --- | --- | --- | --- |
| DataFrame shape | `phospy.validation.common.dataframes` (`require_dataframe`, `require_non_empty_dataframe`, `require_aligned_dataframe_shape`) | Table schemas, dataset constructors, and workflow boundary validators call shared primitives | Workflow scientific interpreters and scoring logic | `tests/unit/test_table_schemas.py`, `tests/unit/test_workflow_boundary_error.py` |
| Finite numeric intensity values | `phospy.tables.datasets.PhosphoIntensityMatrix` (composing `require_numeric_dataframe` + `require_finite_numeric_dataframe`) | `AnalysisReadyPhosphoDataset.__post_init__` via `PhosphoIntensityMatrix` and `TotalProteinMatrix` | Public request dataclasses and workflow scoring stages | `tests/unit/test_table_schemas.py`, `tests/unit/test_public_contract_errors.py` |
| Unique site/sample labels | Dataset table wrappers (`PhosphoIntensityMatrix`, `SampleMetadataTable`, `TotalProteinMatrix`) using shared uniqueness checks | Dataset construction boundary (`AnalysisReadyPhosphoDataset.__post_init__`) | Workflow-local helper modules | `tests/unit/test_table_schemas.py`, `tests/unit/test_public_contract_dataset.py` |
| Required site metadata | `phospy.tables.datasets.SiteMetadataTable` (`require_columns`) | Dataset construction boundary (`AnalysisReadyPhosphoDataset.__post_init__`) | Signalome/Kinase workflow validators as a primary owner | `tests/unit/test_table_schemas.py`, `tests/unit/test_public_contract_dataset.py` |
| `site_sequence` contract | `SiteMetadataTable` + `phospy.validation.datasets.site_metadata.validate_site_sequence_column` | Dataset construction boundary; preprocessing must derive-or-fail before construction | Downstream workflow scoring/interpreters | `tests/unit/test_table_schemas.py`, `tests/integration/test_dataset_builder_boundary_honesty_integration.py` |
| Intensity scale establishment | `phospy.validation.transformations.state.IntensityScaleStateValidator` | `AnalysisReadyPhosphoDataset.__post_init__` and builder transformation establishment | Generic dataframe validators and workflow request dataclasses | `tests/unit/test_validation_ownership.py`, `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_establishes_intensity_scale_state_via_supported_path` |
| Processing state coherence | `AnalysisReadyPhosphoDataset.__post_init__` (analysis-ready state coherence checks) | Dataset construction boundary (`processing_state` type/coherence checks) | Workflow validators and table wrappers | `tests/integration/test_dataset_builder_integration.py`, `tests/unit/test_public_contract_dataset.py` |
| Design matrix validity | `phospy.validation.workflows.differential.ExperimentalDesignContractValidator` and `DifferentialAnalysisInterpreter` (rank/DoF guards) | Differential workflow validator/interpreter pipeline before execution | Dataset schemas and generic `validation/common` | `tests/unit/test_differential_workflow_components.py`, `tests/unit/test_differential_analysis.py` |
| Contrast validity | `ExperimentalDesignContractValidator` (contrast condition/name validity) and `DifferentialAnalysisInterpreter` (estimability) | Differential workflow validator/interpreter before executor | Dataset constructors and table schemas | `tests/unit/test_differential_workflow_components.py`, `tests/unit/test_differential_analysis.py` |
| Replicate policy | `phospy.workflows.differential.replicates.TechnicalReplicateResolver` | Differential workflow validation stage before design-matrix assembly | Dataset table wrappers and shared generic validators | `tests/integration/test_differential_with_technical_replicates.py`, `tests/unit/test_validator_boundaries.py` |
| Kinase reference compatibility | `phospy.validation.references.compatibility.ReferenceCompatibilityValidator` | Reference resolution/runtime boundary (resolver/workflow entry) | Kinase workflow scientific execution stages | `tests/unit/test_validation_ownership.py`, `tests/unit/test_validator_boundaries.py` |
| Signalome protein identity | `phospy.validation.datasets.site_metadata.enforce_required_non_empty_string_column` + `enforce_site_identity_rows` | Signalome workflow validator (`SignalomeWorkflowValidator`) | Dataset generic dataframe helpers and signalome science modules | `tests/unit/test_validator_boundaries.py::test_signalome_validator_requires_explicit_site_metadata_protein_id_column`, `tests/unit/test_workflow_boundary_error.py` |
| Localisation eligibility | `phospy.validation.datasets.site_metadata.enforce_localisation_requirement` | Workflow validators (`KinaseWorkflowValidator`, `SignalomeWorkflowValidator`) with request-configured policy | Dataset structural wrappers and reference validators | `tests/unit/test_validator_boundaries.py`, `tests/integration/test_localisation_policy_integration.py` |
| Peptide evidence ambiguity | `phospy.science.evidence.multi_site.resolve_observation_site_rows` + `phospy.science.evidence.dataset_resolution.PeptideEvidenceDatasetResolver` (+ request mode/policy gate in `DatasetBuildRequestValidator`) | Dataset build preprocessing lane for `site_resolution_mode='peptide_evidence'` | Table schemas for analysis-ready dataset and downstream workflow validators | `tests/unit/test_dataset_peptide_evidence_resolution.py`, `tests/unit/test_evidence_multi_site_handling.py` |

## Ownership Notes

- Dataset validation remains private: there is no public dataset `validate()` API.
- Workflow validators should compose shared/domain validators instead of duplicating row-level scientific identity checks.
- Public request dataclasses remain transport objects; boundary ownership sits in validators and domain contracts.
- Policy enums are owned by behavioural domain modules:
  - preprocessing: `phospy.science.datasets.preprocessing.policy_models`
  - differential: `phospy.science.differential.policy_models`
  - shared scoring: `phospy.science.scoring.policy_models`
  - shared policy-enum infrastructure: `phospy.policies.policy_base`
  - root module `phospy.policy_models` has been removed.
