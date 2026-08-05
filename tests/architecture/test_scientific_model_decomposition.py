from __future__ import annotations

import ast
from pathlib import Path

import phospy.science.activities.diagnostics as activity_diagnostics
import phospy.science.activities.inputs as activity_inputs
import phospy.science.activities.method_models as activity_method_models
import phospy.science.activities.models as activity_models
import phospy.science.activities.results as activity_results
import phospy.science.datasets.preprocessing.batch_correction as batch_correction
import phospy.science.datasets.preprocessing.batch_correction_engine as batch_engine
import phospy.science.datasets.preprocessing.batch_correction_models as batch_models
import phospy.science.datasets.preprocessing.batch_correction_provenance as batch_provenance
import phospy.science.datasets.preprocessing.batch_correction_provenance_validation as batch_provenance_validation
import phospy.science.references.manifest as reference_manifest
import phospy.science.references.manifest_files as reference_manifest_files
import phospy.science.references.manifest_model as reference_manifest_model
import phospy.science.references.manifest_policy as reference_manifest_policy
import phospy.science.references.redistribution as reference_redistribution
import phospy.science.transformations.models as transformation_models
import phospy.science.transformations.policy as transformation_policy
import phospy.science.transformations.provenance as transformation_provenance
import phospy.science.transformations.scale_state as transformation_scale_state
import phospy.science.transformations.scale_values as transformation_scale_values
import phospy.validation.datasets.batch_correction_controls as validation_controls
import phospy.validation.datasets.batch_correction_provenance as validation_provenance

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"


def test_transformation_models_route_preserves_identity_with_split_owners() -> None:
    assert (
        transformation_models.IntensityScaleKind
        is transformation_policy.IntensityScaleKind
    )
    assert (
        transformation_models.QuantitativeMeaning
        is transformation_policy.QuantitativeMeaning
    )
    assert (
        transformation_models.QuantitativeMeaningTransitionProvenance
        is transformation_provenance.QuantitativeMeaningTransitionProvenance
    )
    assert (
        transformation_models.IntensityScaleEstablishmentProvenance
        is transformation_provenance.IntensityScaleEstablishmentProvenance
    )
    assert (
        transformation_models.MatrixIntensityScaleState
        is transformation_scale_values.MatrixIntensityScaleState
    )
    assert (
        transformation_models.IntensityTransformationEvent
        is transformation_scale_values.IntensityTransformationEvent
    )
    assert (
        transformation_models.IntensityScaleState
        is transformation_scale_state.IntensityScaleState
    )


def test_activity_models_route_preserves_identity_with_split_owners() -> None:
    assert (
        activity_models.ActivityMethodMetadata
        is activity_method_models.ActivityMethodMetadata
    )
    assert (
        activity_models.ActivityMethodSummary
        is activity_method_models.ActivityMethodSummary
    )
    assert activity_models.KinaseActivityInputs is activity_inputs.KinaseActivityInputs
    assert (
        activity_models.PredMatOverlapSummary is activity_inputs.PredMatOverlapSummary
    )
    assert (
        activity_models.ActivityMethodDiagnostics
        is activity_diagnostics.ActivityMethodDiagnostics
    )
    assert activity_models.KinaseActivityResult is activity_results.KinaseActivityResult


def test_reference_manifest_route_preserves_identity_with_split_owners() -> None:
    assert (
        reference_manifest.RedistributionStatus
        is reference_manifest_policy.RedistributionStatus
    )
    assert (
        reference_manifest.RedistributionEvidenceType
        is reference_manifest_policy.RedistributionEvidenceType
    )
    assert (
        reference_manifest.ReferenceFileManifest
        is reference_manifest_files.ReferenceFileManifest
    )
    assert (
        reference_manifest.SequenceWindowDefinition
        is reference_manifest_files.SequenceWindowDefinition
    )
    assert (
        reference_manifest.RedistributionEvidence
        is reference_redistribution.RedistributionEvidence
    )
    assert (
        reference_manifest.ReferenceManifest
        is reference_manifest_model.ReferenceManifest
    )


def test_batch_correction_routes_preserve_identity_with_split_owners() -> None:
    assert batch_correction.BatchCorrectionPolicy is batch_models.BatchCorrectionPolicy
    assert (
        batch_correction.BatchCorrectionDiagnostics
        is batch_models.BatchCorrectionDiagnostics
    )
    assert batch_correction.BatchCorrectionReport is batch_models.BatchCorrectionReport
    assert batch_correction.BatchCorrectionResult is batch_models.BatchCorrectionResult
    assert (
        batch_correction.LinearResidualizeBatchCorrectionEngine
        is batch_engine.LinearResidualizeBatchCorrectionEngine
    )
    assert batch_correction.BatchCorrectionEngine is batch_engine.BatchCorrectionEngine
    assert (
        batch_provenance.validate_applied_native_sps_ruv_correction_provenance
        is batch_provenance_validation.validate_applied_native_sps_ruv_correction_provenance
        is validation_provenance.validate_applied_native_sps_ruv_correction_provenance
    )
    assert (
        batch_provenance.normalize_applied_selected_site_key_rows
        is batch_provenance_validation.normalize_applied_selected_site_key_rows
        is validation_controls.normalize_applied_selected_site_key_rows
    )


def test_legacy_aggregate_routes_do_not_define_model_classes_or_helpers() -> None:
    for relative_path in (
        "phospy/science/transformations/models.py",
        "phospy/science/activities/models.py",
        "phospy/science/references/manifest.py",
        "phospy/science/datasets/preprocessing/batch_correction.py",
    ):
        tree = ast.parse((SRC_ROOT / relative_path).read_text(encoding="utf-8"))
        definitions = [
            node.name
            for node in tree.body
            if isinstance(node, ast.ClassDef | ast.FunctionDef)
        ]
        assert definitions == []


def test_decomposed_scientific_model_classes_have_single_definition() -> None:
    expected_owner_by_class = {
        "ActivityMethodDiagnostics": "phospy/science/activities/diagnostics.py",
        "ActivityMethodMetadata": "phospy/science/activities/method_models.py",
        "ActivityMethodSummary": "phospy/science/activities/method_models.py",
        "BatchCorrectionDiagnostics": (
            "phospy/science/datasets/preprocessing/batch_correction_models.py"
        ),
        "BatchCorrectionPolicy": (
            "phospy/science/datasets/preprocessing/batch_correction_models.py"
        ),
        "BatchCorrectionReport": (
            "phospy/science/datasets/preprocessing/batch_correction_models.py"
        ),
        "BatchCorrectionResult": (
            "phospy/science/datasets/preprocessing/batch_correction_models.py"
        ),
        "IntensityScaleState": "phospy/science/transformations/scale_state.py",
        "IntensityTransformationEvent": (
            "phospy/science/transformations/scale_values.py"
        ),
        "KinaseActivityInputs": "phospy/science/activities/inputs.py",
        "KinaseActivityResult": "phospy/science/activities/results.py",
        "MatrixIntensityScaleState": "phospy/science/transformations/scale_values.py",
        "QuantitativeMeaningScaleRule": "phospy/science/transformations/policy.py",
        "QuantitativeMeaningTransitionProvenance": (
            "phospy/science/transformations/provenance.py"
        ),
        "ReferenceFileManifest": "phospy/science/references/manifest_files.py",
        "ReferenceManifest": "phospy/science/references/manifest_model.py",
        "RedistributionEvidence": "phospy/science/references/redistribution.py",
        "SequenceWindowDefinition": "phospy/science/references/manifest_files.py",
    }
    observed: dict[str, list[str]] = {name: [] for name in expected_owner_by_class}
    for path in (SRC_ROOT / "phospy").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        relative = path.relative_to(SRC_ROOT).as_posix()
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name in observed:
                observed[node.name].append(relative)

    assert observed == {
        class_name: [owner]
        for class_name, owner in sorted(expected_owner_by_class.items())
    }
