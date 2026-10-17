from __future__ import annotations

import ast
from pathlib import Path

import phospy.kinase.activity.diagnostics as activity_diagnostics
import phospy.kinase.activity.inputs as activity_inputs
import phospy.kinase.activity.models as activity_method_models
import phospy.kinase.activity.results as activity_results
import phospy.preparation_records.batch_correction as batch_models
import phospy.preprocessing.batch_correction.contracts.results as active_batch_results
import phospy.preprocessing.batch_correction.execution.linear as batch_engine
import phospy.provenance.models as provenance_models
import phospy.provenance.models.batch_correction as batch_provenance_models
import phospy.quantitative.models as transformation_models
import phospy.quantitative.policy as transformation_policy
import phospy.quantitative.provenance as transformation_provenance
import phospy.quantitative.scale_state as transformation_scale_state
import phospy.quantitative.scale_values as transformation_scale_values
import phospy.references.manifest as reference_manifest
import phospy.references.manifest_files as reference_manifest_files
import phospy.references.manifest_model as reference_manifest_model
import phospy.references.manifest_policy as reference_manifest_policy
import phospy.references.redistribution as reference_redistribution

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"


def test_quantitative_models_route_preserves_identity_with_owner_leaves() -> None:
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


def test_activity_models_are_defined_by_their_concrete_owner_leaves() -> None:
    assert activity_method_models.ActivityMethodMetadata.__module__ == (
        activity_method_models.__name__
    )
    assert activity_method_models.ActivityMethodSummary.__module__ == (
        activity_method_models.__name__
    )
    assert activity_inputs.KinaseActivityInputs.__module__ == activity_inputs.__name__
    assert activity_inputs.PredMatOverlapSummary.__module__ == activity_inputs.__name__
    assert activity_diagnostics.ActivityMethodDiagnostics.__module__ == (
        activity_diagnostics.__name__
    )
    assert activity_results.KinaseActivityResult.__module__ == activity_results.__name__


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
    assert batch_models.BatchCorrectionPolicy.__module__ == batch_models.__name__
    assert batch_models.BatchCorrectionDiagnostics.__module__ == batch_models.__name__
    assert batch_models.BatchCorrectionReport.__module__ == batch_models.__name__
    assert active_batch_results.BatchCorrectionResult.__module__ == (
        active_batch_results.__name__
    )
    assert (
        batch_engine.LinearResidualizeBatchCorrectionEngine.__module__
        == batch_engine.__name__
    )
    assert not hasattr(batch_engine, "BatchCorrectionEngine")
    assert (
        provenance_models.BatchCorrectionProvenance
        is batch_provenance_models.BatchCorrectionProvenance
    )


def test_legacy_aggregate_routes_do_not_define_model_classes_or_helpers() -> None:
    for relative_path in (
        "phospy/quantitative/models.py",
        "phospy/references/manifest.py",
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
        "ActivityMethodDiagnostics": "phospy/kinase/activity/diagnostics.py",
        "ActivityMethodMetadata": "phospy/kinase/activity/models.py",
        "ActivityMethodSummary": "phospy/kinase/activity/models.py",
        "BatchCorrectionDiagnostics": (
            "phospy/preparation_records/batch_correction.py"
        ),
        "BatchCorrectionPolicy": ("phospy/preparation_records/batch_correction.py"),
        "BatchCorrectionReport": ("phospy/preparation_records/batch_correction.py"),
        "BatchCorrectionResult": (
            "phospy/preprocessing/batch_correction/contracts/results.py"
        ),
        "BatchCorrectionProvenance": ("phospy/provenance/models/batch_correction.py"),
        "BatchCorrectionRejectedEntity": (
            "phospy/provenance/models/batch_correction.py"
        ),
        "CorrectedPreprocessingOutput": (
            "phospy/preprocessing/batch_correction/contracts/output.py"
        ),
        "IntensityScaleState": "phospy/quantitative/scale_state.py",
        "IntensityTransformationEvent": ("phospy/quantitative/scale_values.py"),
        "KinaseActivityInputs": "phospy/kinase/activity/inputs.py",
        "KinaseActivityResult": "phospy/kinase/activity/results.py",
        "MatrixIntensityScaleState": "phospy/quantitative/scale_values.py",
        "QuantitativeMeaningScaleRule": "phospy/quantitative/policy.py",
        "QuantitativeMeaningTransitionProvenance": (
            "phospy/quantitative/provenance.py"
        ),
        "ReferenceFileManifest": "phospy/references/manifest_files.py",
        "ReferenceManifest": "phospy/references/manifest_model.py",
        "RedistributionEvidence": "phospy/references/redistribution.py",
        "SequenceWindowDefinition": "phospy/references/manifest_files.py",
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
