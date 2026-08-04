from __future__ import annotations

import ast
from pathlib import Path

import phospy
import phospy.advanced as advanced_api
import phospy.api as public_api
import phospy.science.datasets.preprocessing.batch_correction_provenance as science_provenance
import phospy.validation.datasets.batch_correction as facade
import phospy.validation.datasets.batch_correction_controls as controls
import phospy.validation.datasets.batch_correction_design as design
import phospy.validation.datasets.batch_correction_provenance as provenance

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET_VALIDATION_ROOT = PROJECT_ROOT / "src" / "phospy" / "validation" / "datasets"
FOCUSED_BATCH_CORRECTION_MODULES = (
    DATASET_VALIDATION_ROOT / "batch_correction_controls.py",
    DATASET_VALIDATION_ROOT / "batch_correction_design.py",
    DATASET_VALIDATION_ROOT / "batch_correction_provenance.py",
)


def test_batch_correction_facade_has_no_validator_implementation() -> None:
    tree = ast.parse(
        (DATASET_VALIDATION_ROOT / "batch_correction.py").read_text(encoding="utf-8"),
    )

    implementation_nodes = tuple(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef)
    )

    assert implementation_nodes == ()


def test_batch_correction_facade_reexports_extracted_objects() -> None:
    assert facade.BatchDesignMetadataValidator is design.BatchDesignMetadataValidator
    assert (
        facade.SampleMetadataAlignmentValidator
        is design.SampleMetadataAlignmentValidator
    )
    assert facade.BatchStructureValidator is design.BatchStructureValidator
    assert facade.ConditionStructureValidator is design.ConditionStructureValidator
    assert facade.ReplicateStructureValidator is design.ReplicateStructureValidator
    assert facade.DesignRankValidator is design.DesignRankValidator
    assert (
        facade.BatchCorrectionAdequacyValidator
        is design.BatchCorrectionAdequacyValidator
    )
    assert facade.ResolvedBatchDesignMetadata is design.ResolvedBatchDesignMetadata
    assert (
        facade.normalize_applied_selected_site_key_rows
        is controls.normalize_applied_selected_site_key_rows
    )
    assert (
        facade.validate_applied_native_sps_ruv_correction_provenance
        is provenance.validate_applied_native_sps_ruv_correction_provenance
    )


def test_science_preprocessing_uses_validation_owned_applied_helpers() -> None:
    assert (
        science_provenance.validate_applied_native_sps_ruv_correction_provenance
        is provenance.validate_applied_native_sps_ruv_correction_provenance
    )
    assert (
        science_provenance.normalize_applied_selected_site_key_rows
        is controls.normalize_applied_selected_site_key_rows
    )


def test_batch_correction_validation_modules_do_not_import_workflows() -> None:
    offenders: list[str] = []
    for path in FOCUSED_BATCH_CORRECTION_MODULES:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            module_name = None
            if isinstance(node, ast.ImportFrom):
                module_name = node.module
            elif isinstance(node, ast.Import):
                module_name = ",".join(alias.name for alias in node.names)
            if module_name is None:
                continue
            if "phospy.workflows" in module_name:
                offenders.append(f"{path.name}:{node.lineno}:{module_name}")

    assert offenders == []


def test_batch_correction_validators_remain_out_of_public_exports() -> None:
    symbols = set(facade.__all__)

    for public_module in (phospy, public_api, advanced_api):
        assert symbols.isdisjoint(getattr(public_module, "__all__", ()))
