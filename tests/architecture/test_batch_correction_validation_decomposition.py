from __future__ import annotations

import ast
from pathlib import Path

import phospy
import phospy.advanced as advanced_api
import phospy.api as public_api
from phospy.preprocessing.batch_correction.provenance.validation import (
    validate_applied_native_sps_ruv_correction_provenance,
)
from phospy.preprocessing.batch_correction.validation.controls import (
    normalize_applied_selected_site_key_rows,
)
from phospy.preprocessing.batch_correction.validation.design import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
    ResolvedBatchDesignMetadata,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CAPABILITY_ROOT = PROJECT_ROOT / "src" / "phospy" / "preprocessing" / "batch_correction"
FOCUSED_BATCH_CORRECTION_MODULES = tuple(
    sorted((CAPABILITY_ROOT / "validation").glob("*.py"))
) + (CAPABILITY_ROOT / "provenance" / "validation.py",)


def test_obsolete_batch_correction_validation_routes_are_absent() -> None:
    old_routes = (
        "src/phospy/validation/datasets/batch_correction.py",
        "src/phospy/validation/datasets/batch_correction_controls.py",
        "src/phospy/validation/datasets/batch_correction_design.py",
        "src/phospy/validation/datasets/batch_correction_provenance.py",
        "src/phospy/validation/workflows/batch_correction",
    )

    assert all(
        not path.is_file() and not tuple(path.rglob("*.py"))
        for route in old_routes
        if (path := PROJECT_ROOT / route)
    )


def test_batch_correction_validators_have_defining_owner_identity() -> None:
    assert BatchCorrectionAdequacyValidator.__module__.endswith(".validation.design")
    assert BatchDesignMetadataValidator.__module__.endswith(".validation.design")
    assert ResolvedBatchDesignMetadata.__module__.endswith(".validation.design")
    assert normalize_applied_selected_site_key_rows.__module__.endswith(
        ".validation.controls"
    )
    assert validate_applied_native_sps_ruv_correction_provenance.__module__.endswith(
        ".provenance.validation"
    )


def test_batch_correction_validation_modules_do_not_import_deferred_owners() -> None:
    offenders: list[str] = []
    forbidden = (
        "phospy.science.datasets",
        "phospy.workflows",
    )
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
            if module_name.startswith(forbidden):
                offenders.append(f"{path.name}:{node.lineno}:{module_name}")

    assert offenders == []


def test_batch_correction_validators_remain_out_of_public_exports() -> None:
    symbols = {
        "BatchCorrectionAdequacyValidator",
        "BatchDesignMetadataValidator",
        "ResolvedBatchDesignMetadata",
        "normalize_applied_selected_site_key_rows",
        "validate_applied_native_sps_ruv_correction_provenance",
    }

    for public_module in (phospy, public_api, advanced_api):
        assert symbols.isdisjoint(getattr(public_module, "__all__", ()))
