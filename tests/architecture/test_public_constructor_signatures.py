from __future__ import annotations

import inspect

import phospy
import phospy.api as public_api
from phospy.api import AnalysisReadyPhosphoDataset, ReferenceBundle
from phospy.api.results import (
    DifferentialAnalysisResult,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    PhosphositeImportResult,
    SignalomeWorkflowResult,
)

_PUBLIC_OWNERSHIP_BEARING_CLASSES = (
    AnalysisReadyPhosphoDataset,
    PhosphositeImportResult,
    DifferentialAnalysisResult,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
    ReferenceBundle,
)
_FORBIDDEN_TRANSFER_PARAMETERS = frozenset(
    {
        "_assume_owned",
        "assume_owned",
        "_owned",
        "owned",
        "copy",
        "copy_input",
        "copy_inputs",
        "copy_data",
        "_skip_validation",
        "skip_validation",
        "_validate",
        "validate",
    }
)


def _forbidden_constructor_parameters(owner: type[object]) -> list[str]:
    try:
        parameters = inspect.signature(owner).parameters
    except (TypeError, ValueError):
        return []
    return [
        name
        for name in parameters
        if name.startswith("_")
        or name in _FORBIDDEN_TRANSFER_PARAMETERS
        or "assume" in name
        or "ownership" in name
    ]


def test_supported_api_has_no_public_ownership_transfer_parameters() -> None:
    offenders = {
        owner.__name__: _forbidden_constructor_parameters(owner)
        for owner in _PUBLIC_OWNERSHIP_BEARING_CLASSES
    }

    assert offenders == {
        owner.__name__: [] for owner in _PUBLIC_OWNERSHIP_BEARING_CLASSES
    }


def test_exported_classes_have_no_private_ownership_or_validation_controls() -> None:
    offenders: dict[str, list[str]] = {}
    for module_name, module in (("phospy", phospy), ("phospy.api", public_api)):
        for symbol_name in getattr(module, "__all__", ()):
            owner = getattr(module, symbol_name)
            if not inspect.isclass(owner):
                continue
            forbidden = _forbidden_constructor_parameters(owner)
            if forbidden:
                offenders[f"{module_name}.{symbol_name}"] = forbidden

    assert offenders == {}


def test_private_owned_factories_are_not_exported() -> None:
    for module in (phospy, public_api):
        exported = set(getattr(module, "__all__", ()))
        assert "_from_owned" not in exported
        assert not any(name.startswith("_") for name in exported)
        assert not any("owned" in name.lower() for name in exported)
