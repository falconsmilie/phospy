from __future__ import annotations

import importlib
import inspect
import re
from pathlib import Path
from types import ModuleType

import pytest

import phospy
import phospy.advanced as advanced_api
import phospy.api as ordinary_api
import phospy.evidence.config as evidence_config
import phospy.kinase.activity.config as kinase_activity_config
import phospy.kinase.prediction.config as kinase_prediction_config
import phospy.kinase.scoring.config as kinase_scoring_config
import phospy.preprocessing.batch_correction.contracts.configuration as batch_config
import phospy.preprocessing.batch_correction.contracts.missingness as missingness_config
import phospy.preprocessing.contracts.configuration as preprocessing_config
import phospy.references.config as references_config
import phospy.signalome.config as signalome_config

ROOT = Path(__file__).resolve().parents[3]
ADR = ROOT / "docs" / "adr" / "adr_0054_public_api_reset_and_compatibility_removal.md"
CONFIG_OWNER_MODULES = (
    evidence_config,
    kinase_activity_config,
    kinase_prediction_config,
    kinase_scoring_config,
    batch_config,
    missingness_config,
    preprocessing_config,
    references_config,
    signalome_config,
)


def _adr_root_names() -> tuple[str, ...]:
    text = ADR.read_text(encoding="utf-8")
    section = re.search(
        r"^### Root convenience surface\n(?P<body>.*?)(?=^### )",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    assert section is not None
    return tuple(re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"', section.group("body")))


def _adr_table_names(decision: str) -> tuple[str, ...]:
    names: list[str] = []
    for line in ADR.read_text(encoding="utf-8").splitlines():
        match = re.match(
            rf"^\| `([A-Za-z_][A-Za-z0-9_]*)` \| .* \| {decision}(?: |\|)",
            line,
        )
        if match is not None:
            names.append(match.group(1))
    return tuple(names)


def _wildcard_names(module_name: str) -> set[str]:
    namespace: dict[str, object] = {}
    exec(f"from {module_name} import *", namespace)
    return set(namespace).difference({"__builtins__"})


@pytest.mark.parametrize(
    ("module", "expected", "count"),
    (
        (phospy, _adr_root_names(), 3),
        (ordinary_api, _adr_table_names("Retain ordinary"), 66),
        (advanced_api, _adr_table_names("Retain advanced"), 124),
    ),
)
def test_aggregate_export_sets_match_adr_0054_independently(
    module: ModuleType,
    expected: tuple[str, ...],
    count: int,
) -> None:
    assert len(expected) == count
    assert len(set(expected)) == count
    assert tuple(module.__all__) == expected
    assert _wildcard_names(module.__name__) == set(expected)


@pytest.mark.parametrize("module", (phospy, ordinary_api, advanced_api))
def test_aggregate_exports_are_direct_concrete_identities(module: ModuleType) -> None:
    value_owner_overrides = {
        "DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1": (
            "phospy.differential.policies"
        ),
        "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION": ("phospy.differential.policies"),
        "PAIRED_DESIGN_POLICY_FIXED_BLOCK": "phospy.differential.policies",
        "PAIRED_DESIGN_POLICY_REJECT": "phospy.differential.policies",
        "DifferentialImputedValuePolicy": "phospy.differential.policies",
        "DifferentialReliabilityProfile": "phospy.differential.policies",
        "MultipleTestingMethod": "phospy.differential.policies",
        "PairedDesignPolicy": "phospy.differential.policies",
        "TechnicalReplicatePolicy": "phospy.differential.policies",
        "EnrichmentIdentifierKind": "phospy.enrichment.contracts.vocabulary",
        "EnrichmentMethod": "phospy.enrichment.contracts.vocabulary",
        "EnrichmentOutsideBackgroundPolicy": ("phospy.enrichment.contracts.vocabulary"),
        "MultipleTestingCorrection": "phospy.statistics.multiple_testing",
        "DatasetBatchCorrectionMethod": (
            "phospy.preprocessing.batch_correction.contracts.configuration"
        ),
        "DatasetPreprocessingBatchCorrectionConfig": (
            "phospy.preprocessing.batch_correction.contracts.configuration"
        ),
        "SpsRuvBatchCorrectionMethod": (
            "phospy.preprocessing.batch_correction.contracts.configuration"
        ),
        "SignalomeAssignmentPolicy": "phospy.signalome.config",
        "SignalomeCandidateScoringPolicy": "phospy.signalome.config",
        "SignalomeKinaseNetworkPolicy": "phospy.signalome.config",
        "SignalomeMode": "phospy.signalome.config",
        "SignalomeScorePreconditioningPolicy": "phospy.signalome.config",
    }
    for name in module.__all__:
        value = getattr(module, name)
        config_owner = next(
            (
                owner
                for owner in CONFIG_OWNER_MODULES
                if hasattr(owner, name) and value is getattr(owner, name)
            ),
            None,
        )
        if config_owner is not None:
            continue
        owner_module_name = value_owner_overrides.get(name)
        owner_qualname = name
        if owner_module_name is None:
            owner_module_name = getattr(value, "__module__", None)
            owner_qualname = getattr(value, "__qualname__", name)
        assert owner_module_name is not None, name
        owner_module = importlib.import_module(owner_module_name)
        owner_value = owner_module
        for part in owner_qualname.split("."):
            owner_value = getattr(owner_value, part)
        assert value is owner_value, name


def test_builder_definitions_have_exact_owner_identities() -> None:
    from phospy.api import (
        AnalysisReadyDatasetBuilder,
        ReferenceBundleBuilder,
    )
    from phospy.references.builder import ReferenceBundleBuilder as OwnerBuilder

    assert ordinary_api.ReferenceBundleBuilder is ReferenceBundleBuilder
    assert ReferenceBundleBuilder is OwnerBuilder
    assert tuple(inspect.signature(ReferenceBundleBuilder).parameters) == ()
    assert (
        ordinary_api.AnalysisReadyDatasetBuilder is AnalysisReadyDatasetBuilder
    )  # expires 05B
    assert phospy.AnalysisReadyDatasetBuilder is AnalysisReadyDatasetBuilder


def test_removed_root_and_ordinary_names_fail_direct_imports() -> None:
    assert phospy.__version__
    assert "__version__" not in phospy.__all__
    removed = (
        ("phospy", "AnalysisReadyPhosphoDataset"),
        ("phospy", "DifferentialAnalysisWorkflow"),
        ("phospy.api", "DesignMatrix"),
        ("phospy.api", "ContrastMatrix"),
    )
    for module_name, name in removed:
        module = importlib.import_module(module_name)
        assert not hasattr(module, name)
        with pytest.raises(ImportError):
            exec(f"from {module_name} import {name}", {})


def test_aggregate_facades_have_no_dynamic_export_dispatch() -> None:
    for module in (ordinary_api, advanced_api):
        assert "__getattr__" not in vars(module)
    assert "__getattr__" in vars(phospy)
