from __future__ import annotations

import ast
import json
from pathlib import Path

import phospy
import phospy.api
import phospy.api.advanced
from phospy.core.provenance.models import tables as provenance_tables
from phospy.core.references.provenance import models as reference_models
from phospy.preprocessing.provenance import models as preprocessing_models

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUN_LINEAGE_MODELS = PROJECT_ROOT / "src/phospy/core/run_lineage/models.py"
PUBLIC_API_BASELINE = PROJECT_ROOT / "architecture/public-api-baseline.json"

INTERNAL_VALIDATION_CONTRACTS = {
    "require_preprocessing_stage_sequence",
    "require_reference_context_organism_coherence",
    "require_table_fingerprint_sequence",
}

EXPECTED_OWNER_IMPORTS = {
    "phospy.preprocessing.provenance.models": {
        "PreprocessingStageProvenance",
        "require_preprocessing_stage_sequence",
    },
    "phospy.core.provenance.models.tables": {
        "TableFingerprint",
        "require_table_fingerprint_sequence",
    },
    "phospy.core.references.provenance.models": {
        "ReferenceProvenance",
        "require_reference_context_organism_coherence",
    },
}

SUPERSEDED_PRIVATE_HELPERS = {
    "_required_preprocessing_stage_tuple",
    "_required_table_fingerprint_tuple",
    "_require_run_provenance_reference_context_organism_coherence",
}


def test_run_lineage_imports_exact_owner_validation_contracts() -> None:
    tree = ast.parse(RUN_LINEAGE_MODELS.read_text(encoding="utf-8"))
    observed: dict[str, set[str]] = {}
    for node in tree.body:
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.module not in EXPECTED_OWNER_IMPORTS:
            continue
        observed.setdefault(node.module, set()).update(
            alias.name for alias in node.names
        )

    assert observed == EXPECTED_OWNER_IMPORTS


def test_run_lineage_has_no_private_helper_imports_from_affected_owners() -> None:
    tree = ast.parse(RUN_LINEAGE_MODELS.read_text(encoding="utf-8"))
    affected_prefixes = (
        "phospy.preprocessing",
        "phospy.core.provenance",
        "phospy.core.references",
    )
    private_imports = sorted(
        f"{node.module}:{alias.name}"
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        if node.module is not None and node.module.startswith(affected_prefixes)
        for alias in node.names
        if alias.name.startswith("_")
    )

    assert private_imports == []


def test_superseded_private_helpers_are_removed_without_aliases() -> None:
    owner_modules = (
        preprocessing_models,
        provenance_tables,
        reference_models,
    )

    for helper_name in SUPERSEDED_PRIVATE_HELPERS:
        assert all(not hasattr(module, helper_name) for module in owner_modules)


def test_internal_validation_contracts_stay_out_of_supported_facades() -> None:
    supported_facades = (phospy, phospy.api, phospy.api.advanced)
    baseline = json.loads(PUBLIC_API_BASELINE.read_text(encoding="utf-8"))
    baseline_names = {
        name
        for surface in baseline["surfaces"].values()
        for name in surface["ordered_names"]
    }

    assert INTERNAL_VALIDATION_CONTRACTS.isdisjoint(baseline_names)
    for contract_name in INTERNAL_VALIDATION_CONTRACTS:
        assert all(not hasattr(facade, contract_name) for facade in supported_facades)
