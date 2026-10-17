from __future__ import annotations

import ast
import importlib
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_decomposed_owner_import_routes_preserve_supported_identity() -> None:
    from phospy.datasets.models import AnalysisReadyPhosphoDataset
    from phospy.datasets.models import (
        AnalysisReadyPhosphoDataset as OwnerAnalysisReadyPhosphoDataset,
    )
    from phospy.preprocessing.contracts import (
        PreprocessingPlan,
        PreprocessingStageResult,
        PreprocessingState,
    )
    from phospy.preprocessing.contracts.plan import (
        PreprocessingPlan as OwnerPreprocessingPlan,
    )
    from phospy.preprocessing.contracts.state import (
        PreprocessingState as OwnerPreprocessingState,
    )
    from phospy.preprocessing.results.execution import (
        PreprocessingStageResult as OwnerPreprocessingStageResult,
    )
    from phospy.provenance.models import (
        EnvironmentProvenance,
        TableFingerprint,
    )
    from phospy.provenance.models.environment import (
        EnvironmentProvenance as OwnerEnvironmentProvenance,
    )
    from phospy.provenance.models.tables import (
        TableFingerprint as OwnerTableFingerprint,
    )
    from phospy.provenance.serialization import table_fingerprint_to_payload
    from phospy.provenance.serialization.tables import (
        table_fingerprint_to_payload as owner_table_fingerprint_to_payload,
    )
    from phospy.references.validation import (
        load_reference_manifest,
        validate_reference_manifest,
    )
    from phospy.references.validation.bundle_semantics import (
        validate_reference_manifest as owner_validate_reference_manifest,
    )
    from phospy.references.validation.manifest_schema import (
        load_reference_manifest as owner_load_reference_manifest,
    )
    from phospy.run_lineage import RunProvenance, to_payload
    from phospy.run_lineage.models import RunProvenance as OwnerRunProvenance
    from phospy.run_lineage.serialization import to_payload as owner_to_payload

    assert EnvironmentProvenance is OwnerEnvironmentProvenance
    assert RunProvenance is OwnerRunProvenance
    assert to_payload is owner_to_payload
    assert TableFingerprint is OwnerTableFingerprint
    assert table_fingerprint_to_payload is owner_table_fingerprint_to_payload
    assert AnalysisReadyPhosphoDataset is OwnerAnalysisReadyPhosphoDataset
    assert PreprocessingPlan is OwnerPreprocessingPlan
    assert PreprocessingStageResult is OwnerPreprocessingStageResult
    assert PreprocessingState is OwnerPreprocessingState
    assert load_reference_manifest is owner_load_reference_manifest
    assert validate_reference_manifest is owner_validate_reference_manifest


def test_provenance_initialisers_exclude_non_owner_exports() -> None:
    provenance = importlib.import_module("phospy.provenance")
    models = importlib.import_module("phospy.provenance.models")
    serialization = importlib.import_module("phospy.provenance.serialization")

    for name in (
        "DerivedQuantitativeDataProvenance",
        "PreprocessingStageProvenance",
        "ReferenceProvenance",
        "RunProvenance",
        "TrustedDatasetConstructionAssertions",
        "batch_correction_provenance_to_payload",
    ):
        assert not hasattr(provenance, name)
    for name in (
        "KinaseLibraryResourceProvenance",
        "ReferenceContextProtocol",
        "ReferenceProvenance",
        "RunProvenance",
        "TrustedDatasetConstructionAssertions",
    ):
        assert not hasattr(models, name)
    for name in (
        "batch_correction_provenance_to_payload",
        "from_payload",
        "to_payload",
    ):
        assert not hasattr(serialization, name)


def test_decomposed_compatibility_routes_do_not_contain_implementation_logic() -> None:
    facade_paths = (
        PROJECT_ROOT / "src/phospy/provenance/models/__init__.py",
        PROJECT_ROOT / "src/phospy/provenance/serialization/__init__.py",
        PROJECT_ROOT / "src/phospy/references/validation/__init__.py",
        PROJECT_ROOT / "src/phospy/preprocessing/contracts/__init__.py",
        PROJECT_ROOT / "src/phospy/datasets/__init__.py",
    )

    offenders: list[str] = []
    for path in facade_paths:
        offenders.extend(
            _implementation_definitions(
                path.read_text(encoding="utf-8"),
                source_name=path.relative_to(PROJECT_ROOT).as_posix(),
            )
        )

    assert offenders == []


def test_decomposed_facade_guard_rejects_implementation_definitions() -> None:
    assert _implementation_definitions(
        "async def hidden_implementation():\n    return 1\n",
        source_name="synthetic_facade.py",
    ) == ["synthetic_facade.py:1"]


def _implementation_definitions(source: str, *, source_name: str) -> list[str]:
    tree = ast.parse(source, filename=source_name)
    return [
        f"{source_name}:{node.lineno}"
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    ]
