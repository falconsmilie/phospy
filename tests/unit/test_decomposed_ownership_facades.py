from __future__ import annotations

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_decomposed_public_import_routes_preserve_identity() -> None:
    from phospy.provenance.models import (
        EnvironmentProvenance,
        ReferenceProvenance,
        RunProvenance,
        TableFingerprint,
        TrustedDatasetConstructionAssertions,
    )
    from phospy.provenance.models.environment import (
        EnvironmentProvenance as OwnerEnvironmentProvenance,
    )
    from phospy.provenance.models.references import (
        ReferenceProvenance as OwnerReferenceProvenance,
    )
    from phospy.provenance.models.tables import (
        TableFingerprint as OwnerTableFingerprint,
    )
    from phospy.provenance.models.trusted_assertions import (
        TrustedDatasetConstructionAssertions as OwnerTrustedDatasetConstructionAssertions,
    )
    from phospy.provenance.models.workflows import RunProvenance as OwnerRunProvenance
    from phospy.provenance.serialization import (
        table_fingerprint_to_payload,
        to_payload,
    )
    from phospy.provenance.serialization.tables import (
        table_fingerprint_to_payload as owner_table_fingerprint_to_payload,
    )
    from phospy.provenance.serialization.workflows import to_payload as owner_to_payload
    from phospy.science.datasets.construction.analysis_ready import (
        AnalysisReadyPhosphoDataset as OwnerAnalysisReadyPhosphoDataset,
    )
    from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
    from phospy.science.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingStageResult,
        PreprocessingState,
    )
    from phospy.science.datasets.preprocessing.plan import (
        PreprocessingPlan as OwnerPreprocessingPlan,
    )
    from phospy.science.datasets.preprocessing.results import (
        PreprocessingStageResult as OwnerPreprocessingStageResult,
    )
    from phospy.science.datasets.preprocessing.trace import (
        PreprocessingState as OwnerPreprocessingState,
    )
    from phospy.science.references.validation import (
        load_reference_manifest,
        validate_reference_manifest,
    )
    from phospy.science.references.validation.bundle_semantics import (
        validate_reference_manifest as owner_validate_reference_manifest,
    )
    from phospy.science.references.validation.manifest_schema import (
        load_reference_manifest as owner_load_reference_manifest,
    )

    assert EnvironmentProvenance is OwnerEnvironmentProvenance
    assert ReferenceProvenance is OwnerReferenceProvenance
    assert RunProvenance is OwnerRunProvenance
    assert TableFingerprint is OwnerTableFingerprint
    assert (
        TrustedDatasetConstructionAssertions
        is OwnerTrustedDatasetConstructionAssertions
    )
    assert to_payload is owner_to_payload
    assert table_fingerprint_to_payload is owner_table_fingerprint_to_payload
    assert AnalysisReadyPhosphoDataset is OwnerAnalysisReadyPhosphoDataset
    assert PreprocessingPlan is OwnerPreprocessingPlan
    assert PreprocessingStageResult is OwnerPreprocessingStageResult
    assert PreprocessingState is OwnerPreprocessingState
    assert load_reference_manifest is owner_load_reference_manifest
    assert validate_reference_manifest is owner_validate_reference_manifest


def test_decomposed_compatibility_routes_do_not_contain_implementation_logic() -> None:
    facade_paths = (
        PROJECT_ROOT / "src/phospy/provenance/models/__init__.py",
        PROJECT_ROOT / "src/phospy/provenance/serialization/__init__.py",
        PROJECT_ROOT / "src/phospy/science/references/validation/__init__.py",
        PROJECT_ROOT / "src/phospy/science/datasets/preprocessing/models.py",
        PROJECT_ROOT / "src/phospy/science/datasets/models.py",
    )

    offenders: list[str] = []
    for path in facade_paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:
            if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                offenders.append(f"{path.relative_to(PROJECT_ROOT)}:{node.lineno}")

    assert offenders == []
