from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = ROOT / "examples"


def test_dataset_builder_example_sets_explicit_localisation_policy() -> None:
    source = (EXAMPLES_DIR / "dataset_builder_demo.py").read_text(encoding="utf-8")

    assert "DatasetBuildRequest(" in source
    assert "organism=Organism.RAT" in source
    assert "preprocessing_config=DatasetPreprocessingConfig(" in source
    assert "DatasetLocalisationConfig(" in source
    assert 'mode="require_threshold"' in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source
    assert "DatasetSiteMatrixConfig" not in source
    assert "TemporaryDirectory" not in source


def test_kinase_example_keeps_bundled_reference_lane_explicit() -> None:
    source = (EXAMPLES_DIR / "kinase_workflow_demo.py").read_text(encoding="utf-8")

    assert "references=ReferencePreset.AUTO" in source
    assert "organism=Organism.RAT" in source
    assert "DatasetLocalisationConfig(" in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source


def test_signalome_example_keeps_explicit_protein_identity_contract() -> None:
    source = (EXAMPLES_DIR / "signalome_workflow_demo.py").read_text(encoding="utf-8")

    assert '"protein_id": ["TSC2", "GSK3A"]' in source
    assert "SignalomeWorkflowRequest(kinase_result=kinase_result)" in source
    assert "SignalomeConfig" not in source
    assert "DatasetLocalisationConfig(" in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source
