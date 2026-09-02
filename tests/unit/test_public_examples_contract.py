from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = ROOT / "examples"


def _read_example(script_name: str) -> str:
    return (EXAMPLES_DIR / script_name).read_text(encoding="utf-8")


def _parse_example(script_name: str) -> ast.Module:
    return ast.parse(_read_example(script_name))


def _call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    return None


def _calls_named(tree: ast.Module, call_name: str) -> tuple[ast.Call, ...]:
    return tuple(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and _call_name(node) == call_name
    )


def _imported_names(tree: ast.Module, module_name: str) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == module_name:
            names.update(alias.name for alias in node.names)
    return names


def _literal_keyword_values(
    calls: tuple[ast.Call, ...],
    keyword_name: str,
) -> set[object]:
    values: set[object] = set()
    for call in calls:
        for keyword in call.keywords:
            if keyword.arg == keyword_name and isinstance(keyword.value, ast.Constant):
                values.add(keyword.value.value)
    return values


def _keyword_names(call: ast.Call) -> set[str]:
    return {keyword.arg for keyword in call.keywords if keyword.arg is not None}


def _has_keyword_call(
    calls: tuple[ast.Call, ...],
    *,
    keyword_name: str,
    called_name: str,
) -> bool:
    for call in calls:
        for keyword in call.keywords:
            value = keyword.value
            if (
                keyword.arg == keyword_name
                and isinstance(value, ast.Call)
                and _call_name(value) == called_name
            ):
                return True
    return False


def test_dataset_builder_example_sets_explicit_localisation_policy() -> None:
    source = _read_example("dataset_builder_demo.py")

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
    source = _read_example("kinase_workflow_demo.py")

    assert "references=ReferencePreset.AUTO" in source
    assert "ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT" in source
    assert "organism=Organism.RAT" in source
    assert "DatasetLocalisationConfig(" in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source


def test_signalome_example_keeps_explicit_protein_identity_contract() -> None:
    source = _read_example("signalome_workflow_demo.py")

    assert '"protein_group_id": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"]' in source
    assert "SignalomeWorkflowRequest(" in source
    assert "kinase_result=kinase_result" in source
    assert "SignalomeConfig.production()" in source
    assert "production_config.validation" in source
    assert "ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT" in source
    assert "DatasetLocalisationConfig(" in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source


def test_protein_aware_differential_example_uses_explicit_opt_in() -> None:
    tree = _parse_example("differential_protein_aware_demo.py")

    assert {
        "DatasetProteinAwarePreparationConfig",
        "DifferentialAnalysisConfig",
        "DifferentialProteinAwareModelConfig",
    } <= _imported_names(tree, "phospy.advanced")
    assert "prepare_model_inputs" in _literal_keyword_values(
        _calls_named(tree, "DatasetProteinAwarePreparationConfig"),
        "policy",
    )
    assert _has_keyword_call(
        _calls_named(tree, "DifferentialAnalysisConfig"),
        keyword_name="protein_aware_model",
        called_name="DifferentialProteinAwareModelConfig",
    )
    assert "protein_covariate_adjusted_moderated_linear_model_v1" in {
        node.value for node in ast.walk(tree) if isinstance(node, ast.Constant)
    }

    request_calls = _calls_named(tree, "DifferentialAnalysisRequest")
    assert request_calls
    for call in request_calls:
        assert _keyword_names(call) == {"dataset", "design", "contrasts", "config"}

    sample_records = _calls_named(tree, "SampleDesignRecord")
    assert len(sample_records) == 6
    assert not any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "split"
        for node in ast.walk(tree)
    )
    assert any(
        isinstance(node, ast.Attribute) and node.attr == "protein_aware_diagnostics"
        for node in ast.walk(tree)
    )
    assert "result_status" in {
        node.value for node in ast.walk(tree) if isinstance(node, ast.Constant)
    }
