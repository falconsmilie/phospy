from __future__ import annotations

import ast
from collections.abc import Callable, Mapping
from dataclasses import replace
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
import pytest

from phospy.api import DatasetBuildRequest, Organism, ReferenceBundle
from phospy.errors import PhosPyInputError
from phospy.provenance.immutability import FrozenJsonMapping, thaw_json_mapping
from phospy.provenance.models import RunProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
    EmpiricalBayesPriorDiagnostics,
)
from phospy.science.tables.activity import ActivityMatrix
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PUBLIC_PANDAS_CONTAINER_PATHS = (
    "src/phospy/contracts/dataset_build.py",
    "src/phospy/contracts/requests.py",
    "src/phospy/contracts/results/base.py",
    "src/phospy/contracts/results/enrichment.py",
    "src/phospy/contracts/results/kinase.py",
    "src/phospy/contracts/results/signalome.py",
    "src/phospy/frames/table_schema.py",
    "src/phospy/science/activities/diagnostics.py",
    "src/phospy/science/activities/inputs.py",
    "src/phospy/science/activities/results.py",
    "src/phospy/science/activities/semantics.py",
    "src/phospy/science/datasets/construction/analysis_ready.py",
    "src/phospy/science/datasets/imputation_metadata.py",
    "src/phospy/science/datasets/preprocessing/protein_aware_preparation.py",
    "src/phospy/science/datasets/processing_state.py",
    "src/phospy/science/differential/linear_model.py",
    "src/phospy/science/differential/models/design.py",
    "src/phospy/science/differential/models/diagnostics.py",
    "src/phospy/science/differential/models/results.py",
    "src/phospy/science/prediction/models.py",
    "src/phospy/science/references/models.py",
    "src/phospy/science/signalomes/models.py",
    "src/phospy/science/tables/activity.py",
    "src/phospy/science/tables/datasets.py",
    "src/phospy/science/tables/kinase.py",
    "src/phospy/science/tables/references.py",
    "src/phospy/science/tables/signalome/assignments.py",
    "src/phospy/science/tables/signalome/context.py",
    "src/phospy/science/tables/signalome/modules.py",
    "src/phospy/science/tables/signalome/network.py",
)
_PANDAS_NUMPY_FIELD_TOKENS = (
    "DatasetInput",
    "pd.DataFrame",
    "pd.Index",
    "pd.Series",
    "np.ndarray",
    "npt.NDArray",
)


def test_analysis_ready_dataset_equality_hash_and_content_contract() -> None:
    dataset = _dataset()
    assert dataset.provenance is not None
    equivalent = _dataset()
    different_table = _dataset(phospho_values=((9.0, 2.0), (1.1, 2.1)))
    different_provenance = _dataset(
        provenance=replace(dataset.provenance, workflow_name="changed_for_test")
    )

    assert (dataset == dataset) is True
    assert (dataset == equivalent) is False
    assert isinstance(hash(dataset), int)
    assert equivalent not in {dataset: "original"}

    assert dataset.scientifically_equals(equivalent)
    assert not dataset.scientifically_equals(different_table, include_provenance=False)
    assert not dataset.scientifically_equals(different_provenance)
    assert dataset.scientifically_equals(
        different_provenance,
        include_provenance=False,
    )


def test_differential_result_equality_hash_and_content_contract() -> None:
    result = _differential_result(_strict_result_table())
    equivalent = _differential_result(_strict_result_table())
    changed_table = _strict_result_table()
    changed_table.loc[changed_table.index[0], "logFC"] = 3.0
    different_table = _differential_result(changed_table)
    provenance_a = _differential_result(
        _strict_result_table(),
        workflow_provenance={"run_id": "a"},
    )
    provenance_b = _differential_result(
        _strict_result_table(),
        workflow_provenance={"run_id": "b"},
    )

    assert (result == result) is True
    assert (result == equivalent) is False
    assert isinstance(hash(result), int)

    assert result.scientifically_equals(equivalent)
    assert not result.scientifically_equals(different_table, include_provenance=False)
    assert not provenance_a.scientifically_equals(provenance_b)
    assert provenance_a.scientifically_equals(
        provenance_b,
        include_provenance=False,
    )


def test_differential_result_workflow_provenance_is_frozen_json_state() -> None:
    source_provenance: dict[str, object] = {
        "run_id": "a",
        "nested": {"items": [1, {"enabled": True}], "tuple_items": ("x", "y")},
    }
    result = _differential_result(
        _strict_result_table(),
        workflow_provenance=source_provenance,
    )

    assert isinstance(result.workflow_provenance, FrozenJsonMapping)
    nested = result.workflow_provenance["nested"]
    assert isinstance(nested, FrozenJsonMapping)
    assert nested["items"] == (1, FrozenJsonMapping({"enabled": True}))
    assert nested["tuple_items"] == ("x", "y")

    source_provenance["run_id"] = "mutated"
    with pytest.raises(TypeError):
        result.workflow_provenance["run_id"] = "mutated"  # type: ignore[index]
    with pytest.raises(AttributeError):
        nested["items"].append(2)  # type: ignore[attr-defined]

    assert result.workflow_provenance["run_id"] == "a"
    assert thaw_json_mapping(
        result.workflow_provenance,
        field_name="differential_result.workflow_provenance",
    ) == {
        "run_id": "a",
        "nested": {"items": [1, {"enabled": True}], "tuple_items": ["x", "y"]},
    }
    payload = result.to_payload()
    payload_provenance = cast(dict[str, object], payload["workflow_provenance"])
    cast(dict[str, object], payload_provenance["nested"])["items"] = []
    assert cast(FrozenJsonMapping, result.workflow_provenance)["nested"] == nested


def test_differential_result_provenance_scientific_equality_uses_frozen_json() -> None:
    list_provenance = {
        "run_id": "a",
        "nested": {"items": [1, {"threshold": 0.2}], "labels": ["x", "y"]},
    }
    tuple_provenance = {
        "run_id": "a",
        "nested": {"items": (1, {"threshold": 0.2}), "labels": ("x", "y")},
    }
    changed_provenance = {
        "run_id": "a",
        "nested": {"items": [1, {"threshold": 0.3}], "labels": ["x", "y"]},
    }

    list_result = _differential_result(
        _strict_result_table(),
        workflow_provenance=list_provenance,
    )
    tuple_result = _differential_result(
        _strict_result_table(),
        workflow_provenance=tuple_provenance,
    )
    changed_result = _differential_result(
        _strict_result_table(),
        workflow_provenance=changed_provenance,
    )

    assert list_result.scientifically_equals(tuple_result)
    assert not list_result.scientifically_equals(changed_result)
    assert list_result.scientifically_equals(
        changed_result,
        include_provenance=False,
    )


def test_differential_result_payload_provenance_round_trip_is_scientifically_equal() -> (
    None
):
    result = _differential_result(
        _strict_result_table(),
        workflow_provenance={"nested": {"items": (1, {"label": "x"})}},
    )
    payload = result.to_payload()
    reconstructed = _differential_result(
        _strict_result_table(),
        workflow_provenance=cast(
            Mapping[str, object],
            payload["workflow_provenance"],
        ),
    )

    assert result.scientifically_equals(reconstructed)
    assert payload["workflow_provenance"] == {"nested": {"items": [1, {"label": "x"}]}}


class _UnsupportedProvenanceObject:
    pass


@pytest.mark.parametrize(
    ("provenance_factory", "path_fragment", "message_fragment"),
    (
        (
            lambda: {"outer": [{"bad": pd.DataFrame({"x": [1]})}]},
            "differential_result.workflow_provenance.'outer'[0].'bad'",
            "got DataFrame",
        ),
        (
            lambda: {"outer": [{"bad": pd.Series([1])}]},
            "differential_result.workflow_provenance.'outer'[0].'bad'",
            "got Series",
        ),
        (
            lambda: {"outer": [{"bad": pd.Index([1])}]},
            "differential_result.workflow_provenance.'outer'[0].'bad'",
            "got Index",
        ),
        (
            lambda: {"outer": [{"bad": np.array([1])}]},
            "differential_result.workflow_provenance.'outer'[0].'bad'",
            "got ndarray",
        ),
        (
            lambda: {"outer": [{"bad": _UnsupportedProvenanceObject()}]},
            "differential_result.workflow_provenance.'outer'[0].'bad'",
            "got _UnsupportedProvenanceObject",
        ),
        (
            lambda: {"outer": {1: "bad"}},
            "differential_result.workflow_provenance.'outer'",
            "JSON object keys must be strings; got int",
        ),
        (
            lambda: {"outer": [float("nan")]},
            "differential_result.workflow_provenance.'outer'[0]",
            "must be a finite JSON number",
        ),
        (
            lambda: {"outer": [{"bad": np.int64(1)}]},
            "differential_result.workflow_provenance.'outer'[0].'bad'",
            "got int64",
        ),
    ),
)
def test_differential_result_rejects_unsupported_workflow_provenance_values(
    provenance_factory: Callable[[], object],
    path_fragment: str,
    message_fragment: str,
) -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        _differential_result(
            _strict_result_table(),
            workflow_provenance=cast(Mapping[str, object], provenance_factory()),
        )

    message = str(exc_info.value)
    assert path_fragment in message
    assert message_fragment in message


def test_public_table_wrappers_and_requests_use_identity_equality() -> None:
    frame = pd.DataFrame(
        {"kinase_a": [1.0], "kinase_b": [2.0]},
        index=pd.Index(["site_a"], name="site_id"),
    )
    table = ActivityMatrix(frame=frame)
    equivalent_table = ActivityMatrix(frame=frame.copy(deep=True))
    changed_table = ActivityMatrix(frame=frame.assign(kinase_a=[9.0]))

    assert (table == table) is True
    assert (table == equivalent_table) is False
    assert isinstance(hash(table), int)
    assert table.scientifically_equals(equivalent_table)
    assert not table.scientifically_equals(changed_table)

    design_frame = pd.DataFrame({"intercept": [1.0, 1.0], "group_b": [0.0, 1.0]})
    contrast_frame = pd.DataFrame({"B_vs_A": [0.0, 1.0]}, index=design_frame.columns)
    design = DesignMatrix(design_frame)
    contrast = ContrastMatrix(contrast_frame)

    assert design.scientifically_equals(DesignMatrix(design_frame.copy(deep=True)))
    assert contrast.scientifically_equals(
        ContrastMatrix(contrast_frame.copy(deep=True))
    )

    request = DatasetBuildRequest(phospho=frame)
    equivalent_request = DatasetBuildRequest(phospho=frame.copy(deep=True))

    assert (request == request) is True
    assert (request == equivalent_request) is False
    assert isinstance(hash(request), int)


def test_reference_bundle_reconstruction_uses_named_content_comparison() -> None:
    bundle = _reference_bundle()
    reconstructed = ReferenceBundle._from_owned(
        organism=bundle.organism,
        kinase_substrate_map=bundle.kinase_substrate_map_dataframe(),
        site_sequences=bundle.site_sequences_dataframe(),
        provenance=bundle.provenance,
        manifest=bundle.manifest,
    )
    changed_sequences = bundle.site_sequences_dataframe()
    changed_sequences.loc["AKT1;T308;", "site_sequence"] = ("C" * 15) + "T" + ("A" * 15)
    different_bundle = ReferenceBundle(
        organism=bundle.organism,
        kinase_substrate_map=bundle.kinase_substrate_map_dataframe(),
        site_sequences=changed_sequences,
    )

    assert (bundle == bundle) is True
    assert (bundle == reconstructed) is False
    assert isinstance(hash(bundle), int)
    assert bundle.scientifically_equals(reconstructed)
    assert not bundle.scientifically_equals(different_bundle, include_provenance=False)


def test_public_pandas_dataclasses_do_not_use_implicit_dataclass_equality() -> None:
    problems: list[str] = []
    unsafe_hashes: list[str] = []

    for relative_path in _PUBLIC_PANDAS_CONTAINER_PATHS:
        path = _REPO_ROOT / relative_path
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in (item for item in tree.body if isinstance(item, ast.ClassDef)):
            dataclass_decorator = _dataclass_decorator(node)
            if dataclass_decorator is None:
                continue
            if _dataclass_flag(dataclass_decorator, "unsafe_hash") is True:
                unsafe_hashes.append(f"{relative_path}:{node.name}")
            if not _class_contains_pandas_payload(node):
                continue
            if _dataclass_flag(dataclass_decorator, "eq") is False:
                continue
            if _defines_method(node, "__eq__"):
                continue
            problems.append(f"{relative_path}:{node.name}")

    assert problems == []
    assert unsafe_hashes == []


def _dataset(
    *,
    phospho_values: tuple[tuple[float, float], tuple[float, float]] = (
        (1.0, 2.0),
        (1.1, 2.1),
    ),
    provenance: RunProvenance | None = None,
) -> AnalysisReadyPhosphoDataset:
    index = protein_site_key_index(
        protein_identifiers=["MAPK14", "AKT1"],
        sites=["Y182", "T308"],
    )
    phospho = pd.DataFrame(
        {
            "sample_a": [phospho_values[0][0], phospho_values[0][1]],
            "sample_b": [phospho_values[1][0], phospho_values[1][1]],
        },
        index=index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": index.tolist(),
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
            ],
            "protein_id": ["MAPK14", "AKT1"],
        },
        index=index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
        provenance=provenance,
    )


def _strict_result_table() -> pd.DataFrame:
    site_keys = protein_site_key_index(
        protein_identifiers=["MAPK14", "GSK3B"],
        sites=["Y182", "S9"],
    )
    context = site_key_context_columns(site_keys)
    return pd.DataFrame(
        {
            "site_key": site_keys.tolist(),
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;"],
            "organism": context["organism"],
            "protein_namespace": context["protein_namespace"],
            "protein_identifier": context["protein_identifier"],
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "protein_id": ["MAPK14", "GSK3B"],
            "logFC": [1.0, -1.0],
            "t": [2.0, -2.0],
            "P.Value": [0.05, 0.10],
            "adj.P.Val": [0.10, 0.10],
        },
        index=site_keys.copy(),
    )


def _prior_diagnostics(index: pd.Index) -> EmpiricalBayesPriorDiagnostics:
    return EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=1.0,
        base_prior_degrees_of_freedom=10.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
    )


def _differential_result(
    table: pd.DataFrame,
    *,
    workflow_provenance: Mapping[str, object] | None = None,
) -> DifferentialAnalysisResult:
    index = table.index.copy()
    return DifferentialAnalysisResult(
        residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="residual_variance",
        ),
        posterior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="posterior_residual_variance",
        ),
        prior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom_series_value=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
        prior_variance=1.0,
        prior_degrees_of_freedom=10.0,
        residual_degrees_of_freedom=4.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_prior_diagnostics(index),
        mean_variance_trend_diagnostics=None,
        contrast_tables={"B_vs_A": table},
        workflow_provenance=workflow_provenance,
    )


def _reference_bundle() -> ReferenceBundle:
    substrate_sites = ["MAPK14;Y182;", "AKT1;T308;"]
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": substrate_sites,
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    ("A" * 15) + "Y" + ("A" * 15),
                    ("A" * 15) + "T" + ("A" * 15),
                ]
            },
            index=pd.Index(substrate_sites, name="site_id"),
        ),
    )


def _dataclass_decorator(
    node: ast.ClassDef,
) -> ast.expr | None:
    for decorator in node.decorator_list:
        target = decorator.func if isinstance(decorator, ast.Call) else decorator
        if isinstance(target, ast.Name) and target.id == "dataclass":
            return decorator
        if isinstance(target, ast.Attribute) and target.attr == "dataclass":
            return decorator
    return None


def _dataclass_flag(
    decorator: ast.expr,
    flag_name: str,
) -> bool | None:
    if not isinstance(decorator, ast.Call):
        return None
    for keyword in decorator.keywords:
        if keyword.arg == flag_name and isinstance(keyword.value, ast.Constant):
            value = keyword.value.value
            if isinstance(value, bool):
                return value
    return None


def _class_contains_pandas_payload(node: ast.ClassDef) -> bool:
    if any(_base_name(base) in {"TableSchema", "SeriesSchema"} for base in node.bases):
        return True
    for statement in node.body:
        if not isinstance(statement, ast.AnnAssign):
            continue
        annotation = ast.unparse(statement.annotation)
        if any(token in annotation for token in _PANDAS_NUMPY_FIELD_TOKENS):
            return True
    return False


def _base_name(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    return ""


def _defines_method(node: ast.ClassDef, method_name: str) -> bool:
    return any(
        isinstance(statement, ast.FunctionDef) and statement.name == method_name
        for statement in node.body
    )
