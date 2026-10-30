from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from tools.architecture.import_graph import (
    ImportGraph,
    build_eager_import_graph,
    build_import_graph,
    physical_package_name,
    strongly_connected_components,
)
from tools.architecture.logical_ownership import (
    architectural_package_unit,
    capability_private_import_problems,
    classify_imports,
    modules_matching_prefixes,
    owner_namespace_prefixes,
)


@pytest.fixture
def nested_owner_tree(tmp_path: Path) -> tuple[Path, Path, Mapping[str, str]]:
    src_root = tmp_path / "src"
    package_root = src_root / "phospy"
    sources = {
        "phospy/__init__.py": "",
        # The future container is inert and separately accounted for. It does
        # not become the owner of any child package.
        "phospy/core/__init__.py": "",
        "phospy/core/errors/__init__.py": "",
        "phospy/core/errors/cycle_a.py": "import phospy.core.errors.cycle_b\n",
        "phospy/core/errors/cycle_b.py": "import phospy.core.errors.cycle_a\n",
        "phospy/core/quantitative/__init__.py": ("import phospy.core.errors\n"),
        "phospy/core/quantitative/_private.py": "_secret = object()\n",
        "phospy/core/frames/__init__.py": (
            "import phospy.core.quantitative\nimport phospy.kinase\n"
        ),
        "phospy/kinase/__init__.py": (
            "from phospy.core.quantitative._private import _secret\n"
        ),
        "phospy/datasets/__init__.py": "",
        "phospy/datasets/builder.py": "import phospy.api\n",
        "phospy/api/__init__.py": "",
        "phospy/api/advanced/__init__.py": "",
    }
    for relative, source in sources.items():
        path = src_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source, encoding="utf-8")

    owners_by_module = {
        "phospy": "public_facades",
        "phospy.core": "namespace_containers",
        "phospy.core.errors": "errors",
        "phospy.core.errors.cycle_a": "errors",
        "phospy.core.errors.cycle_b": "errors",
        "phospy.core.quantitative": "quantitative",
        "phospy.core.quantitative._private": "quantitative",
        "phospy.core.frames": "frames",
        "phospy.kinase": "kinase",
        "phospy.datasets": "datasets",
        "phospy.datasets.builder": "datasets",
        "phospy.api": "public_facades",
        "phospy.api.advanced": "public_facades",
    }
    return package_root, src_root, owners_by_module


def _graph(
    nested_owner_tree: tuple[Path, Path, Mapping[str, str]],
) -> tuple[ImportGraph, Path, Path, Mapping[str, str]]:
    package_root, src_root, owners_by_module = nested_owner_tree
    return (
        build_import_graph(package_root=package_root, src_root=src_root),
        package_root,
        src_root,
        owners_by_module,
    )


def test_nested_package_units_follow_owner_namespaces_not_physical_depth(
    nested_owner_tree: tuple[Path, Path, Mapping[str, str]],
) -> None:
    graph, _package_root, _src_root, owners_by_module = _graph(nested_owner_tree)
    unit_edges = {
        (
            architectural_package_unit(
                record.source_module,
                owners_by_module=owners_by_module,
            ),
            architectural_package_unit(
                record.target,
                owners_by_module=owners_by_module,
            ),
        )
        for record in graph.records
    }

    assert physical_package_name("phospy.core.frames.schema") == "phospy.core"
    assert graph.package_edges >= {("phospy.core", "phospy.kinase")}
    assert ("phospy.core.frames", "phospy.core.quantitative") in unit_edges
    assert ("phospy.core.quantitative", "phospy.core.errors") in unit_edges
    assert (
        architectural_package_unit(
            "phospy.core.frames",
            owners_by_module=owners_by_module,
        )
        == "phospy.core.frames"
    )
    assert (
        architectural_package_unit(
            "phospy.core.quantitative._private",
            owners_by_module=owners_by_module,
        )
        == "phospy.core.quantitative"
    )
    assert (
        architectural_package_unit(
            "phospy.api",
            owners_by_module=owners_by_module,
        )
        == "phospy.api"
    )
    assert (
        architectural_package_unit(
            "phospy.api.advanced",
            owners_by_module=owners_by_module,
        )
        == "phospy.api.advanced"
    )


def test_nested_owner_policy_keeps_allowed_forbidden_and_private_edges_visible(
    nested_owner_tree: tuple[Path, Path, Mapping[str, str]],
    tmp_path: Path,
) -> None:
    graph, _package_root, _src_root, owners_by_module = _graph(nested_owner_tree)
    classified = classify_imports(
        graph,
        project_root=tmp_path,
        owners_by_module=owners_by_module,
    )

    def allowed(source: str, target: str) -> bool:
        return next(
            item.allowed
            for item in classified
            if item.record.source_module == source and item.record.target == target
        )

    assert allowed("phospy.core.quantitative", "phospy.core.errors")
    assert not allowed("phospy.core.frames", "phospy.kinase")
    assert not allowed("phospy.core.frames", "phospy.core.quantitative")
    assert not allowed("phospy.datasets.builder", "phospy.api")
    assert not allowed("phospy.kinase", "phospy.core.quantitative._private")

    private_problems = capability_private_import_problems(
        graph,
        project_root=tmp_path,
        owners_by_module=owners_by_module,
    )
    assert len(private_problems) == 1
    assert "prohibited private owner import" in private_problems[0]
    assert "phospy.core.quantitative._private._secret" in private_problems[0]


def test_runtime_prefixes_distinguish_forbidden_and_allowed_core_siblings(
    nested_owner_tree: tuple[Path, Path, Mapping[str, str]],
) -> None:
    _package_root, _src_root, owners_by_module = nested_owner_tree
    forbidden_prefixes = owner_namespace_prefixes(
        {"quantitative", "kinase"},
        owners_by_module=owners_by_module,
    )

    assert forbidden_prefixes == {
        "phospy.core.quantitative",
        "phospy.kinase",
    }
    assert "phospy.core" not in forbidden_prefixes
    assert owner_namespace_prefixes(
        {"errors"},
        owners_by_module=owners_by_module,
    ) == {"phospy.core.errors"}
    assert modules_matching_prefixes(
        (
            "phospy.core",
            "phospy.core.errors",
            "phospy.core.errors.cycle_a",
            "phospy.core.quantitative",
            "phospy.core.quantitative.scale_state",
            "phospy.kinase",
        ),
        prefixes=forbidden_prefixes,
    ) == (
        "phospy.core.quantitative",
        "phospy.core.quantitative.scale_state",
        "phospy.kinase",
    )

    with pytest.raises(ValueError, match="ambiguous runtime prefix"):
        owner_namespace_prefixes(
            {"namespace_containers"},
            owners_by_module=owners_by_module,
        )


def test_eager_cycle_detection_remains_module_exact_under_nested_owners(
    nested_owner_tree: tuple[Path, Path, Mapping[str, str]],
) -> None:
    package_root, src_root, _owners_by_module = nested_owner_tree
    eager_graph = build_eager_import_graph(
        package_root=package_root,
        src_root=src_root,
    )
    cycles = {
        component
        for component in strongly_connected_components(eager_graph.module_edges)
        if len(component) > 1
    }

    assert cycles == {
        (
            "phospy.core.errors.cycle_a",
            "phospy.core.errors.cycle_b",
        )
    }
