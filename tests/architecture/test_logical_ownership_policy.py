from __future__ import annotations

import ast
import json
import types
import typing
from pathlib import Path

from import_graph_support import (
    ImportGraph,
    ImportRecord,
    build_import_graph,
    package_edges,
)
from logical_ownership_policy import (
    DATASET_PATH_READER_EDGE,
    DATASET_PREPROCESSING_ALLOWED_IMPORTS,
    PERMITTED_LOGICAL_OWNER_EDGES,
    PREPARATION_RECORD_ALLOWED_SYMBOLS,
    REVIEWED_MIXED_MODULES,
    classify_imports,
    is_allowed_import,
    load_module_owners,
    logical_owner_edges,
    logical_owner_graph_problems,
    mixed_module_declaration_problems,
    module_inventory_problems,
)

import phospy
import phospy.advanced
import phospy.api
from phospy.statistics.multiple_testing import MultipleTestingCorrection

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
POLICY_ROOT = PROJECT_ROOT / "architecture"
OWNER_MAP = POLICY_ROOT / "module-owners.csv"
API_BASELINE = POLICY_ROOT / "public-api-baseline.json"


def test_every_live_python_module_has_exactly_one_reviewed_owner() -> None:
    graph = _live_graph()
    live_paths = {
        path.relative_to(PROJECT_ROOT).as_posix(): module
        for module, path in graph.module_paths.items()
    }
    entries = load_module_owners(OWNER_MAP)

    assert module_inventory_problems(live_paths=live_paths, entries=entries) == []


def test_permitted_logical_owner_graph_is_acyclic() -> None:
    assert logical_owner_graph_problems(PERMITTED_LOGICAL_OWNER_EDGES) == []


def test_preparation_record_owner_has_closed_paths_and_readers() -> None:
    expected_leaves = {
        "phospy.preparation_records.batch_correction",
        "phospy.preparation_records.group_evidence",
        "phospy.preparation_records.imputation",
        "phospy.preparation_records.policies",
        "phospy.preparation_records.protein_aware",
        "phospy.preparation_records.report_schema",
        "phospy.preparation_records.validation.sample_groups",
    }
    assert set(PREPARATION_RECORD_ALLOWED_SYMBOLS) == expected_leaves
    assert all(PREPARATION_RECORD_ALLOWED_SYMBOLS.values())

    for target in expected_leaves:
        for source_owner in {"datasets", "preprocessing"}:
            assert is_allowed_import(
                source_owner=source_owner,
                target_owner="preparation_records",
                source_module=f"phospy.{source_owner}.reviewed_leaf",
                target_module=target,
            )

    differential_targets = {
        "phospy.differential.planning.protein_aware_inputs": {
            "phospy.preparation_records.policies",
            "phospy.preparation_records.protein_aware",
        },
        "phospy.differential.result_assembly": {
            "phospy.preparation_records.batch_correction",
        },
    }
    for source, permitted_targets in differential_targets.items():
        for target in expected_leaves:
            assert is_allowed_import(
                source_owner="differential",
                target_owner="preparation_records",
                source_module=source,
                target_module=target,
            ) is (target in permitted_targets)
    for target in expected_leaves:
        assert is_allowed_import(
            source_owner="io",
            target_owner="preparation_records",
            source_module="phospy.io.bundles._shared.processing_state",
            target_module=target,
        ) is (target == "phospy.preparation_records.policies")

    assert not is_allowed_import(
        source_owner="differential",
        target_owner="preparation_records",
        source_module="phospy.differential.unapproved_reader",
        target_module="phospy.preparation_records.policies",
    )
    assert not is_allowed_import(
        source_owner="io",
        target_owner="preparation_records",
        source_module="phospy.io.bundles._signalome.reconstruction",
        target_module="phospy.preparation_records.policies",
    )

    for source_owner in {"kinase", "run_lineage", "signalome"}:
        assert not is_allowed_import(
            source_owner=source_owner,
            target_owner="preparation_records",
            source_module=f"phospy.{source_owner}.unapproved_reader",
            target_module="phospy.preparation_records.policies",
        )


def test_logical_owner_cycle_fixture_is_rejected_within_one_physical_package() -> None:
    module_edges = {
        "phospy.science.dataset_side": frozenset({"phospy.science.preprocessing_side"}),
        "phospy.science.preprocessing_side": frozenset({"phospy.science.dataset_side"}),
    }
    owners_by_module = {
        "phospy.science.dataset_side": "datasets",
        "phospy.science.preprocessing_side": "preprocessing",
    }

    # Both modules occupy phospy.science, so the physical package graph has no
    # cross-package edge and its cycle check cannot see this owner-level cycle.
    assert package_edges(module_edges) == set()
    fixture_edges = logical_owner_edges(
        module_edges,
        owners_by_module=owners_by_module,
    )

    assert logical_owner_graph_problems(fixture_edges) == [
        "logical-owner dependency cycle: datasets, preprocessing"
    ]


def test_unregistered_declaration_in_reviewed_mixed_module_fails() -> None:
    review = REVIEWED_MIXED_MODULES["src/phospy/errors/workflows.py"]
    problems = mixed_module_declaration_problems(
        actual_declarations=review.declarations | {"class:UnreviewedProvenance"},
        review=review,
    )

    assert problems == [
        "unregistered mixed-module declaration: class:UnreviewedProvenance"
    ]


def test_dataset_processing_state_builder_has_authoritative_dataset_owner() -> None:
    entry = next(
        entry
        for entry in load_module_owners(OWNER_MAP)
        if entry.source_path == "src/phospy/datasets/construction/state_assembly.py"
    )

    assert (entry.owner, entry.rationale) == (
        "datasets",
        "Dataset construction lifecycle",
    )


def test_capability_specific_provenance_modules_have_semantic_owners() -> None:
    entries = {
        entry.source_path: entry.owner for entry in load_module_owners(OWNER_MAP)
    }
    expected = {
        "src/phospy/datasets/derived_quantitative.py": "datasets",
        "src/phospy/references/provenance/models.py": "references",
        "src/phospy/datasets/trusted_assertions.py": "datasets",
        "src/phospy/preprocessing/provenance/models.py": "preprocessing",
        "src/phospy/references/provenance/context.py": "references",
        "src/phospy/references/provenance/identifiers.py": "references",
        "src/phospy/provenance/serialization/batch_correction.py": "provenance",
        "src/phospy/references/provenance/serialization.py": "references",
        "src/phospy/preprocessing/provenance/serialization.py": "preprocessing",
        "src/phospy/run_lineage/models.py": "run_lineage",
        "src/phospy/run_lineage/serialization.py": "run_lineage",
        "src/phospy/kinase/references/provenance.py": "kinase",
    }

    assert {path: entries[path] for path in expected} == expected


def test_dataset_preprocessing_edges_have_exact_symbol_permissions() -> None:
    observed: dict[tuple[str, str], frozenset[str]] = {}
    for (
        source_module,
        target_module,
    ), permitted_symbols in DATASET_PREPROCESSING_ALLOWED_IMPORTS.items():
        source_path = SRC_ROOT / Path(*source_module.split(".")).with_suffix(".py")
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        imported = frozenset(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            if node.module == target_module
            for alias in node.names
        )
        observed[source_module, target_module] = imported
        assert imported == permitted_symbols

    assert observed == DATASET_PREPROCESSING_ALLOWED_IMPORTS
    for forbidden_target in (
        "phospy.preprocessing.batch_correction.execution.linear",
        "phospy.preprocessing.batch_correction.workflow",
        "phospy.preprocessing.batch_correction.metadata.resolver",
    ):
        assert not is_allowed_import(
            source_owner="datasets",
            target_owner="preprocessing",
            source_module="phospy.datasets.construction.contracts",
            target_module=forbidden_target,
        )


def test_reader_permission_is_the_exact_module_edge() -> None:
    assert DATASET_PATH_READER_EDGE == (
        "phospy.datasets.builder",
        "phospy.io.readers.dataset_inputs",
    )
    assert is_allowed_import(
        source_owner="datasets",
        target_owner="io",
        source_module=DATASET_PATH_READER_EDGE[0],
        target_module=DATASET_PATH_READER_EDGE[1],
    )
    assert not is_allowed_import(
        source_owner="datasets",
        target_owner="io",
        source_module="phospy.datasets.construction.input_resolution",
        target_module="phospy.io.readers.dataset_inputs",
    )


def test_owner_boundaries_and_omitted_count_kinds_are_locked() -> None:
    owner_entries = load_module_owners(OWNER_MAP)
    owners = {entry.source_path: entry.owner for entry in owner_entries}
    owners_by_module = {entry.module: entry.owner for entry in owner_entries}
    assert owners["src/phospy/evidence/quantification_depth.py"] == "evidence"
    assert owners["src/phospy/references/tables.py"] == "references"
    assert "src/phospy/science/tables/references.py" not in owners
    depth_source = (PACKAGE_ROOT / "evidence/quantification_depth.py").read_text(
        encoding="utf-8"
    )
    for declaration in (
        "QUANTIFICATION_DEPTH_KIND_PSM_COUNT",
        "QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT",
        "SUPPORTED_QUANTIFICATION_DEPTH_KINDS",
        "QuantificationDepthKind",
    ):
        assert declaration in depth_source

    graph = _live_graph()
    forbidden_owner_pairs = {
        ("references", "evidence"),
        ("references", "io"),
        ("references", "kinase"),
        ("evidence", "references"),
        ("evidence", "differential"),
    }
    forbidden = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if (
            owners_by_module[record.source_module],
            owners_by_module[record.target],
        )
        in forbidden_owner_pairs
    ]
    assert forbidden == []


def test_stored_reference_contract_has_one_references_owner() -> None:
    table_path = PACKAGE_ROOT / "references/tables.py"
    declarations = _top_level_declarations(table_path)
    moved = {
        ("class", "KinaseSubstrateReference"),
        ("function", "_classify_duplicate_and_conflicting_pair_records"),
        ("function", "_build_pair_classification_records"),
    }

    assert moved.issubset(declarations)
    assert not (PACKAGE_ROOT / "science/tables/references.py").exists()
    assert "src/phospy/science/tables/references.py" not in REVIEWED_MIXED_MODULES


def test_caveat_helper_imports_are_symbol_exact() -> None:
    kinase_source = (PACKAGE_ROOT / "kinase/caveats.py").read_text(encoding="utf-8")
    signalome_source = (PACKAGE_ROOT / "signalome/caveats.py").read_text(
        encoding="utf-8"
    )
    former_helpers = {
        "build_localisation_policy_details",
        "is_permissive_localisation_requirement",
        "_resolve_localisation_column_name",
        "build_reference_context_compatibility_caveat",
    }
    helper_module = "phospy.workflows.result_caveat_helpers"
    for source in (kinase_source, signalome_source):
        tree = ast.parse(source)
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module == helper_module
            for alias in node.names
        }
        assert former_helpers.isdisjoint(imported)
    assert "build_direct_trusted_dataset_construction_caveat" in kinase_source
    assert helper_module not in signalome_source


def test_taxonomy_is_an_errors_only_leaf_and_old_organism_routes_are_unused() -> None:
    taxonomy_root = PACKAGE_ROOT / "taxonomy"
    taxonomy_import_offenders: list[str] = []
    old_route_offenders: list[str] = []
    old_routes = {
        "phospy.provenance.organisms",
        "phospy.science.references",
        "phospy.references.models",
        "phospy.science.sites.organisms",
    }

    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom) or node.module is None:
                continue
            if path.is_relative_to(taxonomy_root):
                if (
                    node.module.startswith("phospy.")
                    and not node.module.startswith("phospy.errors")
                    and not node.module.startswith("phospy.taxonomy")
                ):
                    taxonomy_import_offenders.append(
                        f"{path.relative_to(PROJECT_ROOT)} -> {node.module}"
                    )
            if node.module in old_routes and any(
                alias.name
                in {"Organism", "normalize_optional_organism", "normalize_organism"}
                for alias in node.names
            ):
                old_route_offenders.append(
                    f"{path.relative_to(PROJECT_ROOT)} -> {node.module}"
                )

    assert taxonomy_import_offenders == []
    assert old_route_offenders == []


def test_current_public_api_matches_baseline() -> None:
    payload = json.loads(API_BASELINE.read_text(encoding="utf-8"))
    surfaces = {
        "phospy": phospy,
        "phospy.api": phospy.api,
        "phospy.advanced": phospy.advanced,
    }
    for module_name, module in surfaces.items():
        baseline = payload["surfaces"][module_name]
        assert list(module.__all__) == baseline["ordered_names"]
        assert len(module.__all__) == baseline["count"]
        assert [
            _runtime_definition_identity(module, name) for name in module.__all__
        ] == baseline["definitions"]


def test_public_api_identity_normalizes_pep_604_union_runtime_type() -> None:
    class Surface:
        UnionAlias = int | str

    assert _runtime_definition_identity(Surface, "UnionAlias") == {
        "name": "UnionAlias",
        "runtime_type": "types.UnionType",
        "kind": "value",
        "value_repr": "int | str",
    }


def test_unmapped_and_multiply_mapped_fixture_modules_fail() -> None:
    live_paths = {
        "src/phospy/a.py": "phospy.a",
        "src/phospy/b.py": "phospy.b",
    }
    entry = load_module_owners(OWNER_MAP)[0]
    fixture_entry = type(entry)(
        source_path="src/phospy/a.py",
        module="phospy.a",
        owner="errors",
        rationale="fixture",
    )
    problems = module_inventory_problems(
        live_paths=live_paths,
        entries=(fixture_entry, fixture_entry),
    )

    assert any(
        problem.startswith("multiply mapped source path") for problem in problems
    )
    assert "unmapped source module: src/phospy/b.py" in problems


def test_contract_path_rules_reject_forbidden_capability_and_adapter_paths() -> None:
    assert not is_allowed_import(
        source_owner="phosphosites",
        target_owner="references",
        source_module="phospy.phosphosites.contracts",
        target_module="phospy.provenance.models.references",
    )
    assert not is_allowed_import(
        source_owner="phosphosites",
        target_owner="kinase",
        source_module="phospy.phosphosites.dataset_identity",
        target_module="phospy.kinase.contracts",
    )
    assert not is_allowed_import(
        source_owner="datasets",
        target_owner="preprocessing",
        source_module="phospy.datasets.models",
        target_module="phospy.preprocessing.results.diagnostics",
    )
    assert is_allowed_import(
        source_owner="datasets",
        target_owner="preprocessing",
        source_module="phospy.datasets.construction.planning",
        target_module="phospy.preprocessing.contracts",
    )
    assert not is_allowed_import(
        source_owner="datasets",
        target_owner="preprocessing",
        source_module="phospy.datasets.construction.planning",
        target_module="phospy.preprocessing.execution",
    )
    assert is_allowed_import(
        source_owner="signalome",
        target_owner="kinase",
        source_module="phospy.signalome.planning",
        target_module="phospy.kinase.contracts",
    )
    assert not is_allowed_import(
        source_owner="signalome",
        target_owner="kinase",
        source_module="phospy.signalome.planning",
        target_module="phospy.kinase.execution",
    )
    assert not is_allowed_import(
        source_owner="references",
        target_owner="kinase",
        source_module="phospy.references.models",
        target_module="phospy.kinase.scoring",
    )
    assert not is_allowed_import(
        source_owner="preprocessing",
        target_owner="differential",
        source_module="phospy.preprocessing.evidence",
        target_module="phospy.differential.quantification_depth",
    )
    assert not is_allowed_import(
        source_owner="kinase",
        target_owner="io",
        source_module="phospy.kinase.results",
        target_module="phospy.io.bundles.kinase",
    )
    assert not is_allowed_import(
        source_owner="kinase",
        target_owner="public_facades",
        source_module="phospy.kinase.public",
        target_module="phospy.api",
    )
    assert not is_allowed_import(
        source_owner="public_facades",
        target_owner="kinase",
        source_module="phospy.api.builders",
        target_module="phospy.science.prediction.scoring",
    )
    contract_descendants = (
        (
            "datasets",
            "preprocessing",
            "phospy.datasets.builders.planning",
            "phospy.preprocessing.results.internal",
        ),
        (
            "differential",
            "datasets",
            "phospy.differential.planning",
            "phospy.datasets.models.internal",
        ),
        (
            "differential",
            "preprocessing",
            "phospy.differential.planning",
            "phospy.preprocessing.contracts.internal",
        ),
        (
            "signalome",
            "kinase",
            "phospy.signalome.planning",
            "phospy.kinase.contracts.internal",
        ),
    )
    for (
        source_owner,
        target_owner,
        source_module,
        target_module,
    ) in contract_descendants:
        assert not is_allowed_import(
            source_owner=source_owner,
            target_owner=target_owner,
            source_module=source_module,
            target_module=target_module,
        )


def test_io_adapter_allow_list_is_exact_and_rejects_internal_lookalikes() -> None:
    assert is_allowed_import(
        source_owner="io",
        target_owner="datasets",
        source_module="phospy.io.bundles.dataset",
        target_module="phospy.datasets.models",
    )
    assert is_allowed_import(
        source_owner="io",
        target_owner="preprocessing",
        source_module="phospy.io.importers.contracts",
        target_module="phospy.preprocessing.contracts.configuration.dataset",
    )
    assert not is_allowed_import(
        source_owner="io",
        target_owner="preprocessing",
        source_module="phospy.io.importers.contracts",
        target_module="phospy.preprocessing.contracts.configuration.total_protein",
    )
    assert not is_allowed_import(
        source_owner="io",
        target_owner="kinase",
        source_module="phospy.io.bundles.kinase",
        target_module="phospy.kinase.scoring.algorithms.models",
    )
    assert not is_allowed_import(
        source_owner="io",
        target_owner="signalome",
        source_module="phospy.io.bundles.signalome",
        target_module="phospy.signalome.clustering.backends.constants",
    )
    assert not is_allowed_import(
        source_owner="io",
        target_owner="datasets",
        source_module="phospy.io.bundles.dataset",
        target_module="phospy.datasets.models.serialization",
    )
    assert not is_allowed_import(
        source_owner="io",
        target_owner="datasets",
        source_module="phospy.io.bundles.dataset",
        target_module="phospy.science.datasets.internal_models",
    )


def test_future_adapter_owner_edges_accept_only_published_paths() -> None:
    allowed = (
        ("public_facades", "taxonomy", "phospy.taxonomy.organisms"),
        (
            "public_facades",
            "statistics",
            "phospy.statistics.multiple_testing",
        ),
        (
            "io",
            "enrichment",
            "phospy.enrichment.contracts.collections",
        ),
        (
            "io",
            "enrichment",
            "phospy.enrichment.contracts.vocabulary",
        ),
        ("io", "evidence", "phospy.evidence.contracts"),
        ("io", "preprocessing", "phospy.preprocessing.results"),
    )
    rejected = (
        ("public_facades", "taxonomy", "phospy.taxonomy.validation"),
        (
            "public_facades",
            "statistics",
            "phospy.statistics.multiple_testing.internal",
        ),
        ("io", "enrichment", "phospy.enrichment.contracts"),
        ("io", "enrichment", "phospy.enrichment.execution"),
        ("io", "evidence", "phospy.evidence.validation"),
        ("io", "preprocessing", "phospy.preprocessing.executor"),
    )

    for source_owner, target_owner, target_module in allowed:
        assert (source_owner, target_owner) in PERMITTED_LOGICAL_OWNER_EDGES
        assert is_allowed_import(
            source_owner=source_owner,
            target_owner=target_owner,
            source_module=(
                "phospy.api" if source_owner == "public_facades" else "phospy.io"
            ),
            target_module=target_module,
        )
    for source_owner, target_owner, target_module in rejected:
        assert not is_allowed_import(
            source_owner=source_owner,
            target_owner=target_owner,
            source_module=(
                "phospy.api" if source_owner == "public_facades" else "phospy.io"
            ),
            target_module=target_module,
        )


def test_public_facade_allow_list_rejects_non_contract_implementation_paths(
    tmp_path: Path,
) -> None:
    source_module = "phospy.api"
    source_path = tmp_path / "src/phospy/api/__init__.py"
    published_contract = "phospy.kinase.results"
    implementation_targets = {
        "phospy.kinase.results.serialization",
        "phospy.science.prediction.scoring",
    }
    targets = {published_contract, *implementation_targets}
    graph = ImportGraph(
        modules=frozenset({source_module, *targets}),
        module_paths={
            source_module: source_path,
            **{
                target: tmp_path / f"src/{target.replace('.', '/')}.py"
                for target in targets
            },
        },
        module_edges={
            source_module: frozenset(targets),
            **{target: frozenset() for target in targets},
        },
        package_edges=frozenset(),
        records=tuple(
            ImportRecord(
                source_module=source_module,
                source_path=source_path,
                target=target,
                line=line,
            )
            for line, target in enumerate(sorted(targets), start=1)
        ),
    )

    classified = classify_imports(
        graph,
        project_root=tmp_path,
        owners_by_module={
            source_module: "public_facades",
            **{target: "kinase" for target in targets},
        },
    )
    by_target = {item.record.target: item for item in classified}

    assert by_target[published_contract].allowed
    assert all(not by_target[target].allowed for target in implementation_targets)


def test_domain_to_facade_fixture_is_forbidden() -> None:
    assert not is_allowed_import(
        source_owner="evidence",
        target_owner="public_facades",
        source_module="phospy.evidence.models",
        target_module="phospy.api",
    )


def test_root_facade_import_forms_are_resolved_and_forbidden(tmp_path: Path) -> None:
    src_root = tmp_path / "src"
    package_root = src_root / "phospy"
    evidence_path = package_root / "evidence/models.py"
    evidence_path.parent.mkdir(parents=True)
    (package_root / "__init__.py").write_text("", encoding="utf-8")
    (package_root / "api.py").write_text("", encoding="utf-8")
    evidence_path.write_text(
        "import phospy\nfrom phospy import api\n",
        encoding="utf-8",
    )
    graph = build_import_graph(package_root=package_root, src_root=src_root)

    classified = classify_imports(
        graph,
        project_root=tmp_path,
        owners_by_module={
            "phospy": "public_facades",
            "phospy.api": "public_facades",
            "phospy.evidence.models": "evidence",
        },
    )
    evidence_imports = [
        item
        for item in classified
        if item.record.source_module == "phospy.evidence.models"
    ]

    assert {item.record.target for item in evidence_imports} == {
        "phospy",
        "phospy.api",
    }
    assert all(not item.allowed for item in evidence_imports)


def test_importlib_callable_aliases_resolve_to_forbidden_facade_imports(
    tmp_path: Path,
) -> None:
    sources = {
        "direct": 'from importlib import import_module\nimport_module("phospy.api")\n',
        "aliased": (
            "from importlib import import_module as load_module\n"
            'load_module("phospy.api")\n'
        ),
    }

    for variant, source in sources.items():
        src_root = tmp_path / variant / "src"
        package_root = src_root / "phospy"
        evidence_path = package_root / "evidence/models.py"
        evidence_path.parent.mkdir(parents=True)
        (package_root / "__init__.py").write_text("", encoding="utf-8")
        (package_root / "api.py").write_text("", encoding="utf-8")
        evidence_path.write_text(source, encoding="utf-8")
        graph = build_import_graph(package_root=package_root, src_root=src_root)

        classified = classify_imports(
            graph,
            project_root=tmp_path / variant,
            owners_by_module={
                "phospy": "public_facades",
                "phospy.api": "public_facades",
                "phospy.evidence.models": "evidence",
            },
        )
        facade_imports = [
            item
            for item in classified
            if item.record.source_module == "phospy.evidence.models"
        ]

        assert [item.record.target for item in facade_imports] == ["phospy.api"]
        assert not facade_imports[0].allowed


def _live_graph() -> ImportGraph:
    return build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)


def _top_level_declarations(path: Path) -> set[tuple[str, str]]:
    tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
    declarations: set[tuple[str, str]] = set()
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            declarations.add(("class", node.name))
        elif isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
            declarations.add(("function", node.name))
        elif isinstance(node, ast.Assign):
            declarations.update(
                ("constant", target.id)
                for target in node.targets
                if isinstance(target, ast.Name)
            )
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            declarations.add(("constant", node.target.id))
    return declarations


def _runtime_definition_identity(module: object, name: str) -> dict[str, object]:
    value = getattr(module, name)
    if _is_union(value):
        return {
            "name": name,
            "runtime_type": "types.UnionType",
            "kind": "value",
            "value_repr": repr(value),
        }
    definition_module = getattr(value, "__module__", None)
    qualname = getattr(value, "__qualname__", None)
    result: dict[str, object] = {
        "name": name,
        "runtime_type": f"{type(value).__module__}.{type(value).__qualname__}",
    }
    if module is phospy.advanced and name == "MultipleTestingCorrection":
        assert value is MultipleTestingCorrection
        result["definition"] = (
            "phospy.statistics.multiple_testing:MultipleTestingCorrection"
        )
        result["kind"] = "callable_or_type"
        return result
    if isinstance(definition_module, str) and isinstance(qualname, str):
        result["definition"] = f"{definition_module}:{qualname}"
        result["kind"] = "callable_or_type"
    else:
        result["kind"] = "value"
        result["value_repr"] = repr(value)
    return result


def _is_union(value: object) -> bool:
    origin = typing.get_origin(value)
    return origin is typing.Union or origin is types.UnionType


def _format_import(record: ImportRecord, source_owner: str, target_owner: str) -> str:
    return (
        f"{record.source_path.relative_to(PROJECT_ROOT)}:{record.line}: "
        f"{record.source_module} ({source_owner}) -> {record.target} ({target_owner})"
    )
