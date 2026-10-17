from __future__ import annotations

import ast
from pathlib import Path

from import_graph_support import build_import_graph
from logical_ownership_policy import CAPABILITY_OWNERS, load_module_owners

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
OWNER_MAP = PROJECT_ROOT / "architecture/module-owners.csv"
TOTAL_PROTEIN_OWNER = PACKAGE_ROOT / "preprocessing/stages/total_protein.py"
TOTAL_PROTEIN_POLICY_NAMES = frozenset(
    {
        "total_protein_evidence_meaning_resolver",
        "total_protein_evidence_caveats",
    }
)


def test_total_protein_evidence_policy_has_one_preprocessing_definition() -> None:
    definitions = {name: [] for name in sorted(TOTAL_PROTEIN_POLICY_NAMES)}
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in definitions:
                definitions[node.name].append(path)

    assert definitions == {
        name: [TOTAL_PROTEIN_OWNER] for name in sorted(TOTAL_PROTEIN_POLICY_NAMES)
    }
    owner_entry = next(
        entry
        for entry in load_module_owners(OWNER_MAP)
        if entry.source_path == TOTAL_PROTEIN_OWNER.relative_to(PROJECT_ROOT).as_posix()
    )
    assert owner_entry.owner == "preprocessing"


def test_quantitative_owner_does_not_reexport_total_protein_policy() -> None:
    offenders: list[str] = []
    for path in sorted((PACKAGE_ROOT / "quantitative").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    bound_name = alias.asname or alias.name
                    if bound_name in TOTAL_PROTEIN_POLICY_NAMES:
                        offenders.append(f"{path}:{node.lineno}:{bound_name}")
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    if (
                        isinstance(target, ast.Name)
                        and target.id in TOTAL_PROTEIN_POLICY_NAMES
                    ):
                        offenders.append(f"{path}:{node.lineno}:{target.id}")
                if any(
                    isinstance(target, ast.Name) and target.id == "__all__"
                    for target in node.targets
                ):
                    exported_names = {
                        item.value
                        for item in ast.walk(node.value)
                        if isinstance(item, ast.Constant)
                        and isinstance(item.value, str)
                    }
                    for name in sorted(exported_names & TOTAL_PROTEIN_POLICY_NAMES):
                        offenders.append(f"{path}:{node.lineno}:__all__:{name}")
            elif (
                isinstance(node, ast.AnnAssign)
                and isinstance(node.target, ast.Name)
                and node.target.id in TOTAL_PROTEIN_POLICY_NAMES
            ):
                offenders.append(f"{path}:{node.lineno}:{node.target.id}")

    assert offenders == []


def test_complete_quantitative_import_graph_has_no_capability_dependency() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    owners_by_module = {
        entry.module: entry.owner for entry in load_module_owners(OWNER_MAP)
    }
    offenders = [
        f"{record.source_path.relative_to(PROJECT_ROOT)}:{record.line} -> "
        f"{record.target} ({owners_by_module[record.target]})"
        for record in graph.records
        if owners_by_module[record.source_module] == "quantitative"
        if owners_by_module[record.target] in CAPABILITY_OWNERS
    ]

    assert offenders == []
