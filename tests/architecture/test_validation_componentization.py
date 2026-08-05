from __future__ import annotations

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"

EXTRACTED_RULE_COMPONENTS = (
    SRC_ROOT / "phospy" / "validation" / "workflows" / "differential_design_rules.py",
    SRC_ROOT / "phospy" / "science" / "datasets" / "preprocessing" / "plan_rules.py",
    SRC_ROOT
    / "phospy"
    / "science"
    / "datasets"
    / "preprocessing"
    / "plan_stage_order.py",
    SRC_ROOT
    / "phospy"
    / "science"
    / "datasets"
    / "preprocessing"
    / "plan_config_resolution.py",
)

COMPLEXITY_BUDGETS = {
    (
        SRC_ROOT / "phospy" / "validation" / "workflows" / "differential.py",
        "ExperimentalDesignContractValidator",
        "run",
    ): 4,
    (
        SRC_ROOT / "phospy" / "science" / "datasets" / "preprocessing" / "plan.py",
        "PreprocessingPlan",
        "__post_init__",
    ): 4,
    (
        SRC_ROOT
        / "phospy"
        / "science"
        / "datasets"
        / "preprocessing"
        / "plan_interpreter.py",
        "PreprocessingPlanInterpreter",
        "run",
    ): 4,
    (
        SRC_ROOT
        / "phospy"
        / "validation"
        / "workflows"
        / "differential_design_rules.py",
        "FixedBlockDesignValidator",
        "run",
    ): 24,
    (
        SRC_ROOT
        / "phospy"
        / "science"
        / "datasets"
        / "preprocessing"
        / "plan_stage_order.py",
        "PreprocessingStageOrderPlanner",
        "run",
    ): 18,
    (
        SRC_ROOT
        / "phospy"
        / "science"
        / "datasets"
        / "preprocessing"
        / "plan_rules.py",
        "PreprocessingBatchCorrectionPlanRuleFamily",
        "run",
    ): 18,
}


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported.add(node.module)
    return imported


def _method_node(path: Path, class_name: str, method_name: str) -> ast.FunctionDef:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    return item
    raise AssertionError(f"{class_name}.{method_name} not found in {path}")


def _branch_complexity(node: ast.FunctionDef) -> int:
    complexity = 1
    for child in ast.walk(node):
        if isinstance(
            child,
            (
                ast.If,
                ast.For,
                ast.AsyncFor,
                ast.While,
                ast.Try,
                ast.ExceptHandler,
                ast.Match,
            ),
        ):
            complexity += 1
        elif isinstance(child, ast.BoolOp):
            complexity += max(0, len(child.values) - 1)
        elif isinstance(child, ast.comprehension):
            complexity += len(child.ifs)
    return complexity


def test_extracted_rule_components_do_not_import_higher_level_workflows() -> None:
    offenders = sorted(
        f"{path.relative_to(PROJECT_ROOT)} imports {module}"
        for path in EXTRACTED_RULE_COMPONENTS
        for module in _imported_modules(path)
        if module == "phospy.workflows" or module.startswith("phospy.workflows.")
    )

    assert offenders == []


def test_modified_validator_and_planning_functions_stay_within_complexity_budget() -> (
    None
):
    observed = {
        f"{path.relative_to(PROJECT_ROOT)}::{class_name}.{method_name}": (
            _branch_complexity(_method_node(path, class_name, method_name)),
            budget,
        )
        for (path, class_name, method_name), budget in COMPLEXITY_BUDGETS.items()
    }
    over_budget = {
        name: {"observed": actual, "budget": budget}
        for name, (actual, budget) in observed.items()
        if actual > budget
    }

    assert over_budget == {}
