from __future__ import annotations

import ast
import importlib
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"


@dataclass(frozen=True, slots=True)
class CompatibilityRoute:
    compat_module: str
    canonical_module: str
    symbols: tuple[str, ...]


TABLE_COMPATIBILITY_ROUTES = (
    CompatibilityRoute(
        compat_module="phospy.tables",
        canonical_module="phospy.science.tables",
        symbols=(
            "ActivityCountSeries",
            "ActivityMatrix",
            "ActivityTargetTable",
            "KinasePredictionMatrix",
            "KinaseNetworkCandidateCorrelationsTable",
            "KinaseNetworkEdgesTable",
            "KinaseNetworkNodesTable",
            "KinaseScoreMatrix",
            "KinaseSubstrateContributionTable",
            "KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS",
            "KinaseSubstrateReference",
            "PhosphoIntensityMatrix",
            "filter_differential_results",
            "rank_differential_results",
            "SampleMetadataTable",
            "SignalomeAssignmentsTable",
            "SignalomeModulesTable",
            "SignalomeProteinSiteContext",
            "SignalomeSiteContext",
            "SiteMetadataTable",
            "SiteSequenceReference",
            "TableSchema",
            "TotalProteinMatrix",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.base",
        canonical_module="phospy.frames.table_schema",
        symbols=(
            "TableSchema",
            "ValidationErrorType",
            "require_canonical_label_index",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.activity",
        canonical_module="phospy.science.tables.activity",
        symbols=(
            "ActivityCountMatrix",
            "ActivityCountSeries",
            "ActivityMatrix",
            "ActivityStatisticsTable",
            "ActivityTargetTable",
            "SeriesSchema",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.datasets",
        canonical_module="phospy.science.tables.datasets",
        symbols=(
            "PhosphoIntensityMatrix",
            "SampleMetadataTable",
            "SiteMetadataTable",
            "TotalProteinMatrix",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.differential",
        canonical_module="phospy.science.tables.differential",
        symbols=(
            "ADJUSTED_P_VALUE_COLUMN",
            "LOG_FOLD_CHANGE_COLUMN",
            "RAW_P_VALUE_COLUMN",
            "filter_differential_results",
            "rank_differential_results",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables._differential_validation",
        canonical_module="phospy.science.tables._differential_validation",
        symbols=(
            "require_boolean",
            "require_column_name",
            "require_differential_result_columns",
            "require_na_position",
            "require_non_negative_threshold",
            "require_numeric_result_column",
            "require_probability_threshold",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.kinase",
        canonical_module="phospy.science.tables.kinase",
        symbols=(
            "KINASE_PROFILE_SCORE_DIAGNOSTIC_COLUMNS",
            "KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT",
            "KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_SCORED",
            "KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED",
            "KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS",
            "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES",
            "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE",
            "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_IN_PROFILE_SUPPORT",
            "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED",
            "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NO_SCORE_COLUMN",
            "KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED",
            "KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED",
            "KinasePredictionMatrix",
            "KinaseProfileScoreDiagnosticTable",
            "KinaseScoreMatrix",
            "KinaseSubstrateContributionTable",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.references",
        canonical_module="phospy.science.tables.references",
        symbols=(
            "KinaseSubstrateReference",
            "ProteinAccessionReference",
            "SiteSequenceReference",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.signalome",
        canonical_module="phospy.science.tables.signalome",
        symbols=(
            "KinaseNetworkCandidateCorrelationsTable",
            "KinaseNetworkEdgesTable",
            "KinaseNetworkNodesTable",
            "SignalomeAssignmentsTable",
            "SignalomeModulesTable",
            "SignalomeProteinSiteContext",
            "SignalomeSiteContext",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.signalome.assignments",
        canonical_module="phospy.science.tables.signalome.assignments",
        symbols=("SignalomeAssignmentsTable",),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.signalome.context",
        canonical_module="phospy.science.tables.signalome.context",
        symbols=(
            "SignalomeProteinSiteContext",
            "SignalomeSiteContext",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.signalome.modules",
        canonical_module="phospy.science.tables.signalome.modules",
        symbols=("SignalomeModulesTable",),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.signalome.network",
        canonical_module="phospy.science.tables.signalome.network",
        symbols=(
            "KinaseNetworkCandidateCorrelationsTable",
            "KinaseNetworkEdgesTable",
            "KinaseNetworkNodesTable",
        ),
    ),
    CompatibilityRoute(
        compat_module="phospy.tables.signalome.common",
        canonical_module="phospy.science.tables.signalome.common",
        symbols=(),
    ),
)


def test_table_compatibility_modules_reexport_canonical_symbols_by_identity() -> None:
    for route in TABLE_COMPATIBILITY_ROUTES:
        compat = importlib.import_module(route.compat_module)
        canonical = importlib.import_module(route.canonical_module)

        assert tuple(compat.__all__) == _expected_symbols(route, canonical)
        for symbol_name in compat.__all__:
            assert getattr(compat, symbol_name) is getattr(canonical, symbol_name)


def test_public_table_import_routes_remain_available() -> None:
    namespace: dict[str, object] = {}
    for route in TABLE_COMPATIBILITY_ROUTES:
        for symbol_name in route.symbols:
            exec(f"from {route.compat_module} import {symbol_name}", namespace)
            exec(f"from {route.canonical_module} import {symbol_name}", namespace)
            assert symbol_name in namespace


def test_table_compatibility_modules_do_not_define_classes_or_functions() -> None:
    offenders: list[str] = []
    for path in (PACKAGE_ROOT / "tables").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef):
                relative = path.relative_to(PROJECT_ROOT).as_posix()
                offenders.append(f"{relative}:{node.lineno}: {node.name}")

    assert offenders == []


def test_science_package_does_not_import_table_compatibility_package() -> None:
    offenders: list[str] = []
    for path in (PACKAGE_ROOT / "science").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported, line_number in _absolute_imports(tree):
            if imported == "phospy.tables" or imported.startswith("phospy.tables."):
                relative = path.relative_to(PROJECT_ROOT).as_posix()
                offenders.append(f"{relative}:{line_number}: {imported}")

    assert offenders == []


def test_generic_table_schema_infrastructure_is_defined_once_under_frames() -> None:
    expected = "src/phospy/frames/table_schema.py"
    class_definitions = _definition_paths(ast.ClassDef, "TableSchema")
    function_definitions = _definition_paths(
        ast.FunctionDef,
        "require_canonical_label_index",
    )

    assert class_definitions == [expected]
    assert function_definitions == [expected]
    assert _assigns_name(
        PACKAGE_ROOT / "frames" / "table_schema.py", "ValidationErrorType"
    )
    assert not _table_compatibility_or_science_tables_assign_name("ValidationErrorType")


def _expected_symbols(
    route: CompatibilityRoute,
    canonical: ModuleType,
) -> tuple[str, ...]:
    canonical_all = getattr(canonical, "__all__", None)
    if canonical_all is not None:
        return tuple(canonical_all)
    return route.symbols


def _absolute_imports(tree: ast.AST) -> tuple[tuple[str, int], ...]:
    imported: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend((alias.name, node.lineno) for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            imported.append((node.module, node.lineno))
        elif isinstance(node, ast.Call):
            target = _static_dynamic_import_target(node)
            if target is not None:
                imported.append((target, node.lineno))
    return tuple(imported)


def _static_dynamic_import_target(node: ast.Call) -> str | None:
    is_importlib_call = (
        isinstance(node.func, ast.Attribute) and node.func.attr == "import_module"
    )
    is_dunder_import = isinstance(node.func, ast.Name) and node.func.id == "__import__"
    if not (is_importlib_call or is_dunder_import) or not node.args:
        return None
    first_arg = node.args[0]
    if isinstance(first_arg, ast.Constant) and isinstance(first_arg.value, str):
        return first_arg.value
    return None


def _definition_paths(node_type: type[ast.AST], name: str) -> list[str]:
    paths: list[str] = []
    for path in PACKAGE_ROOT.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(
            isinstance(node, node_type) and node.name == name for node in ast.walk(tree)
        ):
            paths.append(path.relative_to(PROJECT_ROOT).as_posix())
    return sorted(paths)


def _assigns_name(path: Path, name: str) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return any(_node_assigns_name(node, name) for node in ast.walk(tree))


def _table_compatibility_or_science_tables_assign_name(name: str) -> bool:
    roots = (
        PACKAGE_ROOT / "tables",
        PACKAGE_ROOT / "science" / "tables",
    )
    for root in roots:
        for path in root.rglob("*.py"):
            if _assigns_name(path, name):
                return True
    return False


def _node_assigns_name(node: ast.AST, name: str) -> bool:
    if isinstance(node, ast.Assign):
        return any(_target_name(target) == name for target in node.targets)
    if isinstance(node, ast.AnnAssign):
        return _target_name(node.target) == name
    return False


def _target_name(target: ast.AST) -> str | None:
    if isinstance(target, ast.Name):
        return target.id
    return None
