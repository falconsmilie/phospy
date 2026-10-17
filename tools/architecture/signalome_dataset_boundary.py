"""Static D08-1 guard for Signalome access to Dataset-owned definitions."""

from __future__ import annotations

import ast
import importlib.util
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

DatasetImport = tuple[str, str, str]

APPROVED_DATASET_IMPORTS: frozenset[DatasetImport] = frozenset(
    {
        (
            "phospy.signalome.results",
            "phospy.datasets.models",
            "AnalysisReadyPhosphoDataset",
        ),
        (
            "phospy.signalome.results",
            "phospy.datasets.processing_state",
            "DatasetPreprocessingReport",
        ),
        (
            "phospy.signalome.result_validation",
            "phospy.datasets.models",
            "AnalysisReadyPhosphoDataset",
        ),
        (
            "phospy.signalome.result_validation",
            "phospy.datasets.internal_view",
            "DatasetInternalView",
        ),
        (
            "phospy.signalome.contracts",
            "phospy.datasets.models",
            "AnalysisReadyPhosphoDataset",
        ),
        (
            "phospy.signalome.contracts",
            "phospy.datasets.internal_view",
            "DatasetInternalView",
        ),
        (
            "phospy.signalome.validator",
            "phospy.datasets.internal_view",
            "DatasetInternalView",
        ),
    }
)

_DATASET_TYPES = frozenset(
    {
        "AnalysisReadyPhosphoDataset",
        "DatasetInternalFrameStore",
        "DatasetInternalView",
    }
)
_PROTECTED_MEMBERS = frozenset(
    {
        "_phospho",
        "_site_metadata",
        "_sample_metadata",
        "_total",
        "_comparisons",
        "_internal_frame_store",
        "_internal_frame_store_for_current_frames",
    }
)
_DATASET_SENSITIVE_NAMES = (
    _DATASET_TYPES | frozenset({"DatasetPreprocessingReport"}) | _PROTECTED_MEMBERS
)
_ACQUISITION_CALLABLES = frozenset(
    {
        "builtins.__import__",
        "builtins.getattr",
        "builtins.hasattr",
        "builtins.vars",
        "importlib.import_module",
        "object.__getattribute__",
        "operator.attrgetter",
    }
)


@dataclass(frozen=True, slots=True)
class DatasetBoundaryViolation:
    category: str
    consumer: str
    line: int
    detail: str


@dataclass(frozen=True, slots=True)
class DatasetBoundaryReport:
    visited_modules: tuple[str, ...]
    observed_imports: tuple[DatasetImport, ...]
    violations: tuple[DatasetBoundaryViolation, ...]

    @property
    def import_counts(self) -> Counter[DatasetImport]:
        return Counter(self.observed_imports)


@dataclass(frozen=True, slots=True)
class _Origin:
    kind: str
    module: str | None = None
    symbol: str | None = None


_UNKNOWN = _Origin("unknown")
_NON_DATASET = _Origin("non_dataset")
_DATASET_MODULE = _Origin("dataset_module")
_DATASET_OBJECT = _Origin("dataset_object")
_PHOSPY_MODULE = _Origin("phospy_module")
_UNKNOWN_MODULE = _Origin("unknown_module")


@dataclass(slots=True)
class _ModuleInfo:
    tree: ast.Module
    definitions: set[str]
    imports: dict[str, tuple[str, str | None, int, bool]]
    all_names: tuple[str, ...] | None
    dynamic_all: bool


class _Index:
    def __init__(self, sources: Mapping[str, str]) -> None:
        self.sources = dict(sources)
        self.modules: dict[str, _ModuleInfo] = {}
        self._resolve_cache: dict[tuple[str, str], _Origin] = {}
        for module, source in self.sources.items():
            tree = ast.parse(source, filename=module)
            definitions: set[str] = set()
            imports: dict[str, tuple[str, str | None, int, bool]] = {}
            values: dict[str, tuple[str, ...]] = {}
            all_names: tuple[str, ...] | None = None
            dynamic_all = False
            for node in tree.body:
                if isinstance(
                    node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
                ):
                    definitions.add(node.name)
                elif isinstance(node, ast.Import):
                    for alias in node.names:
                        bound = alias.asname or alias.name.split(".")[0]
                        imports[bound] = (alias.name, None, 0, alias.asname is not None)
                elif isinstance(node, ast.ImportFrom):
                    resolved = _resolve_from_module(module, node.module, node.level)
                    for alias in node.names:
                        if alias.name == "*":
                            continue
                        imports[alias.asname or alias.name] = (
                            resolved,
                            alias.name,
                            node.level,
                            alias.asname is not None,
                        )
                elif isinstance(node, (ast.Assign, ast.AnnAssign)):
                    targets = (
                        node.targets if isinstance(node, ast.Assign) else [node.target]
                    )
                    value = node.value
                    literal = _literal_name_sequence(value, values)
                    for target in targets:
                        if isinstance(target, ast.Name):
                            definitions.add(target.id)
                            if literal is not None:
                                values[target.id] = literal
                            if target.id == "__all__":
                                if literal is None:
                                    dynamic_all = True
                                else:
                                    all_names = literal
            self.modules[module] = _ModuleInfo(
                tree=tree,
                definitions=definitions,
                imports=imports,
                all_names=all_names,
                dynamic_all=dynamic_all,
            )

    def resolve(self, module: str, symbol: str) -> _Origin:
        return self._resolve(module, symbol, ())

    def _resolve(
        self,
        module: str,
        symbol: str,
        trail: tuple[tuple[str, str], ...],
    ) -> _Origin:
        key = (module, symbol)
        if key in self._resolve_cache:
            return self._resolve_cache[key]
        if key in trail:
            return _Origin("export_cycle", module, symbol)
        if module == "phospy.datasets" or module.startswith("phospy.datasets."):
            origin = _Origin("dataset_symbol", module, symbol)
            self._resolve_cache[key] = origin
            return origin
        info = self.modules.get(module)
        if info is None:
            return _Origin("external_symbol", module, symbol)
        imported = info.imports.get(symbol)
        if imported is not None:
            target_module, target_symbol, _level, _aliased = imported
            if target_symbol is None:
                if target_module == "phospy.datasets" or target_module.startswith(
                    "phospy.datasets."
                ):
                    return _Origin("dataset_module", target_module)
                return _Origin("module", target_module)
            origin = self._resolve(target_module, target_symbol, trail + (key,))
            self._resolve_cache[key] = origin
            return origin
        if symbol in info.definitions:
            origin = _Origin("definition", module, symbol)
            self._resolve_cache[key] = origin
            return origin
        return _Origin("unresolved_symbol", module, symbol)

    def exported_names(self, module: str) -> tuple[str, ...] | None:
        info = self.modules.get(module)
        if info is None or info.dynamic_all:
            return None
        if info.all_names is not None:
            return info.all_names
        return tuple(sorted((info.definitions | set(info.imports)) - {"__all__"}))

    def module_exposes_dataset(self, module: str) -> bool:
        if module == "phospy.datasets" or module.startswith("phospy.datasets."):
            return True
        names = self.exported_names(module)
        if names is None:
            return False
        return any(
            self.resolve(module, name).kind.startswith("dataset") for name in names
        )


def audit_signalome_sources(
    sources: Mapping[str, str],
    *,
    local_phospy_sources: Mapping[str, str] | None = None,
    require_complete_set: bool = False,
) -> DatasetBoundaryReport:
    """Audit in-memory Signalome sources without importing production modules."""

    index_sources = dict(local_phospy_sources or {})
    index_sources.update(sources)
    index = _Index(index_sources)
    observed: list[DatasetImport] = []
    violations: list[DatasetBoundaryViolation] = []
    for consumer in sorted(sources):
        visitor = _SignalomeVisitor(
            consumer=consumer,
            index=index,
            observed=observed,
            violations=violations,
        )
        visitor.visit(index.modules[consumer].tree)
    counts = Counter(observed)
    for triple, count in sorted(counts.items()):
        if count > 1:
            violations.append(
                DatasetBoundaryViolation(
                    "duplicate_import",
                    triple[0],
                    0,
                    f"approved Dataset import appears {count} times: {triple!r}",
                )
            )
    if require_complete_set:
        missing = APPROVED_DATASET_IMPORTS - set(counts)
        extra = set(counts) - APPROVED_DATASET_IMPORTS
        if missing or extra or any(count != 1 for count in counts.values()):
            violations.append(
                DatasetBoundaryViolation(
                    "import_set_mismatch",
                    "phospy.signalome",
                    0,
                    f"missing={sorted(missing)!r}; extra={sorted(extra)!r}; counts={counts!r}",
                )
            )
    return DatasetBoundaryReport(
        visited_modules=tuple(sorted(sources)),
        observed_imports=tuple(observed),
        violations=tuple(violations),
    )


def audit_signalome_tree(src_root: Path) -> DatasetBoundaryReport:
    """Audit every Python source below ``src/phospy/signalome``."""

    package_root = src_root / "phospy"
    signalome_root = package_root / "signalome"
    all_sources = {
        _module_name(path, src_root): path.read_text(encoding="utf-8")
        for path in sorted(package_root.rglob("*.py"))
    }
    signalome_sources = {
        _module_name(path, src_root): path.read_text(encoding="utf-8")
        for path in sorted(signalome_root.rglob("*.py"))
    }
    return audit_signalome_sources(
        signalome_sources,
        local_phospy_sources=all_sources,
        require_complete_set=True,
    )


class _SignalomeVisitor(ast.NodeVisitor):
    def __init__(
        self,
        *,
        consumer: str,
        index: _Index,
        observed: list[DatasetImport],
        violations: list[DatasetBoundaryViolation],
    ) -> None:
        self.consumer = consumer
        self.index = index
        self.observed = observed
        self.violations = violations
        self.bindings: dict[str, _Origin] = {
            "__import__": _Origin("callable", "builtins", "__import__"),
            "getattr": _Origin("callable", "builtins", "getattr"),
            "hasattr": _Origin("callable", "builtins", "hasattr"),
            "vars": _Origin("callable", "builtins", "vars"),
        }
        self.strings: dict[str, str] = {}

    def _violate(self, category: str, node: ast.AST, detail: str) -> None:
        self.violations.append(
            DatasetBoundaryViolation(
                category, self.consumer, getattr(node, "lineno", 0), detail
            )
        )

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            bound = alias.asname or alias.name.split(".")[0]
            if alias.name == "phospy":
                self.bindings[bound] = _PHOSPY_MODULE
            elif alias.name == "importlib":
                self.bindings[bound] = _Origin("module", "importlib")
            elif alias.name == "operator":
                self.bindings[bound] = _Origin("module", "operator")
            elif alias.name == "object":
                self.bindings[bound] = _Origin("module", "object")
            elif alias.name == "phospy.datasets" or alias.name.startswith(
                "phospy.datasets."
            ):
                self._violate("whole_module_import", node, alias.name)
                self.bindings[bound] = _DATASET_MODULE
            elif self.index.module_exposes_dataset(alias.name):
                self._violate("reexport_module_import", node, alias.name)
                self.bindings[bound] = _DATASET_MODULE
            else:
                self.bindings[bound] = _Origin("module", alias.name)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        module = _resolve_from_module(self.consumer, node.module, node.level)
        if any(alias.name == "*" for alias in node.names):
            if self.index.module_exposes_dataset(module):
                self._violate("wildcard_import", node, module)
            return
        for alias in node.names:
            bound = alias.asname or alias.name
            origin = self.index.resolve(module, alias.name)
            callable_name = f"{module}.{alias.name}"
            if callable_name in _ACQUISITION_CALLABLES:
                self.bindings[bound] = _Origin("callable", module, alias.name)
                continue
            if origin.kind == "export_cycle":
                self._violate("export_cycle", node, f"{module}.{alias.name}")
                continue
            if origin.kind == "dataset_symbol":
                triple = (self.consumer, module, alias.name)
                direct = origin.module == module and origin.symbol == alias.name
                if node.level != 0:
                    self._violate("relative_import", node, repr(triple))
                elif alias.asname is not None:
                    self._violate("aliased_import", node, repr(triple))
                elif not direct:
                    self._violate("reexport_import", node, repr(triple))
                elif triple not in APPROVED_DATASET_IMPORTS:
                    category = (
                        "wrong_consumer"
                        if any(
                            approved[1:] == triple[1:]
                            for approved in APPROVED_DATASET_IMPORTS
                        )
                        else "unapproved_dataset_symbol"
                    )
                    self._violate(category, node, repr(triple))
                else:
                    self.observed.append(triple)
                self.bindings[bound] = (
                    _DATASET_OBJECT
                    if alias.name in _DATASET_TYPES
                    else _Origin("dataset_type", origin.module, origin.symbol)
                )
                continue
            if module == "phospy.datasets" or module.startswith("phospy.datasets."):
                self._violate("whole_module_import", node, f"{module}.{alias.name}")
                self.bindings[bound] = _DATASET_MODULE
                continue
            if origin.kind == "dataset_module":
                self._violate("reexport_module_import", node, f"{module}.{alias.name}")
                self.bindings[bound] = _DATASET_MODULE
            else:
                self.bindings[bound] = _origin_to_taint(origin)

    def visit_Assign(self, node: ast.Assign) -> None:
        value_origin = self._eval(node.value)
        string_value = self._string(node.value)
        for target in node.targets:
            self._bind_target(target, value_origin, string_value)

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        annotation_origin = self._annotation_origin(node.annotation)
        value_origin = self._eval(node.value) if node.value is not None else _UNKNOWN
        origin = (
            _DATASET_OBJECT
            if annotation_origin.kind == "dataset_object"
            else value_origin
        )
        self._bind_target(node.target, origin, self._string(node.value))

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        for arg in (*node.args.posonlyargs, *node.args.args, *node.args.kwonlyargs):
            if arg.annotation is not None:
                self.bindings[arg.arg] = self._annotation_origin(arg.annotation)
            elif arg.arg in {"dataset", "dataset_view"}:
                self.bindings[arg.arg] = _DATASET_OBJECT
        for statement in node.body:
            self.visit(statement)

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        for statement in node.body:
            self.visit(statement)

    def visit_If(self, node: ast.If) -> None:
        self._eval(node.test)
        for statement in (*node.body, *node.orelse):
            self.visit(statement)

    def visit_For(self, node: ast.For) -> None:
        origin = self._eval(node.iter)
        self._bind_target(node.target, origin, None)
        for statement in (*node.body, *node.orelse):
            self.visit(statement)

    visit_AsyncFor = visit_For

    def visit_While(self, node: ast.While) -> None:
        self._eval(node.test)
        for statement in (*node.body, *node.orelse):
            self.visit(statement)

    def visit_With(self, node: ast.With) -> None:
        for item in node.items:
            origin = self._eval(item.context_expr)
            if item.optional_vars is not None:
                self._bind_target(item.optional_vars, origin, None)
        for statement in node.body:
            self.visit(statement)

    visit_AsyncWith = visit_With

    def visit_Expr(self, node: ast.Expr) -> None:
        self._eval(node.value)

    def visit_Return(self, node: ast.Return) -> None:
        if node.value is not None:
            self._eval(node.value)

    def visit_Raise(self, node: ast.Raise) -> None:
        if node.exc is not None:
            self._eval(node.exc)

    def _bind_target(
        self, target: ast.AST, origin: _Origin, string_value: str | None
    ) -> None:
        if isinstance(target, ast.Name):
            self.bindings[target.id] = origin
            if string_value is not None:
                self.strings[target.id] = string_value
            else:
                self.strings.pop(target.id, None)
        elif isinstance(target, (ast.Tuple, ast.List)):
            for element in target.elts:
                self._bind_target(element, origin, None)

    def _annotation_origin(self, node: ast.AST) -> _Origin:
        names = {part for part in ast.walk(node) if isinstance(part, ast.Name)}
        attrs = {
            part.attr for part in ast.walk(node) if isinstance(part, ast.Attribute)
        }
        if ({part.id for part in names} | attrs) & _DATASET_TYPES:
            return _DATASET_OBJECT
        return _UNKNOWN

    def _eval(self, node: ast.AST | None) -> _Origin:
        if node is None:
            return _UNKNOWN
        if isinstance(node, ast.Name):
            return self.bindings.get(node.id, _UNKNOWN)
        if isinstance(node, ast.Constant):
            return _NON_DATASET
        if isinstance(node, ast.Attribute):
            base = self._eval(node.value)
            if base.kind == "phospy_module" and node.attr == "datasets":
                self._violate("whole_module_import", node, "phospy.datasets")
                return _DATASET_MODULE
            if (
                base.kind == "module"
                and f"{base.module}.{node.attr}" in _ACQUISITION_CALLABLES
            ):
                module, symbol = f"{base.module}.{node.attr}".rsplit(".", 1)
                return _Origin("callable", module, symbol)
            if base.kind in {"dataset_module", "unknown_module"}:
                self._violate("module_indirection", node, node.attr)
                return _DATASET_OBJECT if node.attr in _DATASET_TYPES else base
            if base.kind == "dataset_object" and _is_protected(node.attr):
                self._violate("private_storage", node, node.attr)
            if base.kind == "dataset_object" and node.attr == "__dict__":
                return _DATASET_OBJECT
            if node.attr in {"dataset", "dataset_view"} and _contract_field_base(
                node.value
            ):
                return _DATASET_OBJECT
            return _UNKNOWN
        if isinstance(node, ast.Subscript):
            base = self._eval(node.value)
            member = self._string(node.slice)
            if base.kind in {"dataset_module", "unknown_module"}:
                self._violate("module_indirection", node, member or "<dynamic>")
                return _DATASET_OBJECT if member in _DATASET_TYPES else base
            if base.kind == "dataset_object" and member and _is_protected(member):
                self._violate("private_storage", node, member)
            return base
        if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
            origins = [self._eval(element) for element in node.elts]
            return _combine_origins(origins)
        if isinstance(node, ast.Dict):
            origins = [self._eval(value) for value in node.values]
            return _combine_origins(origins)
        if isinstance(node, ast.BinOp):
            self._eval(node.left)
            self._eval(node.right)
            return _NON_DATASET
        if isinstance(node, ast.Call):
            return self._eval_call(node)
        if isinstance(node, ast.IfExp):
            self._eval(node.test)
            return _combine_origins([self._eval(node.body), self._eval(node.orelse)])
        if isinstance(node, ast.UnaryOp):
            return self._eval(node.operand)
        for child in ast.iter_child_nodes(node):
            self._eval(child)
        return _UNKNOWN

    def _eval_call(self, node: ast.Call) -> _Origin:
        function = self._eval(node.func)
        callable_name = (
            f"{function.module}.{function.symbol}"
            if function.kind == "callable"
            else None
        )
        if callable_name in {"importlib.import_module", "builtins.__import__"}:
            target = self._string(node.args[0]) if node.args else None
            prefix = self._string_prefix(node.args[0]) if node.args else None
            if target == "phospy.datasets" or (
                target is not None and target.startswith("phospy.datasets.")
            ):
                self._violate("dynamic_import", node, target)
                return _DATASET_MODULE
            if prefix is not None and prefix.startswith("phospy.datasets"):
                self._violate("dynamic_import", node, prefix + "<dynamic>")
                return _UNKNOWN_MODULE
            return _NON_DATASET if target is not None else _UNKNOWN_MODULE
        if callable_name in {
            "builtins.getattr",
            "builtins.hasattr",
            "object.__getattribute__",
        }:
            base = self._eval(node.args[0]) if node.args else _UNKNOWN
            member = self._string(node.args[1]) if len(node.args) > 1 else None
            self._check_reflection(node, base, member)
            return _UNKNOWN
        if callable_name == "builtins.vars":
            base = self._eval(node.args[0]) if node.args else _UNKNOWN
            if base.kind in {"dataset_module", "unknown_module"}:
                self._violate("module_indirection", node, "vars")
                return base
            return _DATASET_OBJECT if base.kind == "dataset_object" else _UNKNOWN
        if callable_name == "operator.attrgetter":
            member = self._string(node.args[0]) if node.args else None
            return _Origin("attrgetter", symbol=member)
        if function.kind == "attrgetter":
            base = self._eval(node.args[0]) if node.args else _UNKNOWN
            self._check_reflection(node, base, function.symbol)
            return _UNKNOWN
        if isinstance(node.func, ast.Name):
            bound = self.bindings.get(node.func.id, _UNKNOWN)
            if (
                bound.kind in {"dataset_object", "dataset_type"}
                or node.func.id in _DATASET_TYPES
            ):
                for argument in (
                    *node.args,
                    *(keyword.value for keyword in node.keywords),
                ):
                    self._eval(argument)
                return _DATASET_OBJECT
        for argument in (*node.args, *(keyword.value for keyword in node.keywords)):
            self._eval(argument)
        return _UNKNOWN

    def _check_reflection(
        self, node: ast.AST, base: _Origin, member: str | None
    ) -> None:
        if base.kind in {"dataset_module", "unknown_module"}:
            self._violate("reflective_acquisition", node, member or "<dynamic>")
        elif base.kind == "dataset_object" and (
            member is None or _is_protected(member)
        ):
            self._violate("private_storage", node, member or "<dynamic>")

    def _string(self, node: ast.AST | None) -> str | None:
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value
        if isinstance(node, ast.Name):
            return self.strings.get(node.id)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            left = self._string(node.left)
            right = self._string(node.right)
            if left is not None and right is not None:
                return left + right
        return None

    def _string_prefix(self, node: ast.AST) -> str | None:
        exact = self._string(node)
        if exact is not None:
            return exact
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return self._string_prefix(node.left)
        if isinstance(node, ast.Name):
            return self.strings.get(node.id)
        return None


def _origin_to_taint(origin: _Origin) -> _Origin:
    if origin.kind == "dataset_symbol":
        return _DATASET_OBJECT
    if origin.kind == "module":
        return origin
    return _UNKNOWN


def _combine_origins(origins: list[_Origin]) -> _Origin:
    for kind in ("dataset_module", "unknown_module", "dataset_object"):
        if any(origin.kind == kind for origin in origins):
            return _Origin(kind)
    return _UNKNOWN


def _is_protected(name: str) -> bool:
    return name in _PROTECTED_MEMBERS or name.startswith("_borrow_")


def _contract_field_base(node: ast.AST) -> bool:
    if isinstance(node, ast.Name):
        return node.id in {"request", "result", "resolved", "kinase_result", "self"}
    if isinstance(node, ast.Attribute):
        return node.attr in {
            "request",
            "result",
            "kinase_result",
        } or _contract_field_base(node.value)
    return False


def _module_name(path: Path, src_root: Path) -> str:
    relative = path.relative_to(src_root).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _resolve_from_module(consumer: str, module: str | None, level: int) -> str:
    if level == 0:
        return module or ""
    package = (
        consumer if consumer.endswith(".__init__") else consumer.rpartition(".")[0]
    )
    if not package:
        package = consumer
    relative = "." * level + (module or "")
    try:
        return importlib.util.resolve_name(relative, package)
    except (ImportError, ValueError):
        return relative


def _literal_name_sequence(
    node: ast.AST | None, values: Mapping[str, tuple[str, ...]]
) -> tuple[str, ...] | None:
    if isinstance(node, (ast.List, ast.Tuple)):
        names: list[str] = []
        for element in node.elts:
            if not isinstance(element, ast.Constant) or not isinstance(
                element.value, str
            ):
                return None
            names.append(element.value)
        return tuple(names)
    if isinstance(node, ast.Name):
        return values.get(node.id)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _literal_name_sequence(node.left, values)
        right = _literal_name_sequence(node.right, values)
        if left is not None and right is not None:
            return left + right
    return None


__all__ = [
    "APPROVED_DATASET_IMPORTS",
    "DatasetBoundaryReport",
    "DatasetBoundaryViolation",
    "audit_signalome_sources",
    "audit_signalome_tree",
]
