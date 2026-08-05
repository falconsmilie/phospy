from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import phospy.frames.validation as frame_validation
import phospy.science.sites.identity_columns as site_identity_column_compat
import phospy.science.sites.identity_rules.contracts as site_identity_contract_records
import phospy.science.sites.identity_rules.dataset_identity as site_dataset_identity
import phospy.science.sites.identity_rules.reference_context as site_reference_context
import phospy.science.sites.identity_rules.result_identity as site_result_identity
import phospy.science.sites.metadata_validation as site_metadata_validation
import phospy.science.sites.sequence_context as site_sequence_context
import phospy.validation.common.dataframes as validation_dataframe_compat
import phospy.validation.datasets.protein_scoped_site_identity as site_key_compat
import phospy.validation.datasets.site_metadata as site_metadata_compat
import phospy.validation.identity_contracts as identity_contract_compat

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
DOCS_ROOT = PROJECT_ROOT / "docs"
OWNERSHIP_MAP_PATH = DOCS_ROOT / "validation-ownership.md"

_CODE_SPAN_PATTERN = re.compile(r"`([^`]+)`")
_PHOSPY_REFERENCE_PATTERN = re.compile(r"\bphospy(?:\.[A-Za-z_]\w*)+")


@dataclass(frozen=True, slots=True)
class OwnershipRow:
    invariant: str
    owner: str
    enforcement_point: str
    should_not_live_in: str
    tests: str


@dataclass(frozen=True, slots=True)
class ConcreteSymbolOwner:
    symbol_name: str
    owner_path: Path
    search_roots: tuple[Path, ...]


@dataclass(frozen=True, slots=True)
class ReexportRoute:
    reexport_module: object
    owner_module: object
    symbols: tuple[str, ...]


FRAME_SEARCH_ROOTS = (
    PACKAGE_ROOT / "frames",
    PACKAGE_ROOT / "validation" / "common",
    PACKAGE_ROOT / "science" / "sites",
    PACKAGE_ROOT / "workflows",
)
SITE_IDENTITY_SEARCH_ROOTS = (
    PACKAGE_ROOT / "science" / "sites",
    PACKAGE_ROOT / "validation",
    PACKAGE_ROOT / "workflows",
)

CONCRETE_VALIDATION_SYMBOL_OWNERS = (
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=PACKAGE_ROOT / "frames" / "validation.py",
            search_roots=FRAME_SEARCH_ROOTS,
        )
        for symbol_name in (
            "require_dataframe",
            "require_non_empty_dataframe",
            "require_aligned_dataframe_shape",
            "require_numeric_dataframe",
            "require_finite_numeric_dataframe",
            "require_unique_index",
            "require_unique_columns",
            "require_exact_index_match",
            "require_no_duplicate_labels",
            "require_unique_row_pairs",
        )
    ),
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=PACKAGE_ROOT / "science" / "sites" / "identity.py",
            search_roots=SITE_IDENTITY_SEARCH_ROOTS,
        )
        for symbol_name in ("build_phosphosite_identity",)
    ),
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=PACKAGE_ROOT / "science" / "sites" / "validation.py",
            search_roots=SITE_IDENTITY_SEARCH_ROOTS,
        )
        for symbol_name in (
            "require_site_key_index",
            "require_site_key_series",
            "require_site_identity_coherence",
        )
    ),
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=PACKAGE_ROOT / "science" / "sites" / "metadata_validation.py",
            search_roots=SITE_IDENTITY_SEARCH_ROOTS,
        )
        for symbol_name in (
            "validate_site_identity_metadata",
            "validate_site_sequence_column",
            "enforce_site_identity_rows",
        )
    ),
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=PACKAGE_ROOT / "science" / "sites" / "sequence_context.py",
            search_roots=SITE_IDENTITY_SEARCH_ROOTS,
        )
        for symbol_name in (
            "SequenceContextContract",
            "enforce_site_sequence_context_contract",
            "enforce_centred_site_sequence_context",
        )
    ),
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=(
                PACKAGE_ROOT / "science" / "sites" / "identity_rules" / "contracts.py"
            ),
            search_roots=SITE_IDENTITY_SEARCH_ROOTS,
        )
        for symbol_name in (
            "PhosphositeIdentityContract",
            "ReferenceContextCompatibilityWarning",
        )
    ),
    ConcreteSymbolOwner(
        symbol_name="validate_reference_context_compatibility",
        owner_path=(
            PACKAGE_ROOT
            / "science"
            / "sites"
            / "identity_rules"
            / "reference_context.py"
        ),
        search_roots=SITE_IDENTITY_SEARCH_ROOTS,
    ),
    *(
        ConcreteSymbolOwner(
            symbol_name=symbol_name,
            owner_path=(
                PACKAGE_ROOT
                / "science"
                / "sites"
                / "identity_rules"
                / "dataset_identity.py"
            ),
            search_roots=SITE_IDENTITY_SEARCH_ROOTS,
        )
        for symbol_name in (
            "enforce_phosphosite_identity_contract",
            "enforce_analysis_ready_site_key_index",
            "enforce_site_key_column",
            "enforce_display_id_column",
            "enforce_unique_site_key_identity",
            "enforce_site_key_column_matches_index",
            "enforce_site_key_matches_metadata",
        )
    ),
    ConcreteSymbolOwner(
        symbol_name="enforce_result_identity_metadata_coherence",
        owner_path=(
            PACKAGE_ROOT / "science" / "sites" / "identity_rules" / "result_identity.py"
        ),
        search_roots=SITE_IDENTITY_SEARCH_ROOTS,
    ),
)

IDENTITY_PRESERVING_REEXPORT_ROUTES = (
    ReexportRoute(
        reexport_module=validation_dataframe_compat,
        owner_module=frame_validation,
        symbols=(
            "require_dataframe",
            "require_non_empty_dataframe",
            "require_aligned_dataframe_shape",
            "require_numeric_dataframe",
            "require_finite_numeric_dataframe",
            "require_unique_index",
            "require_unique_columns",
            "require_exact_index_match",
            "require_no_duplicate_labels",
            "require_unique_row_pairs",
        ),
    ),
    ReexportRoute(
        reexport_module=identity_contract_compat,
        owner_module=site_identity_contract_records,
        symbols=(
            "PhosphositeIdentityContract",
            "ReferenceContextCompatibilityWarning",
            "SequenceContextContract",
        ),
    ),
    ReexportRoute(
        reexport_module=identity_contract_compat,
        owner_module=site_reference_context,
        symbols=("validate_reference_context_compatibility",),
    ),
    ReexportRoute(
        reexport_module=identity_contract_compat,
        owner_module=site_dataset_identity,
        symbols=(
            "enforce_phosphosite_identity_contract",
            "enforce_analysis_ready_site_key_index",
            "enforce_site_key_column",
            "enforce_display_id_column",
            "enforce_unique_site_key_identity",
            "enforce_site_key_column_matches_index",
            "enforce_site_key_matches_metadata",
        ),
    ),
    ReexportRoute(
        reexport_module=identity_contract_compat,
        owner_module=site_result_identity,
        symbols=("enforce_result_identity_metadata_coherence",),
    ),
    ReexportRoute(
        reexport_module=site_key_compat,
        owner_module=site_dataset_identity,
        symbols=(
            "enforce_analysis_ready_site_key_index",
            "enforce_site_key_column",
            "enforce_display_id_column",
            "enforce_unique_site_key_identity",
            "enforce_site_key_column_matches_index",
            "enforce_site_key_matches_metadata",
        ),
    ),
    ReexportRoute(
        reexport_module=site_metadata_compat,
        owner_module=site_metadata_validation,
        symbols=(
            "validate_site_identity_metadata",
            "validate_site_sequence_column",
            "enforce_site_identity_rows",
        ),
    ),
    ReexportRoute(
        reexport_module=site_metadata_compat,
        owner_module=site_sequence_context,
        symbols=(
            "SequenceContextContract",
            "enforce_site_sequence_context_contract",
            "enforce_centred_site_sequence_context",
        ),
    ),
    ReexportRoute(
        reexport_module=site_identity_column_compat,
        owner_module=site_dataset_identity,
        symbols=("enforce_display_id_column", "enforce_site_key_column"),
    ),
)


def test_duplicate_symbol_fixture_rejects_duplicate_concrete_owner(
    tmp_path: Path,
) -> None:
    owner_root = tmp_path / "owner"
    duplicate_root = tmp_path / "duplicate"
    owner_root.mkdir()
    duplicate_root.mkdir()
    owner_path = owner_root / "identity_contracts.py"
    duplicate_path = duplicate_root / "copied_contracts.py"
    owner_path.write_text(
        "class PhosphositeIdentityContract:\n    pass\n",
        encoding="utf-8",
    )
    duplicate_path.write_text(
        "class PhosphositeIdentityContract:\n    pass\n",
        encoding="utf-8",
    )

    errors = _concrete_symbol_owner_errors(
        (
            ConcreteSymbolOwner(
                symbol_name="PhosphositeIdentityContract",
                owner_path=owner_path,
                search_roots=(owner_root, duplicate_root),
            ),
        )
    )

    assert errors == [
        (
            "PhosphositeIdentityContract: expected one concrete definition in "
            f"{owner_path}; observed {owner_path}, {duplicate_path}"
        )
    ]


def test_concrete_validation_symbols_have_single_owner_definition() -> None:
    assert _concrete_symbol_owner_errors(CONCRETE_VALIDATION_SYMBOL_OWNERS) == []


def test_identity_preserving_reexport_fixture_detects_copied_object() -> None:
    owner_object = object()
    copied_object = object()
    owner_module = SimpleNamespace(Symbol=owner_object)
    reexport_module = SimpleNamespace(Symbol=copied_object)

    assert _reexport_identity_mismatches(
        ReexportRoute(
            reexport_module=reexport_module,
            owner_module=owner_module,
            symbols=("Symbol",),
        )
    ) == ["Symbol"]


def test_deliberate_validation_compatibility_reexports_preserve_identity() -> None:
    mismatches = [
        mismatch
        for route in IDENTITY_PRESERVING_REEXPORT_ROUTES
        for mismatch in _reexport_identity_mismatches(route)
    ]

    assert mismatches == []


def test_stale_ownership_map_module_path_fixture_rejects_missing_module() -> None:
    rows = _ownership_rows_from_text(
        """
| Invariant | Owner | Enforcement point | Should not live in | Tests |
| --- | --- | --- | --- | --- |
| Fake invariant | `phospy.science.sites.missing_owner.Fake` | `phospy.frames.validation.require_dataframe` | Workflow validators | `tests/architecture/test_validation_ownership_governance.py` |
"""
    )

    assert _stale_module_reference_errors(rows) == [
        "Fake invariant owner: phospy.science.sites.missing_owner.Fake"
    ]


def test_validation_ownership_map_module_paths_and_test_files_are_current() -> None:
    rows = _ownership_rows()

    assert rows
    assert _stale_module_reference_errors(rows) == []
    assert _stale_test_reference_errors(rows) == []


def test_validation_ownership_map_has_unique_reviewable_rows() -> None:
    rows = _ownership_rows()
    invariants = [row.invariant for row in rows]
    duplicated = sorted(
        invariant for invariant in set(invariants) if invariants.count(invariant) > 1
    )
    rows_without_tests = [row.invariant for row in rows if not _code_spans(row.tests)]

    assert duplicated == []
    assert rows_without_tests == []


def test_validation_ownership_map_assigns_frame_and_science_owners() -> None:
    rows = {row.invariant: row for row in _ownership_rows()}

    assert _first_code_span(rows["DataFrame shape"].owner) == "phospy.frames.validation"
    assert (
        _first_code_span(rows["Unique site/sample labels"].owner)
        == "phospy.frames.validation"
    )
    assert (
        _first_code_span(
            rows[
                "Strict phosphosite identifier format and site-ID/metadata coherence"
            ].owner
        )
        == "phospy.science.sites.validation"
    )
    assert (
        _first_code_span(rows["Reusable phosphosite identity contracts"].owner)
        == "phospy.science.sites.identity_rules.contracts"
    )
    assert (
        _first_code_span(rows["Analysis-ready `site_key` row identity"].owner)
        == "phospy.science.sites.identity_rules.dataset_identity"
    )
    assert (
        _first_code_span(rows["Analysis-ready `site_sequence` evidence"].owner)
        == "phospy.science.sites.metadata_validation.validate_site_sequence_column"
    )
    assert (
        _first_code_span(rows["Workflow-specific centered sequence context"].owner)
        == "phospy.science.sites.sequence_context"
    )


def test_ownership_map_records_dataset_construction_boundary_for_dataset_invariants() -> (
    None
):
    rows = {row.invariant: row for row in _ownership_rows()}
    dataset_boundary_rows = (
        "Finite numeric intensity values",
        "Unique site/sample labels",
        "Analysis-ready `site_key` row identity",
        "Required site metadata",
        "Analysis-ready `site_sequence` evidence",
        "Processing state coherence",
    )

    for invariant in dataset_boundary_rows:
        assert (
            "_AnalysisReadyDatasetConstructionService._validate_analysis_ready_tables"
            in rows[invariant].enforcement_point
        )
        assert (
            "private construction service boundary"
            in rows[invariant].enforcement_point.lower()
        )


def test_ownership_map_records_workflow_validators_as_composers() -> None:
    rows = {row.invariant: row for row in _ownership_rows()}
    shared_scientific_rows = (
        "Strict phosphosite identifier format and site-ID/metadata coherence",
        "Reusable phosphosite identity contracts",
        "Analysis-ready `site_key` row identity",
        "Analysis-ready `site_sequence` evidence",
        "Workflow-specific centered sequence context",
        "Localisation eligibility",
    )

    for invariant in shared_scientific_rows:
        row = rows[invariant]
        primary_owner = _first_code_span(row.owner)
        assert not primary_owner.startswith("phospy.workflows")
        assert "compose" in f"{row.owner} {row.enforcement_point}".lower()


def _concrete_symbol_owner_errors(
    specs: tuple[ConcreteSymbolOwner, ...],
) -> list[str]:
    errors: list[str] = []
    for spec in specs:
        observed = _definition_paths(
            symbol_name=spec.symbol_name,
            roots=spec.search_roots,
        )
        expected = (spec.owner_path.resolve(),)
        if observed != expected:
            observed_text = ", ".join(str(path) for path in observed) or "(none)"
            errors.append(
                f"{spec.symbol_name}: expected one concrete definition in "
                f"{expected[0]}; observed {observed_text}"
            )
    return errors


def _definition_paths(*, symbol_name: str, roots: tuple[Path, ...]) -> tuple[Path, ...]:
    paths: list[Path] = []
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in tree.body:
                if (
                    isinstance(
                        node,
                        ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef,
                    )
                    and node.name == symbol_name
                ):
                    paths.append(path.resolve())
    return tuple(paths)


def _reexport_identity_mismatches(route: ReexportRoute) -> list[str]:
    mismatches: list[str] = []
    for symbol_name in route.symbols:
        if getattr(route.reexport_module, symbol_name, None) is not getattr(
            route.owner_module,
            symbol_name,
            None,
        ):
            mismatches.append(symbol_name)
    return mismatches


def _ownership_rows() -> tuple[OwnershipRow, ...]:
    return _ownership_rows_from_text(OWNERSHIP_MAP_PATH.read_text(encoding="utf-8"))


def _ownership_rows_from_text(source: str) -> tuple[OwnershipRow, ...]:
    rows: list[OwnershipRow] = []
    in_table = False
    for line in source.splitlines():
        if line.startswith("| Invariant | Owner | Enforcement point |"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != 5:
            continue
        if set(cells) == {"---"}:
            continue
        rows.append(
            OwnershipRow(
                invariant=cells[0],
                owner=cells[1],
                enforcement_point=cells[2],
                should_not_live_in=cells[3],
                tests=cells[4],
            )
        )
    return tuple(rows)


def _stale_module_reference_errors(rows: tuple[OwnershipRow, ...]) -> list[str]:
    errors: list[str] = []
    for row in rows:
        cells = (
            ("owner", row.owner),
            ("enforcement point", row.enforcement_point),
            ("should not live in", row.should_not_live_in),
        )
        for cell_name, cell in cells:
            for reference in _phospy_references(cell):
                if _resolve_module_reference(reference) is None:
                    errors.append(f"{row.invariant} {cell_name}: {reference}")
    return sorted(set(errors))


def _stale_test_reference_errors(rows: tuple[OwnershipRow, ...]) -> list[str]:
    errors: list[str] = []
    for row in rows:
        for test_reference in _code_spans(row.tests):
            test_path = test_reference.split("::", maxsplit=1)[0]
            if not test_path.startswith("tests/"):
                continue
            if not (PROJECT_ROOT / test_path).exists():
                errors.append(f"{row.invariant} tests: {test_reference}")
    return sorted(set(errors))


def _phospy_references(cell: str) -> tuple[str, ...]:
    references: list[str] = []
    for code_span in _code_spans(cell):
        references.extend(
            match.group(0) for match in _PHOSPY_REFERENCE_PATTERN.finditer(code_span)
        )
    return tuple(references)


def _resolve_module_reference(reference: str) -> Path | None:
    parts = reference.split(".")
    if not parts or parts[0] != "phospy":
        return None
    candidate_lengths = (len(parts), len(parts) - 1)
    for end in candidate_lengths:
        if end <= 1:
            continue
        relative_parts = parts[1:end]
        module_file = PACKAGE_ROOT.joinpath(*relative_parts).with_suffix(".py")
        package_init = PACKAGE_ROOT.joinpath(*relative_parts, "__init__.py")
        if module_file.exists():
            return module_file
        if package_init.exists():
            return package_init
    return None


def _code_spans(cell: str) -> tuple[str, ...]:
    return tuple(_CODE_SPAN_PATTERN.findall(cell))


def _first_code_span(cell: str) -> str:
    spans = _code_spans(cell)
    assert spans
    return spans[0]
