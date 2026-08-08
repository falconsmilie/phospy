from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from fnmatch import fnmatchcase
from pathlib import Path
from typing import Any

CANONICAL_TEXT_BYTE_POLICY = "utf-8 LF with final newline"
TEXT_FIXTURE_SUFFIXES = frozenset(
    (
        ".csv",
        ".json",
        ".md",
        ".txt",
        ".tsv",
        ".gmt",
        ".R",
    )
)


@dataclass(frozen=True)
class FixtureByteReference:
    path: Path
    expected_sha256: str | None
    source_path: Path


@dataclass(frozen=True)
class _GitAttributesRule:
    pattern: str
    attributes: tuple[str, ...]


def is_manifest_text_fixture(path: Path) -> bool:
    return path.suffix in TEXT_FIXTURE_SUFFIXES


def display_path(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def assert_canonical_text_bytes(
    path: Path | str,
    data: bytes,
    *,
    byte_policy: str = CANONICAL_TEXT_BYTE_POLICY,
) -> None:
    display = path.as_posix() if isinstance(path, Path) else path
    try:
        data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise AssertionError(
            f"{display}: invalid UTF-8 fixture bytes at byte {exc.start}; "
            f"expected {byte_policy}"
        ) from exc

    if b"\r\n" in data:
        raise AssertionError(
            f"{display}: CRLF line endings violate fixture byte policy; "
            f"expected {byte_policy}"
        )
    if b"\r" in data:
        raise AssertionError(
            f"{display}: lone CR byte violates fixture byte policy; "
            f"expected {byte_policy}"
        )
    if not data.endswith(b"\n"):
        raise AssertionError(
            f"{display}: missing final newline violates fixture byte policy; "
            f"expected {byte_policy}"
        )


def assert_text_fixture_matches_sha256(
    path: Path,
    *,
    expected_sha256: str | None,
    repo_root: Path,
    byte_policy: str = CANONICAL_TEXT_BYTE_POLICY,
) -> None:
    display = display_path(path, repo_root)
    data = path.read_bytes()
    assert_canonical_text_bytes(display, data, byte_policy=byte_policy)
    if expected_sha256 is None:
        return
    actual_sha256 = hashlib.sha256(data).hexdigest()
    if actual_sha256 != expected_sha256:
        raise AssertionError(
            f"{display}: digest mismatch for raw fixture bytes; "
            f"expected sha256={expected_sha256} actual sha256={actual_sha256}"
        )


def iter_manifest_governed_text_fixture_references(
    repo_root: Path,
) -> tuple[FixtureByteReference, ...]:
    references: dict[str, FixtureByteReference] = {}
    for reference in (
        *_iter_manifest_file_references(repo_root),
        *_iter_fixture_index_references(repo_root),
    ):
        if not is_manifest_text_fixture(reference.path):
            continue
        key = display_path(reference.path, repo_root)
        existing = references.get(key)
        if existing is None:
            references[key] = reference
            continue
        if (
            existing.expected_sha256 is not None
            and reference.expected_sha256 is not None
            and existing.expected_sha256 != reference.expected_sha256
        ):
            raise AssertionError(
                f"{key}: contradictory fixture digests referenced by "
                f"{display_path(existing.source_path, repo_root)} and "
                f"{display_path(reference.source_path, repo_root)}"
            )
        if existing.expected_sha256 is None and reference.expected_sha256 is not None:
            references[key] = reference

    return tuple(
        sorted(
            references.values(),
            key=lambda reference: display_path(reference.path, repo_root),
        )
    )


def iter_importer_fixture_index_references(
    repo_root: Path,
) -> tuple[FixtureByteReference, ...]:
    index_path = (
        repo_root
        / "tests"
        / "fixtures"
        / "release_validation_regression"
        / "importer_edge_cases"
        / "fixture_index.json"
    )
    return tuple(
        reference
        for reference in _iter_fixture_index_references(repo_root, index_path)
        if is_manifest_text_fixture(reference.path)
    )


def assert_lf_gitattributes_coverage(
    repo_root: Path,
    fixture_paths: tuple[Path, ...],
) -> None:
    attributes_path = repo_root / ".gitattributes"
    rules = _read_gitattributes_rules(attributes_path)
    missing: list[str] = []
    for path in fixture_paths:
        if not is_manifest_text_fixture(path):
            continue
        attributes = _attributes_for_path(
            display_path(path, repo_root),
            rules,
        )
        if attributes.get("text") != "set" or attributes.get("eol") != "lf":
            missing.append(display_path(path, repo_root))

    if missing:
        raise AssertionError(
            "manifest-governed text fixtures missing LF .gitattributes coverage: "
            + ", ".join(sorted(missing))
        )


def _iter_manifest_file_references(
    repo_root: Path,
) -> tuple[FixtureByteReference, ...]:
    fixture_root = repo_root / "tests" / "fixtures"
    references: list[FixtureByteReference] = []
    for manifest_path in sorted(fixture_root.rglob("MANIFEST.json")):
        references.append(
            FixtureByteReference(
                path=manifest_path,
                expected_sha256=None,
                source_path=manifest_path,
            )
        )
        manifest = _read_json(manifest_path)
        for item in manifest.get("files", ()):
            if not isinstance(item, dict):
                continue
            relative_path = item.get("relative_path")
            if not isinstance(relative_path, str):
                continue
            expected_sha256 = item.get("sha256")
            references.append(
                FixtureByteReference(
                    path=manifest_path.parent / Path(relative_path),
                    expected_sha256=(
                        expected_sha256 if isinstance(expected_sha256, str) else None
                    ),
                    source_path=manifest_path,
                )
            )
    return tuple(references)


def _iter_fixture_index_references(
    repo_root: Path,
    index_path: Path | None = None,
) -> tuple[FixtureByteReference, ...]:
    fixture_root = repo_root / "tests" / "fixtures"
    index_paths = (
        (index_path,)
        if index_path is not None
        else fixture_root.rglob("fixture_index.json")
    )
    references: list[FixtureByteReference] = []
    for resolved_index_path in sorted(index_paths):
        index = _read_json(resolved_index_path)
        for item in index.get("referenced_fixture_files", ()):
            if not isinstance(item, dict):
                continue
            relative_path = item.get("relative_path")
            if not isinstance(relative_path, str):
                continue
            expected_sha256 = item.get("sha256")
            references.append(
                FixtureByteReference(
                    path=repo_root / Path(relative_path),
                    expected_sha256=(
                        expected_sha256 if isinstance(expected_sha256, str) else None
                    ),
                    source_path=resolved_index_path,
                )
            )
    return tuple(references)


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(payload, dict), f"{path.as_posix()} must contain a JSON object"
    return payload


def _read_gitattributes_rules(path: Path) -> tuple[_GitAttributesRule, ...]:
    rules: list[_GitAttributesRule] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        pattern, *attributes = stripped.split()
        rules.append(_GitAttributesRule(pattern=pattern, attributes=tuple(attributes)))
    return tuple(rules)


def _attributes_for_path(
    relative_path: str,
    rules: tuple[_GitAttributesRule, ...],
) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for rule in rules:
        if not _gitattributes_pattern_matches(rule.pattern, relative_path):
            continue
        for attribute in rule.attributes:
            if attribute == "text":
                attributes["text"] = "set"
            elif attribute == "-text":
                attributes["text"] = "unset"
            elif attribute.startswith("eol="):
                attributes["eol"] = attribute.split("=", maxsplit=1)[1]
    return attributes


def _gitattributes_pattern_matches(pattern: str, relative_path: str) -> bool:
    normalized_pattern = pattern.lstrip("/")
    if normalized_pattern.endswith("/**"):
        prefix = normalized_pattern[:-3].rstrip("/")
        return relative_path == prefix or relative_path.startswith(f"{prefix}/")
    if "/**/*" in normalized_pattern:
        prefix, suffix = normalized_pattern.split("/**/*", maxsplit=1)
        return relative_path.startswith(f"{prefix.rstrip('/')}/") and (
            relative_path.endswith(suffix)
        )
    return fnmatchcase(relative_path, normalized_pattern)
