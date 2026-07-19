from __future__ import annotations

import argparse
import hashlib
import json
import platform
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib

SOURCE_CHECK_REPORT_SCHEMA = "phospy.source-check/v1"


class SourceCheckReportError(RuntimeError):
    """Raised when a source-check report cannot be produced from evidence."""


def write_source_check_report(
    *,
    suite_id: str,
    report_classes: Sequence[str],
    junit_xml: Path,
    source_identity_path: Path,
    output_path: Path,
    repository_root: Path,
    command: str | None = None,
    python_version: str | None = None,
    package_name: str | None = None,
    package_version: str | None = None,
) -> Path:
    payload = build_source_check_report(
        suite_id=suite_id,
        report_classes=report_classes,
        junit_xml=junit_xml,
        source_identity_path=source_identity_path,
        repository_root=repository_root,
        command=command,
        python_version=python_version,
        package_name=package_name,
        package_version=package_version,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def build_source_check_report(
    *,
    suite_id: str,
    report_classes: Sequence[str],
    junit_xml: Path,
    source_identity_path: Path,
    repository_root: Path,
    command: str | None = None,
    python_version: str | None = None,
    package_name: str | None = None,
    package_version: str | None = None,
) -> dict[str, Any]:
    normalized_suite_id = _required_text(suite_id, field_name="suite id")
    normalized_classes = _unique_texts(report_classes, field_name="report class")
    junit_summary = _junit_summary(junit_xml)
    status = (
        "success"
        if junit_summary["failures"] == 0 and junit_summary["errors"] == 0
        else "failure"
    )
    package = _package_metadata(
        repository_root,
        package_name=package_name,
        package_version=package_version,
    )
    return {
        "schema": SOURCE_CHECK_REPORT_SCHEMA,
        "status": status,
        "source_identity_digest": _prefixed_file_sha256(source_identity_path),
        "package": package,
        "suite": {
            "id": normalized_suite_id,
            "report_classes": normalized_classes,
        },
        "python": python_version or platform.python_version(),
        "command": command,
        "evidence": {
            "junit_xml": junit_xml.name,
            "junit_xml_sha256": _prefixed_file_sha256(junit_xml),
            "summary": junit_summary,
        },
    }


def _junit_summary(path: Path) -> dict[str, int]:
    _require(path.is_file(), f"JUnit XML report is missing: {path.as_posix()}")
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise SourceCheckReportError(
            f"JUnit XML report is malformed: {path.as_posix()}: {exc}"
        ) from exc

    elements = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    _require(elements, f"JUnit XML report has no testsuite elements: {path.as_posix()}")
    totals = {"tests": 0, "failures": 0, "errors": 0, "skipped": 0}
    for element in elements:
        totals["tests"] += _int_attr(element, "tests")
        totals["failures"] += _int_attr(element, "failures")
        totals["errors"] += _int_attr(element, "errors")
        totals["skipped"] += _int_attr(element, "skipped")
    return totals


def _package_metadata(
    repository_root: Path,
    *,
    package_name: str | None,
    package_version: str | None,
) -> dict[str, str]:
    pyproject_package = _pyproject_package(repository_root)
    name = package_name or pyproject_package.get("name")
    version = package_version or pyproject_package.get("version")
    _require(
        isinstance(name, str) and bool(name.strip()),
        "package name is required",
    )
    _require(
        isinstance(version, str) and bool(version.strip()),
        "package version is required",
    )
    return {"name": name.strip(), "version": version.strip()}


def _pyproject_package(repository_root: Path) -> dict[str, str | None]:
    pyproject = repository_root / "pyproject.toml"
    if not pyproject.is_file():
        return {"name": None, "version": None}
    with pyproject.open("rb") as handle:
        payload = tomllib.load(handle)
    project = payload.get("project")
    if not isinstance(project, Mapping):
        return {"name": None, "version": None}
    name = project.get("name")
    version = project.get("version")
    return {
        "name": name if isinstance(name, str) else None,
        "version": version if isinstance(version, str) else None,
    }


def _int_attr(element: ET.Element, name: str) -> int:
    raw = element.attrib.get(name, "0")
    try:
        return int(raw)
    except ValueError as exc:
        raise SourceCheckReportError(
            f"JUnit XML testsuite attribute {name!r} is not an integer: {raw!r}"
        ) from exc


def _prefixed_file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def _unique_texts(values: Sequence[str], *, field_name: str) -> list[str]:
    normalized: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = _required_text(value, field_name=field_name)
        if text in seen:
            raise SourceCheckReportError(f"duplicate {field_name}: {text!r}")
        normalized.append(text)
        seen.add(text)
    _require(normalized != [], f"at least one {field_name} is required")
    return normalized


def _required_text(value: object, *, field_name: str) -> str:
    _require(
        isinstance(value, str) and bool(value.strip()),
        f"{field_name} must be non-empty text",
    )
    return value.strip()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceCheckReportError(message)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    path = write_source_check_report(
        suite_id=args.suite_id,
        report_classes=tuple(args.report_class),
        junit_xml=Path(args.junit_xml),
        source_identity_path=Path(args.source_identity),
        output_path=Path(args.output),
        repository_root=Path(args.repository_root),
        command=args.command,
        python_version=args.python_version,
        package_name=args.package_name,
        package_version=args.package_version,
    )
    print(path.as_posix())
    return 0


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a source-check JSON report from a completed JUnit report."
    )
    parser.add_argument("--suite-id", required=True, help="Policy source-suite ID.")
    parser.add_argument(
        "--report-class",
        action="append",
        required=True,
        help="Report class represented by this source-suite report.",
    )
    parser.add_argument(
        "--junit-xml", required=True, help="Completed JUnit XML report."
    )
    parser.add_argument(
        "--source-identity",
        required=True,
        help="Source identity record used for this source-check run.",
    )
    parser.add_argument("--output", required=True, help="Output JSON report path.")
    parser.add_argument(
        "--repository-root",
        default=".",
        help="Source tree root used only for package metadata fallback.",
    )
    parser.add_argument(
        "--command", default=None, help="Command that produced the report."
    )
    parser.add_argument(
        "--python-version",
        default=None,
        help="Python version for the source-check report.",
    )
    parser.add_argument("--package-name", default=None, help="Package name override.")
    parser.add_argument(
        "--package-version",
        default=None,
        help="Package version override.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
