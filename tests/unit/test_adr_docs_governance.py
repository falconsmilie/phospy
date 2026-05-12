"""ADR documentation governance checks."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ADR_ROOT = ROOT / "docs" / "adr"
ADR_INDEX = ADR_ROOT / "index.md"

_ADR_ID_PATTERN = re.compile(r"^- \*\*ADR ID:\*\*\s*(ADR-\d{4})\s*$", re.MULTILINE)
_ADR_STATUS_PATTERN = re.compile(r"^- \*\*Status:\*\*\s*(.+?)\s*$", re.MULTILINE)
_INDEX_ROW_PATTERN = re.compile(
    r"^\|\s*(ADR-\d{4})\s*\|.*\|\s*([A-Za-z]+)\s*\|.*\|\s*\[([^\]]+)\]\(([^)]+)\)\s*\|$",
    re.MULTILINE,
)
_RELATIVE_LINK_PATTERN = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
_SOURCE_REF_PATTERN = re.compile(r"src/phospy/[A-Za-z0-9_./-]+")

_ALLOWED_STATUSES = {"Accepted", "Superseded", "Amended", "Deprecated", "Draft"}


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _adr_files() -> tuple[Path, ...]:
    return tuple(sorted(ADR_ROOT.glob("adr_*.md")))


def _parse_adr_control(path: Path) -> tuple[str, str]:
    text = _read_text(path)
    id_match = _ADR_ID_PATTERN.search(text)
    status_match = _ADR_STATUS_PATTERN.search(text)
    assert id_match is not None, f"missing ADR ID control field: {path.as_posix()}"
    assert status_match is not None, (
        f"missing ADR Status control field: {path.as_posix()}"
    )
    return id_match.group(1), status_match.group(1).strip()


def _iter_relative_links(path: Path) -> tuple[str, ...]:
    text = _read_text(path)
    links: list[str] = []
    for match in _RELATIVE_LINK_PATTERN.finditer(text):
        target = match.group(1).strip()
        lowered = target.lower()
        if lowered.startswith(("http://", "https://", "mailto:")):
            continue
        if target.startswith("#"):
            continue
        links.append(target)
    return tuple(links)


def _resolve_markdown_target(base: Path, target: str) -> Path:
    clean_target = target.split("#", 1)[0].strip()
    return (base / clean_target).resolve()


def test_adr_index_status_and_link_consistency() -> None:
    index_text = _read_text(ADR_INDEX)
    rows = list(_INDEX_ROW_PATTERN.finditer(index_text))
    assert rows, "ADR index table rows were not parsed"

    adr_by_filename = {path.name: path for path in _adr_files()}
    indexed_filenames: set[str] = set()

    for row in rows:
        adr_id = row.group(1)
        status = row.group(2)
        link_label = row.group(3)
        link_target = row.group(4)
        link_name = Path(link_target).name

        assert status in _ALLOWED_STATUSES, (
            f"index uses unsupported ADR status '{status}' for {adr_id}"
        )
        assert link_label == link_name, (
            f"index label/target mismatch for {adr_id}: {link_label} vs {link_name}"
        )
        assert link_name in adr_by_filename, (
            f"index references missing ADR file for {adr_id}: {link_target}"
        )
        indexed_filenames.add(link_name)

        adr_path = adr_by_filename[link_name]
        file_adr_id, file_status = _parse_adr_control(adr_path)
        assert file_adr_id == adr_id, (
            f"index/file ADR ID mismatch for {link_name}: "
            f"index={adr_id} file={file_adr_id}"
        )
        assert status == file_status, (
            f"index/file status mismatch for {link_name}: "
            f"index={status} file={file_status}"
        )

    assert indexed_filenames == set(adr_by_filename), (
        "ADR index is missing ADR files or contains stale entries: "
        f"indexed={sorted(indexed_filenames)} files={sorted(adr_by_filename)}"
    )


def test_adr_markdown_links_are_not_broken() -> None:
    paths = (ADR_INDEX,) + _adr_files()
    for path in paths:
        for target in _iter_relative_links(path):
            resolved = _resolve_markdown_target(path.parent, target)
            assert resolved.exists(), (
                f"broken markdown link in {path.as_posix()}: "
                f"{target} -> {resolved.as_posix()}"
            )


def test_adr_source_file_references_exist() -> None:
    for path in _adr_files():
        text = _read_text(path)
        for match in _SOURCE_REF_PATTERN.finditer(text):
            raw = match.group(0).rstrip(".,)")
            resolved = (ROOT / raw).resolve()
            assert resolved.exists(), (
                f"stale source path in {path.as_posix()}: "
                f"{raw} -> {resolved.as_posix()}"
            )


def test_adr_0012_is_not_active_governance() -> None:
    adr_0012 = ADR_ROOT / "adr_0012_rewrite_roadmap_and_fresh_start_plan.md"
    _, status = _parse_adr_control(adr_0012)
    assert status != "Accepted", "ADR-0012 must not be an active accepted decision"
