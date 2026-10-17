from __future__ import annotations

from pathlib import Path

from scripts import verify_distribution_git_blobs as MODULE

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/verify_distribution_git_blobs.py"


def _problems(observed: dict[str, bytes]) -> tuple[str, ...]:
    return MODULE.compare_production_files(
        {
            "src/phospy/__init__.py": b"root\n",
            "src/phospy/data/resource.json": b'{"schema":1}\n',
        },
        observed,
    )


def test_distribution_git_blob_comparison_accepts_exact_paths_and_bytes() -> None:
    assert (
        _problems(
            {
                "src/phospy/__init__.py": b"root\n",
                "src/phospy/data/resource.json": b'{"schema":1}\n',
            }
        )
        == ()
    )


def test_distribution_git_blob_comparison_rejects_missing_member() -> None:
    assert _problems({"src/phospy/__init__.py": b"root\n"}) == (
        "missing production member: src/phospy/data/resource.json",
    )


def test_distribution_git_blob_comparison_rejects_additional_member() -> None:
    problems = _problems(
        {
            "src/phospy/__init__.py": b"root\n",
            "src/phospy/data/resource.json": b'{"schema":1}\n',
            "src/phospy/data/untracked.json": b"{}\n",
        }
    )
    assert problems == ("additional production member: src/phospy/data/untracked.json",)


def test_distribution_git_blob_comparison_rejects_renamed_member() -> None:
    problems = _problems(
        {
            "src/phospy/__init__.py": b"root\n",
            "src/phospy/data/RESOURCE.json": b'{"schema":1}\n',
        }
    )
    assert problems == (
        "renamed production member: src/phospy/data/resource.json -> src/phospy/data/RESOURCE.json",
    )


def test_distribution_git_blob_comparison_rejects_normalised_newlines() -> None:
    problems = _problems(
        {
            "src/phospy/__init__.py": b"root\r\n",
            "src/phospy/data/resource.json": b'{"schema":1}\n',
        }
    )
    assert problems == ("byte-different production member: src/phospy/__init__.py",)


def test_distribution_git_blob_comparison_rejects_other_byte_difference() -> None:
    problems = _problems(
        {
            "src/phospy/__init__.py": b"changed\n",
            "src/phospy/data/resource.json": b'{"schema":1}\n',
        }
    )
    assert problems == ("byte-different production member: src/phospy/__init__.py",)


def test_distribution_git_blob_verifier_declares_archive_scope() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    assert "outer archive reproducibility" in source
    assert "timestamps" in source
    assert "compression" in source
