from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_user_facing_docs_are_utf8_readable() -> None:
    for path in sorted((ROOT / "docs").rglob("*.md")):
        _ = path.read_text(encoding="utf-8")
