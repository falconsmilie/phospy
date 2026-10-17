from __future__ import annotations

from pathlib import Path

from scripts.check_repository_hygiene import hygiene_problems

ROOT = Path(__file__).resolve().parents[2]


def test_repository_hygiene_policy() -> None:
    assert hygiene_problems(ROOT) == []
