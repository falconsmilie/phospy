from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _user_facing_doc_paths() -> tuple[Path, ...]:
    docs: list[Path] = [
        ROOT / "README.md",
        ROOT / "docs" / "index.md",
        ROOT / "docs" / "api.md",
        ROOT / "docs" / "cli.md",
        ROOT / "docs" / "validation.md",
        ROOT / "docs" / "output_bundles.md",
    ]
    docs.extend(sorted((ROOT / "docs" / "getting-started").glob("*.md")))
    docs.extend(sorted((ROOT / "docs" / "concepts").glob("*.md")))
    docs.extend(sorted((ROOT / "docs" / "learning-paths").glob("*.md")))
    return tuple(dict.fromkeys(path.resolve() for path in docs))


_BANNED_WORDING_PATTERNS = (
    (re.compile(r"\blegacy\b", re.IGNORECASE), "legacy"),
    (re.compile(r"\brewrite\b", re.IGNORECASE), "rewrite"),
    (re.compile(r"\bmigration\b", re.IGNORECASE), "migration"),
    (re.compile(r"\bdeprecated\b", re.IGNORECASE), "deprecated"),
    (re.compile(r"\bbackwards?\b", re.IGNORECASE), "backwards"),
    (re.compile(r"\bprerelease\b", re.IGNORECASE), "prerelease"),
    (
        re.compile(r"\bcanonical\s+namespace\b", re.IGNORECASE),
        "canonical namespace",
    ),
)

_REMOVED_ALIAS_TOKENS = (
    "clustering_backend",
    "max_exact_clustering_sites",
    "prediction-ensemble-size",
    "ensemble_size",
    "ratio_to_total",
    "signalome_cutoff",
    "kinase_network_policy",
)


def test_user_facing_docs_reject_legacy_wording_and_removed_aliases() -> None:
    failures: list[str] = []
    for path in _user_facing_doc_paths():
        text = path.read_text(encoding="utf-8")
        for pattern, label in _BANNED_WORDING_PATTERNS:
            if pattern.search(text) is None:
                continue
            relative = path.relative_to(ROOT)
            failures.append(f"{relative}: found banned wording '{label}'")
        lower_text = text.lower()
        for token in _REMOVED_ALIAS_TOKENS:
            if token.lower() not in lower_text:
                continue
            relative = path.relative_to(ROOT)
            failures.append(f"{relative}: found removed alias token '{token}'")

    assert not failures, "\n".join(failures)
