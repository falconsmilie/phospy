from __future__ import annotations

from pathlib import Path

import phospy.validation.common.dataframes as common_dataframes


def test_validation_common_dataframes_has_no_phosphosite_specific_imports() -> None:
    path = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "phospy"
        / "validation"
        / "common"
        / "dataframes.py"
    )
    source = path.read_text(encoding="utf-8")
    assert "phospy.science.sites" not in source
    assert "def require_canonical_site_index(" not in source
    assert "def require_canonical_site_series(" not in source
    assert "def require_site_identity_coherence(" not in source


def test_validation_common_dataframes_does_not_export_phosphosite_validators() -> None:
    assert not hasattr(common_dataframes, "require_canonical_site_index")
    assert not hasattr(common_dataframes, "require_canonical_site_series")
    assert not hasattr(common_dataframes, "require_site_identity_coherence")
