from __future__ import annotations

from tools.testing.provenance_capability import (
    REQUIRED_SKYLAKEX_FEATURES,
    build_provenance_capability_disposition,
)


def test_exact_provenance_capability_requires_every_skylakex_feature() -> None:
    record = build_provenance_capability_disposition(
        cpu_features={feature: True for feature in REQUIRED_SKYLAKEX_FEATURES},
        system="Windows",
        machine="AMD64",
        openblas_coretype=None,
    )

    assert record["eligible"] is True
    assert record["missing_features"] == []
    assert record["disposition"] == "eligible_for_controlled_exact_byte_fixture"


def test_exact_provenance_capability_records_missing_hardware() -> None:
    record = build_provenance_capability_disposition(
        cpu_features={"AVX512F": True, "AVX512CD": False, "AVX512_SKX": False},
        system="Windows",
        machine="AMD64",
        openblas_coretype=None,
    )

    assert record["eligible"] is False
    assert record["missing_features"] == ["AVX512CD", "AVX512_SKX"]
    assert record["disposition"] == "ineligible_missing_required_cpu_features"


def test_exact_provenance_capability_rejects_preexisting_override() -> None:
    record = build_provenance_capability_disposition(
        cpu_features={feature: True for feature in REQUIRED_SKYLAKEX_FEATURES},
        system="Windows",
        machine="AMD64",
        openblas_coretype="SkylakeX",
    )

    assert record["eligible"] is False
    assert record["disposition"] == "ineligible_preflight_override_present"
