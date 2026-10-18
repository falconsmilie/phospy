"""Record whether a Windows runner can safely use the SkylakeX OpenBLAS core."""

from __future__ import annotations

import argparse
import json
import os
import platform
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

REQUIRED_SKYLAKEX_FEATURES = ("AVX512F", "AVX512CD", "AVX512_SKX")


def numpy_cpu_features() -> dict[str, bool]:
    """Return NumPy's runtime hardware feature detection without overrides."""

    try:
        detected = np._core._multiarray_umath.__cpu_features__  # type: ignore[attr-defined]
    except AttributeError:
        return {}
    return {str(name): bool(value) for name, value in detected.items()}


def build_provenance_capability_disposition(
    *,
    cpu_features: Mapping[str, bool],
    system: str,
    machine: str,
    openblas_coretype: str | None,
    required_features: Sequence[str] = REQUIRED_SKYLAKEX_FEATURES,
) -> dict[str, object]:
    """Build the capability decision made before any OpenBLAS override."""

    availability = {
        feature: bool(cpu_features.get(feature, False)) for feature in required_features
    }
    missing = [feature for feature, available in availability.items() if not available]
    preflight_clean = openblas_coretype is None
    eligible = system == "Windows" and preflight_clean and not missing
    if not preflight_clean:
        disposition = "ineligible_preflight_override_present"
    elif system != "Windows":
        disposition = "ineligible_non_windows_runner"
    elif missing:
        disposition = "ineligible_missing_required_cpu_features"
    else:
        disposition = "eligible_for_controlled_exact_byte_fixture"
    return {
        "schema": "phospy-exact-provenance-capability-v1",
        "platform": {"system": system, "machine": machine},
        "feature_source": "NumPy runtime __cpu_features__",
        "required_features": list(required_features),
        "feature_availability": availability,
        "missing_features": missing,
        "openblas_coretype_before_preflight": openblas_coretype,
        "eligible": eligible,
        "disposition": disposition,
    }


def record_provenance_capability(output: Path) -> dict[str, object]:
    record = build_provenance_capability_disposition(
        cpu_features=numpy_cpu_features(),
        system=platform.system(),
        machine=platform.machine(),
        openblas_coretype=os.environ.get("OPENBLAS_CORETYPE"),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    record = record_provenance_capability(args.output)
    print(json.dumps(record, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = [
    "REQUIRED_SKYLAKEX_FEATURES",
    "build_provenance_capability_disposition",
    "numpy_cpu_features",
    "record_provenance_capability",
]
