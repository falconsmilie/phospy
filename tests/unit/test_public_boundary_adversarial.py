from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
VERIFIER_PATH = ROOT / "scripts" / "verify_distribution_artifact.py"


def _load_verifier() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "verify_distribution_artifact_public_boundary_under_test",
        VERIFIER_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _identity_context(verifier: ModuleType, tmp_path: Path):
    return verifier.VerificationContext(
        artifact_kind="wheel",
        artifact_path=tmp_path / "phospy-1.6.0-py3-none-any.whl",
        build_manifest_path=tmp_path / "build-manifest.json",
        repository_root=ROOT,
        report_json=tmp_path / "report.json",
    )


def test_public_boundary_adversarial_registry_reports_required_outcomes(
    tmp_path: Path,
) -> None:
    verifier = _load_verifier()

    details = verifier._check_public_boundary_integrity(
        _identity_context(verifier, tmp_path)
    )

    assert tuple(details["outcomes"]) == (
        "public-signature-boundary",
        "dataset-provenance-binding",
        "public-dataframe-ownership",
        "public-json-immutability",
    )
    assert {
        name: outcome["status"] for name, outcome in details["outcomes"].items()
    } == {name: "pass" for name in details["outcomes"]}
    assert (
        "AnalysisReadyPhosphoDataset.phospho/site_metadata"
        in (details["outcomes"]["public-dataframe-ownership"]["probed_classes"])
    )
    assert (
        "ResultCaveat.details"
        in (details["outcomes"]["public-json-immutability"]["probed_fields"])
    )


def test_signature_probe_detects_injected_private_ownership_switch() -> None:
    verifier = _load_verifier()

    class UnsafeExport:
        def __init__(self, phospho: object, _assume_owned: bool = False) -> None:
            self.phospho = phospho
            self.assume_owned = _assume_owned

    assert verifier._forbidden_public_parameters(UnsafeExport) == ["_assume_owned"]


def test_dataset_provenance_probe_rejects_unsafe_stale_fingerprint_acceptance() -> None:
    verifier = _load_verifier()
    dataset = verifier._trusted_dataset(
        identity_details={"fixture": "unsafe-stale-provenance-double"}
    )

    def unsafe_constructor(**kwargs: object) -> object:
        return object()

    with pytest.raises(verifier.VerificationError, match="DatasetValidationError"):
        verifier._expect_public_dataset_constructor_rejects_stale_provenance(
            dataset,
            constructor=unsafe_constructor,
        )


def test_dataframe_ownership_probe_detects_caller_owned_alias() -> None:
    verifier = _load_verifier()

    class UnsafeResultDouble:
        def __init__(self, table: pd.DataFrame) -> None:
            self.table = table

    case = verifier._PublicFrameOwnerCase(
        name="UnsafeResultDouble.table",
        make_numeric_source=lambda: pd.DataFrame({"value": [1.0, 2.0]}),
        construct_from_numeric=UnsafeResultDouble,
        observe_numeric=lambda owner: owner.table,
    )

    with pytest.raises(verifier.VerificationError, match="UnsafeResultDouble.table"):
        verifier._assert_numeric_frame_owner_case_isolated(case)


def test_json_immutability_probe_detects_shallow_caveat_details() -> None:
    verifier = _load_verifier()

    class UnsafeCaveatDouble:
        def __init__(self, details: dict[str, object]) -> None:
            self.details = dict(details)

    case = verifier._JsonImmutabilityCase(
        name="UnsafeCaveatDouble.details",
        construct=UnsafeCaveatDouble,
        observe=lambda owner: owner.details,
    )

    with pytest.raises(verifier.VerificationError, match="UnsafeCaveatDouble.details"):
        verifier._assert_json_owner_case_isolated(case)
