from __future__ import annotations

import argparse
import copy
import hashlib
import importlib
import inspect
import json
import math
import platform
import sys
import time
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from importlib import metadata, resources
from pathlib import Path, PurePosixPath
from typing import Any

REPORT_SCHEMA = "phospy.artifact-verification/v1"
BUILD_MANIFEST_SCHEMA = "phospy.build-manifest/v1"
SUPPORTED_ARTIFACT_KINDS = ("wheel", "sdist")
SUMMARY_CHECK_KEYS_BY_DETAIL = {
    "artifact-manifest-binding": ("artifact_manifest_binding",),
    "installed-package-identity": ("installed_import_origin", "package_metadata"),
    "packaged-scientific-resources": ("scientific_resources",),
    "public-boundary-integrity": (
        "public_signature_boundary",
        "dataset_provenance_binding",
        "public_dataframe_ownership",
        "public_json_immutability",
    ),
    "corrected-construction-and-provenance-path": (
        "trusted_construction",
        "provenance_immutability",
    ),
    "corrected-derived-and-ownership-path": (
        "derived_lineage",
        "dataframe_ownership",
    ),
    "corrected-differential-path": ("differential_execution",),
}
PUBLIC_BOUNDARY_OUTCOME_NAMES = (
    "public-signature-boundary",
    "dataset-provenance-binding",
    "public-dataframe-ownership",
    "public-json-immutability",
)
_FORBIDDEN_PUBLIC_PARAMETER_EXACT = frozenset(
    {
        "_assume_owned",
        "assume_owned",
        "_owned",
        "owned",
        "copy",
        "copy_input",
        "copy_inputs",
        "copy_data",
        "no_copy",
        "skip_copy",
        "ownership_token",
        "transfer_token",
        "skip_validation",
        "_skip_validation",
        "disable_validation",
        "bypass_validation",
        "skip_fingerprint",
        "disable_fingerprint",
        "bypass_fingerprint",
        "trust_fingerprint",
        "suppress_warnings",
        "ignore_warnings",
        "disable_warnings",
        "silence_warnings",
        "validator",
        "request_validator",
        "source_validator",
        "config_validator",
        "interpreter",
        "executor",
        "path_reader",
        "source_reader",
        "batch_correction_runner",
        "internal_view",
    }
)
_FORBIDDEN_PUBLIC_PARAMETER_FRAGMENTS = (
    "assume_own",
    "owned_input",
    "ownership",
    "copy_bypass",
    "validation_bypass",
    "fingerprint_bypass",
    "suppress_warning",
    "ignore_warning",
    "disable_warning",
    "silence_warning",
    "internal_view",
)
_FORBIDDEN_PUBLIC_EXPORT_FRAGMENTS = (
    "Validator",
    "Interpreter",
    "Executor",
    "InternalView",
    "Ownership",
    "OwnershipToken",
    "OwnedFactory",
    "AssumeOwned",
)
_PUBLIC_FRAME_OBJECT_PAYLOAD_COLUMN = "ownership_payload"
_PUBLIC_FRAME_OBJECT_PAYLOAD_STATE = {
    "list": ("list-start",),
    "dict": ("dict-start",),
    "array": (1.0, 2.0),
    "set": ("set-start",),
    "nested_array": (3.0, 4.0),
    "nested_set": ("nested-set-start",),
    "nested_list": ("nested-list-start",),
}


class VerificationError(RuntimeError):
    """Raised when installed-artifact verification detects an invalid artifact."""


@dataclass(slots=True)
class VerificationContext:
    artifact_kind: str
    artifact_path: Path
    build_manifest_path: Path
    repository_root: Path
    report_json: Path
    package_root: Path | None = None
    distribution_version: str | None = None
    source_identity_digest: str | None = None
    artifact_filename: str | None = None
    artifact_sha256: str | None = None
    build_manifest_package_version: str | None = None


CheckFunction = Callable[[VerificationContext], Mapping[str, object]]


@dataclass(frozen=True, slots=True)
class CheckResult:
    name: str
    status: str
    duration_seconds: float
    details: Mapping[str, object]
    message: str | None = None

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "name": self.name,
            "status": self.status,
            "duration_seconds": round(float(self.duration_seconds), 6),
            "details": _json_safe(self.details),
        }
        if self.message is not None:
            payload["message"] = self.message
        return payload


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    context = VerificationContext(
        artifact_kind=args.artifact_kind,
        artifact_path=Path(args.artifact_path).resolve(),
        build_manifest_path=Path(args.build_manifest).resolve(),
        repository_root=Path(args.repository_root).resolve(),
        report_json=Path(args.report_json).resolve(),
    )
    report = _build_report(context)
    _write_report(context.report_json, report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "artifact": report["artifact"],
                "report_json": context.report_json.as_posix(),
            },
            sort_keys=True,
        )
    )
    return 0 if report["status"] == "success" else 1


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify an installed PhosPy distribution artifact."
    )
    parser.add_argument(
        "--artifact-kind",
        required=True,
        choices=SUPPORTED_ARTIFACT_KINDS,
        help="Distribution artifact kind under verification.",
    )
    parser.add_argument(
        "--artifact-path",
        required=True,
        help="Path to the wheel or sdist artifact that was installed.",
    )
    parser.add_argument(
        "--build-manifest",
        required=True,
        help="Path to the build manifest binding source identity to artifacts.",
    )
    parser.add_argument(
        "--repository-root",
        required=True,
        help="Source repository root; imports from this tree are rejected.",
    )
    parser.add_argument(
        "--report-json",
        required=True,
        help="Path where the machine-readable verification report is written.",
    )
    return parser.parse_args(argv)


def _build_report(context: VerificationContext) -> dict[str, object]:
    check_results: list[CheckResult] = []
    identity_failed = False

    for check_name, check_function in _CHECKS:
        if identity_failed:
            check_results.append(
                CheckResult(
                    name=check_name,
                    status="skipped",
                    duration_seconds=0.0,
                    details={},
                    message=(
                        "installed-package-identity failed; runtime checks were "
                        "not run against an ambiguous import"
                    ),
                )
            )
            continue

        started_at = time.perf_counter()
        try:
            details = check_function(context)
            _assert_loaded_phospy_modules_not_from_repository(context)
            check_results.append(
                CheckResult(
                    name=check_name,
                    status="pass",
                    duration_seconds=time.perf_counter() - started_at,
                    details=details,
                )
            )
        except VerificationError as exc:
            if check_name == "installed-package-identity":
                identity_failed = True
            check_results.append(
                CheckResult(
                    name=check_name,
                    status="fail",
                    duration_seconds=time.perf_counter() - started_at,
                    details={},
                    message=str(exc),
                )
            )
        except Exception as exc:
            if check_name == "installed-package-identity":
                identity_failed = True
            check_results.append(
                CheckResult(
                    name=check_name,
                    status="fail",
                    duration_seconds=time.perf_counter() - started_at,
                    details={},
                    message=f"{type(exc).__name__}: {exc}",
                )
            )

    status = (
        "success" if all(item.status == "pass" for item in check_results) else "failure"
    )
    dependency_snapshot = _dependency_snapshot()
    return {
        "schema": REPORT_SCHEMA,
        "status": status,
        "source_identity_digest": context.source_identity_digest,
        "artifact": {
            "kind": context.artifact_kind,
            "filename": context.artifact_filename or context.artifact_path.name,
            "sha256": context.artifact_sha256,
        },
        "package": {
            "name": "phospy",
            "version": context.distribution_version
            or context.build_manifest_package_version,
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
            "isolated": bool(sys.flags.isolated),
            "dependency_snapshot_sha256": _payload_sha256(dependency_snapshot),
            "dependency_snapshot": dependency_snapshot,
        },
        "checks": _summarize_check_statuses(check_results),
        "check_details": [item.to_payload() for item in check_results],
        "repository_root": context.repository_root.as_posix(),
        "working_directory": Path.cwd().resolve().as_posix(),
        "distribution": {
            "name": "phospy",
            "version": context.distribution_version,
            "package_root": (
                None
                if context.package_root is None
                else context.package_root.as_posix()
            ),
        },
    }


def _write_report(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _check_artifact_manifest_binding(
    context: VerificationContext,
) -> Mapping[str, object]:
    _require_artifact_filename_matches_kind(
        context.artifact_path, context.artifact_kind
    )
    _require(
        context.artifact_path.is_file(),
        f"artifact under verification is missing: {context.artifact_path.as_posix()}",
    )
    _require(
        context.build_manifest_path.is_file(),
        f"build manifest is missing: {context.build_manifest_path.as_posix()}",
    )

    artifact_sha256 = _file_sha256(context.artifact_path)
    manifest = _read_json_mapping(context.build_manifest_path)
    _require_text_equals(
        manifest.get("schema"),
        BUILD_MANIFEST_SCHEMA,
        field_name="build manifest schema",
    )
    source_identity_digest = _required_text(
        manifest.get("source_identity_digest"),
        field_name="build manifest source_identity_digest",
    )
    _require_sha256_digest(source_identity_digest, field_name="source_identity_digest")
    package_version = _required_text(
        manifest.get("package_version"),
        field_name="build manifest package_version",
    )
    artifacts = manifest.get("artifacts")
    _require(
        isinstance(artifacts, list) and len(artifacts) > 0,
        "build manifest artifacts must be a non-empty array",
    )

    matching_artifacts: list[Mapping[str, object]] = []
    for artifact in artifacts:
        _require(
            isinstance(artifact, Mapping),
            "build manifest artifacts entries must be objects",
        )
        if (
            artifact.get("kind") == context.artifact_kind
            and artifact.get("filename") == context.artifact_path.name
        ):
            matching_artifacts.append(artifact)
    _require(
        len(matching_artifacts) == 1,
        "build manifest must contain exactly one entry for "
        f"{context.artifact_kind} artifact {context.artifact_path.name!r}",
    )
    artifact_manifest = matching_artifacts[0]
    expected_sha256 = _required_text(
        artifact_manifest.get("sha256"),
        field_name=f"build manifest artifact {context.artifact_path.name}.sha256",
    ).lower()
    _require_sha256_hex(expected_sha256, field_name="artifact sha256")
    if artifact_sha256 != expected_sha256:
        raise VerificationError(
            "artifact digest does not match build manifest: "
            f"{context.artifact_path.name}: expected {expected_sha256}, "
            f"got {artifact_sha256}"
        )

    context.source_identity_digest = source_identity_digest
    context.artifact_filename = context.artifact_path.name
    context.artifact_sha256 = artifact_sha256
    context.build_manifest_package_version = package_version
    return {
        "build_manifest": context.build_manifest_path.as_posix(),
        "source_identity_digest": source_identity_digest,
        "package_version": package_version,
        "artifact": {
            "kind": context.artifact_kind,
            "filename": context.artifact_path.name,
            "sha256": artifact_sha256,
        },
    }


def _check_installed_package_identity(
    context: VerificationContext,
) -> Mapping[str, object]:
    _require_isolated_interpreter()
    _require_external_working_directory(context.repository_root)

    phospy = importlib.import_module("phospy")
    distribution = metadata.distribution("phospy")
    metadata_version = metadata.version("phospy")
    context.distribution_version = metadata_version
    runtime_version = getattr(phospy, "__version__", None)
    _require(
        isinstance(runtime_version, str) and bool(runtime_version.strip()),
        "phospy runtime version is missing; expected phospy.__version__",
    )
    if runtime_version != metadata_version:
        raise VerificationError(
            "phospy runtime version does not match package metadata: "
            f"runtime={runtime_version!r}; metadata={metadata_version!r}"
        )

    package_file = _module_file(phospy, module_name="phospy")
    _reject_if_inside_repository(
        package_file,
        context.repository_root,
        description="phospy package",
    )
    distribution_init_paths = _distribution_init_paths(distribution)
    _require(
        any(
            _same_resolved_path(package_file, candidate)
            for candidate in distribution_init_paths
        ),
        "phospy.__file__ does not match importlib.metadata distribution files: "
        f"{package_file.as_posix()}",
    )

    context.package_root = package_file.parent
    if context.build_manifest_package_version is not None:
        _require(
            metadata_version == context.build_manifest_package_version,
            "package metadata version does not match build manifest: "
            f"metadata={metadata_version!r}; "
            f"manifest={context.build_manifest_package_version!r}",
        )
    _assert_loaded_phospy_modules_not_from_repository(context)
    return {
        "metadata_version": metadata_version,
        "runtime_version": runtime_version,
        "package_file": package_file.as_posix(),
        "package_root": package_file.parent.as_posix(),
        "distribution_init_files": [
            path.as_posix() for path in sorted(distribution_init_paths)
        ],
    }


def _check_packaged_scientific_resources(
    context: VerificationContext,
) -> Mapping[str, object]:
    from phospy.api import Organism, ReferencePreset
    from phospy.science.references.resolution import ReferenceResolver
    from phospy.science.references.resources import (
        available_bundled_reference_lanes,
        load_bundled_motif_scores,
        load_bundled_motif_sizes,
    )

    bundle = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=Organism.RAT,
    )
    _require(bundle.manifest is not None, "bundled reference manifest was not loaded")
    _require(
        bundle.manifest.reference_id == "l6_native",
        f"unexpected bundled reference id: {bundle.manifest.reference_id!r}",
    )
    kinase_substrate_map = bundle.kinase_substrate_map_dataframe()
    site_sequences = bundle.site_sequences_dataframe()
    _require(
        not kinase_substrate_map.empty,
        "normal bundled kinase-substrate loading path returned an empty table",
    )
    _require(
        not site_sequences.empty,
        "normal bundled site-sequence loading path returned an empty table",
    )

    motif_scores = load_bundled_motif_scores(Organism.RAT)
    motif_sizes = load_bundled_motif_sizes(Organism.RAT)
    _require(
        motif_scores is not None and not motif_scores.empty,
        "bundled motif_scores.csv resource did not load",
    )
    _require(
        motif_sizes is not None and not motif_sizes.empty,
        "bundled motif_sizes.csv resource did not load",
    )
    lanes = available_bundled_reference_lanes()
    _require(
        any(
            lane.organism is Organism.RAT and lane.bundle_id == "l6_native"
            for lane in lanes
        ),
        "bundled reference inventory does not include rat/l6_native",
    )

    bundle_root = (
        resources.files("phospy")
        .joinpath("data")
        .joinpath("reference_bundles")
        .joinpath("rat")
        .joinpath("l6_native")
    )
    manifest_resource = bundle_root.joinpath("manifest.json")
    _require(
        manifest_resource.is_file(),
        "packaged scientific resource manifest is missing: rat/l6_native/manifest.json",
    )
    manifest_text = manifest_resource.read_text(encoding="utf-8")
    manifest_payload = json.loads(manifest_text)
    digest_report = _validate_manifest_resource_digests(
        bundle_root,
        manifest_payload,
    )
    return {
        "bundle": "rat/l6_native",
        "reference_id": bundle.manifest.reference_id,
        "kinase_substrate_rows": int(len(kinase_substrate_map.index)),
        "site_sequence_rows": int(len(site_sequences.index)),
        "motif_score_shape": [
            int(motif_scores.shape[0]),
            int(motif_scores.shape[1]),
        ],
        "motif_size_count": int(motif_sizes.size),
        "manifest_sha256": hashlib.sha256(manifest_text.encode("utf-8")).hexdigest(),
        **digest_report,
    }


def _check_public_boundary_integrity(
    context: VerificationContext,
) -> Mapping[str, object]:
    outcomes: dict[str, object] = {}
    probes: tuple[tuple[str, Callable[[], Mapping[str, object]]], ...] = (
        ("public-signature-boundary", _probe_public_signature_boundary),
        ("dataset-provenance-binding", _probe_dataset_provenance_binding),
        ("public-dataframe-ownership", _probe_public_dataframe_ownership),
        ("public-json-immutability", _probe_public_json_immutability),
    )
    for outcome_name, probe in probes:
        try:
            details = probe()
        except VerificationError as exc:
            raise VerificationError(f"{outcome_name}: {exc}") from exc
        outcomes[outcome_name] = {"status": "pass", **dict(details)}
    _require(
        tuple(outcomes) == PUBLIC_BOUNDARY_OUTCOME_NAMES,
        "public-boundary-integrity outcome registry changed unexpectedly",
    )
    return {
        "outcomes": outcomes,
        "required_outcome_names": list(PUBLIC_BOUNDARY_OUTCOME_NAMES),
    }


def _probe_public_signature_boundary() -> Mapping[str, object]:
    exports = _public_boundary_exports()
    forbidden_exports = _forbidden_public_exports(exports)
    _require(
        forbidden_exports == {},
        "private boundary helpers exported from supported public namespaces: "
        + "; ".join(
            f"{name} ({reason})" for name, reason in sorted(forbidden_exports.items())
        ),
    )

    signature_offenders: dict[str, list[str]] = {}
    inspected_signatures: list[str] = []
    for exported_name, exported in sorted(exports.items()):
        if not (inspect.isclass(exported) or inspect.isfunction(exported)):
            continue
        forbidden = _forbidden_public_parameters(exported)
        inspected_signatures.append(exported_name)
        if forbidden:
            signature_offenders[exported_name] = forbidden
    _require(
        signature_offenders == {},
        "public signatures expose private boundary controls: "
        + "; ".join(
            f"{name}({', '.join(parameters)})"
            for name, parameters in sorted(signature_offenders.items())
        ),
    )
    return {
        "exported_symbol_count": len(exports),
        "inspected_signature_count": len(inspected_signatures),
        "inspected_signatures": inspected_signatures,
    }


def _public_boundary_exports() -> dict[str, object]:
    import phospy
    import phospy.api as public_api

    exports: dict[str, object] = {}
    for module_name, module in (("phospy", phospy), ("phospy.api", public_api)):
        names = getattr(module, "__all__", None)
        _require(
            isinstance(names, Sequence) and not isinstance(names, (str, bytes)),
            f"{module_name}.__all__ must be a sequence of public names",
        )
        for name in names:
            _require(
                isinstance(name, str) and name.strip(),
                f"{module_name}.__all__ contains a non-text export name",
            )
            _require(
                not name.startswith("_"),
                f"{module_name} exports private symbol {name!r}",
            )
            _require(hasattr(module, name), f"{module_name}.{name} is missing")
            exports[f"{module_name}.{name}"] = getattr(module, name)
    return exports


def _forbidden_public_exports(exports: Mapping[str, object]) -> dict[str, str]:
    offenders: dict[str, str] = {}
    for exported_name in exports:
        symbol_name = exported_name.rsplit(".", maxsplit=1)[1]
        for fragment in _FORBIDDEN_PUBLIC_EXPORT_FRAGMENTS:
            if fragment in symbol_name:
                offenders[exported_name] = fragment
                break
    return offenders


def _forbidden_public_parameters(owner: object) -> list[str]:
    try:
        parameters = inspect.signature(owner).parameters
    except (TypeError, ValueError):
        return []
    return [
        name
        for name in parameters
        if name not in {"self", "cls"} and _is_forbidden_public_parameter(name)
    ]


def _is_forbidden_public_parameter(name: str) -> bool:
    normalized = name.lower()
    if name.startswith("_") or normalized in _FORBIDDEN_PUBLIC_PARAMETER_EXACT:
        return True
    if any(
        fragment in normalized for fragment in _FORBIDDEN_PUBLIC_PARAMETER_FRAGMENTS
    ):
        return True
    if any(
        fragment in normalized for fragment in ("validator", "interpreter", "executor")
    ):
        return True
    disabling_tokens = ("skip", "disable", "bypass", "ignore", "suppress", "silence")
    if "warning" in normalized and any(
        token in normalized for token in disabling_tokens
    ):
        return True
    if "validation" in normalized and any(
        token in normalized for token in disabling_tokens
    ):
        return True
    if "fingerprint" in normalized and any(
        token in normalized for token in (*disabling_tokens, "trust")
    ):
        return True
    if "check" in normalized and any(token in normalized for token in disabling_tokens):
        return True
    return False


def _probe_dataset_provenance_binding() -> Mapping[str, object]:
    dataset = _trusted_dataset(
        identity_details={
            "columns": ["site_key"],
            "fixture": "public-boundary-integrity",
        }
    )
    message = _expect_public_dataset_constructor_rejects_stale_provenance(dataset)
    return {
        "mutated_table": "dataset.phospho",
        "stale_provenance_rejected": True,
        "rejection_message": message,
    }


def _probe_public_dataframe_ownership() -> Mapping[str, object]:
    probed_cases: list[str] = []
    object_cell_cases: list[str] = []
    for case in _public_frame_owner_cases():
        _assert_numeric_frame_owner_case_isolated(case)
        probed_cases.append(case.name)
        if (
            case.make_object_source is not None
            and case.construct_from_object is not None
            and case.observe_object_payload is not None
        ):
            _assert_object_frame_owner_case_isolated(case)
            object_cell_cases.append(case.name)
    return {
        "probed_classes": probed_cases,
        "nested_object_cell_classes": object_cell_cases,
    }


def _probe_public_json_immutability() -> Mapping[str, object]:
    probed_fields: list[str] = []
    for case in _json_immutability_cases():
        _assert_json_owner_case_isolated(case)
        probed_fields.append(case.name)
    return {"probed_fields": probed_fields}


def _check_corrected_construction_and_provenance_path(
    context: VerificationContext,
) -> Mapping[str, object]:
    from phospy.api import Organism, ReferencePreset
    from phospy.errors import ReferenceCompatibilityError
    from phospy.provenance import from_payload, to_payload
    from phospy.provenance.hashing import hash_json_payload
    from phospy.science.references.resolution import ReferenceResolver

    identity_details: dict[str, object] = {
        "columns": ["site_key"],
        "schema": {"version": 1},
    }
    dataset = _trusted_dataset(identity_details=identity_details)
    _require(dataset.provenance is not None, "trusted construction lacks provenance")
    _require(
        dataset.trusted_construction_assertions is not None,
        "trusted construction assertions were not retained",
    )
    _require(
        dataset.trusted_construction_assertions.all_required_assertions_present,
        "trusted construction assertions are incomplete",
    )

    before_payload = to_payload(dataset.provenance)
    before_hash = hash_json_payload(before_payload)
    identity_details["columns"].append("mutated-after-construction")
    identity_schema = _mapping(identity_details["schema"])
    identity_schema["version"] = 99
    _try_mutate_public_provenance_details(dataset)

    after_payload = to_payload(dataset.provenance)
    _require(
        after_payload == before_payload,
        "provenance changed after mutating caller-owned or public nested payloads",
    )
    _require(
        hash_json_payload(after_payload) == before_hash,
        "provenance digest changed after nested mutation attempt",
    )

    serialized_payload = copy.deepcopy(before_payload)
    serialized_workflow = _mapping(serialized_payload["workflow_parameters"])
    serialized_construction = _mapping(serialized_workflow["construction"])
    serialized_construction["method"] = "tampered-after-serialization"
    _require(
        to_payload(dataset.provenance) == before_payload,
        "serialized provenance payload is not detached from internal state",
    )
    _require(
        to_payload(from_payload(before_payload)) == before_payload,
        "serialized provenance does not round-trip exactly",
    )

    construction_payload = _mapping(
        _mapping(before_payload["workflow_parameters"])["construction"]
    )
    trusted_payload = _mapping(construction_payload["trusted_construction_assertions"])
    _require(
        trusted_payload["schema_version"] == 3,
        "trusted construction assertion schema version was not recorded",
    )
    _require(
        trusted_payload["missing_assertions"] == [],
        "trusted construction assertion evidence is incomplete",
    )
    _require(
        construction_payload["trusted_construction_assertion_fingerprint"]
        == dataset.trusted_construction_assertions.assertion_fingerprint,
        "trusted construction assertion fingerprint was not embedded in provenance",
    )

    _expect_raises(
        ReferenceCompatibilityError,
        lambda: ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.RAT,
        ),
        contains="dataset.organism",
    )
    stale_direct_constructor_message = (
        _expect_direct_constructor_rejects_stale_provenance(dataset)
    )
    return {
        "trusted_assertion_fingerprint": (
            dataset.trusted_construction_assertions.assertion_fingerprint
        ),
        "provenance_hash": before_hash,
        "round_trip_verified": True,
        "organism_reference_contradiction_rejected": True,
        "stale_direct_constructor_provenance_rejected": True,
        "stale_direct_constructor_message": stale_direct_constructor_message,
    }


def _check_corrected_derived_and_ownership_path(
    context: VerificationContext,
) -> Mapping[str, object]:
    from phospy.api import (
        AnalysisReadyPhosphoDataset,
        Organism,
        PhosphositeImportResult,
    )
    from phospy.errors import DatasetValidationError
    from phospy.provenance import DerivedSampleMapping
    from phospy.provenance.derived_quantitative import (
        DerivedQuantitativeDataProvenance,
        build_derived_quantitative_run_provenance,
    )
    from phospy.science.datasets.derived_quantitative import (
        DerivedAnalysisReadyPhosphoDataset,
    )

    fixture = _derived_fixture()
    _require(
        "_assume_owned"
        not in inspect.signature(AnalysisReadyPhosphoDataset).parameters,
        "AnalysisReadyPhosphoDataset exposes public ownership transfer",
    )
    _require(
        "_assume_owned" not in inspect.signature(PhosphositeImportResult).parameters,
        "PhosphositeImportResult exposes public ownership transfer",
    )
    _expect_raises(
        TypeError,
        lambda: AnalysisReadyPhosphoDataset(
            phospho=fixture.phospho,
            site_metadata=fixture.site_metadata,
            organism=Organism.RAT,
            intensity_scale_state=_supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=_supported_linear_processing_state(has_total_matrix=False),
            _assume_owned=True,
        ),
        contains="_assume_owned",
    )
    _expect_raises(
        TypeError,
        lambda: PhosphositeImportResult(
            phospho_matrix_candidate=fixture.phospho,
            site_metadata_candidate=fixture.site_metadata,
            sample_column_mapping={"bio_a": "bio_a"},
            _assume_owned=True,
        ),
        contains="_assume_owned",
    )
    provenance = build_derived_quantitative_run_provenance(lineage=fixture.lineage)
    dataset = DerivedAnalysisReadyPhosphoDataset(
        phospho=fixture.phospho,
        site_metadata=fixture.site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=_supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=_supported_linear_processing_state(has_total_matrix=False),
        parent_state=fixture.parent_state,
        derived_lineage=fixture.lineage,
        provenance=provenance,
    )

    retained_lineage = DerivedQuantitativeDataProvenance.from_payload(
        _mapping(dataset.provenance.workflow_parameters)["derived_quantitative_data"]
    )
    _require(
        retained_lineage.to_payload() == dataset.derived_lineage.to_payload(),
        "derived dataset provenance did not retain source lineage",
    )
    _require(
        dataset.derived_lineage.sample_groups()
        == (
            ("bio_a", ("bio_a_t1", "bio_a_t2")),
            ("bio_b", ("bio_b_t1", "bio_b_t2")),
        ),
        "derived dataset source sample lineage was not retained",
    )

    bad_mapping = (
        DerivedSampleMapping(
            output_sample_id="bio_a",
            input_sample_ids=("fabricated_input", "bio_a_t2"),
            condition="A",
            biological_replicate_id="bio_a",
            technical_replicate_ids=("t1", "t2"),
        ),
        fixture.lineage.sample_mapping[1],
    )
    bad_lineage = DerivedQuantitativeDataProvenance(
        derivation_type=fixture.lineage.derivation_type,
        parent_dataset_type=fixture.lineage.parent_dataset_type,
        derived_dataset_type=fixture.lineage.derived_dataset_type,
        parent_dataset_fingerprints=fixture.lineage.parent_dataset_fingerprints,
        derived_dataset_fingerprints=fixture.lineage.derived_dataset_fingerprints,
        sample_mapping=bad_mapping,
        aggregation_method=fixture.lineage.aggregation_method,
        input_intensity_scale=fixture.lineage.input_intensity_scale,
        output_intensity_scale=fixture.lineage.output_intensity_scale,
        quantitative_meaning=fixture.lineage.quantitative_meaning,
        missingness_policy=fixture.lineage.missingness_policy,
        matrices_transformed=fixture.lineage.matrices_transformed,
        implementation=fixture.lineage.implementation,
        implementation_version=fixture.lineage.implementation_version,
    )
    bad_provenance = build_derived_quantitative_run_provenance(lineage=bad_lineage)
    _expect_raises(
        DatasetValidationError,
        lambda: DerivedAnalysisReadyPhosphoDataset(
            phospho=fixture.phospho.copy(deep=True),
            site_metadata=fixture.site_metadata.copy(deep=True),
            organism=Organism.RAT,
            intensity_scale_state=_supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=_supported_linear_processing_state(has_total_matrix=False),
            parent_state=fixture.parent_state,
            derived_lineage=bad_lineage,
            provenance=bad_provenance,
        ),
        contains="input_sample_ids",
    )

    original_phospho_value = float(dataset.phospho.iloc[0, 0])
    original_site_key = str(dataset.site_metadata.iloc[0]["site_key"])
    fixture.phospho.iloc[0, 0] = 999.0
    fixture.site_metadata.iloc[0, fixture.site_metadata.columns.get_loc("site_key")] = (
        "mutated"
    )
    _require(
        float(dataset.phospho.iloc[0, 0]) == original_phospho_value,
        "mutating caller-owned phospho input changed the derived dataset",
    )
    _require(
        str(dataset.site_metadata.iloc[0]["site_key"]) == original_site_key,
        "mutating caller-owned site metadata changed the derived dataset",
    )

    exported_phospho = dataset.phospho
    exported_site_metadata = dataset.site_metadata
    exported_phospho.iloc[0, 0] = 777.0
    exported_site_metadata.iloc[
        0, exported_site_metadata.columns.get_loc("site_key")
    ] = "export-mutated"
    _require(
        float(dataset.phospho.iloc[0, 0]) == original_phospho_value,
        "mutating public phospho export changed the derived dataset",
    )
    _require(
        str(dataset.site_metadata.iloc[0]["site_key"]) == original_site_key,
        "mutating public site metadata export changed the derived dataset",
    )
    return {
        "lineage_hash": dataset.derived_lineage.lineage_hash_value,
        "sample_groups": [
            [output, list(inputs)]
            for output, inputs in dataset.derived_lineage.sample_groups()
        ],
        "fabricated_input_lineage_rejected": True,
        "caller_owned_inputs_are_isolated": True,
        "public_exports_are_isolated": True,
        "public_ownership_transfer_parameters_absent": True,
        "public_ownership_aliasing_rejected": True,
    }


def _check_corrected_differential_path(
    context: VerificationContext,
) -> Mapping[str, object]:
    import numpy as np

    from phospy.api import (
        Contrast,
        DifferentialAnalysisRequest,
        DifferentialAnalysisWorkflow,
        ExperimentalDesign,
        SampleDesignRecord,
    )
    from phospy.science.differential.linear_model import decompose_differential_design

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_differential_dataset(),
            design=ExperimentalDesign(
                samples=(
                    SampleDesignRecord(sample_id="A_1", condition="A"),
                    SampleDesignRecord(sample_id="A_2", condition="A"),
                    SampleDesignRecord(sample_id="B_1", condition="B"),
                    SampleDesignRecord(sample_id="B_2", condition="B"),
                )
            ),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
        )
    )
    _require(
        result.policy_provenance is not None,
        "differential result lacks policy provenance",
    )
    table = result.table_for("B_vs_A")
    first_site = table.index[0]
    _assert_close(
        float(table.loc[first_site, "logFC"]),
        1.0,
        label="B_vs_A first-site logFC",
    )

    design = result.policy_provenance.design
    diagnostics = result.diagnostics.to_payload()
    _require(design.rank == 2, f"unexpected differential design rank: {design.rank}")
    _assert_close(
        float(design.residual_degrees_of_freedom),
        2.0,
        label="design residual degrees of freedom",
    )
    _require(diagnostics["rank"] == design.rank, "diagnostic rank does not match fit")
    _assert_close(
        float(diagnostics["residual_degrees_of_freedom"]),
        float(design.residual_degrees_of_freedom),
        label="diagnostic residual degrees of freedom",
    )
    _require(
        diagnostics["decomposition_method"] == design.decomposition_method,
        "diagnostic decomposition method does not match fit provenance",
    )

    design_values = np.array(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0],
        ],
        dtype=float,
    )
    decomposition = decompose_differential_design(design_values)
    _require(
        design.decomposition_method == decomposition.decomposition_method,
        "recorded design decomposition does not match deterministic fit",
    )
    _require(
        design.solver == decomposition.solver,
        "recorded design solver does not match deterministic fit",
    )
    _assert_close(
        float(design.condition_number),
        float(decomposition.condition_number),
        label="design condition number",
    )
    _assert_sequence_close(
        design.singular_values,
        decomposition.singular_values,
        label="design singular values",
    )

    covariance = decomposition.coefficient_covariance
    contrast_covariance = decomposition.contrast_covariance(
        np.array([[-1.0], [1.0]], dtype=float)
    )
    _require(
        bool(np.isfinite(covariance).all()),
        "differential coefficient covariance contains non-finite values",
    )
    _require(
        bool(np.allclose(covariance, covariance.T, rtol=0.0, atol=1.0e-12)),
        "differential coefficient covariance is not symmetric",
    )
    _require(
        bool(np.isfinite(contrast_covariance).all()),
        "differential contrast covariance contains non-finite values",
    )
    _require(
        bool(
            np.allclose(
                contrast_covariance,
                contrast_covariance.T,
                rtol=0.0,
                atol=1.0e-12,
            )
        ),
        "differential contrast covariance is not symmetric",
    )
    _require(
        bool(np.isfinite(result.posterior_residual_variance_series()).all()),
        "posterior residual variance contains non-finite values",
    )
    return {
        "contrast": "B_vs_A",
        "first_site": str(first_site),
        "first_site_logfc": float(table.loc[first_site, "logFC"]),
        "rank": int(design.rank),
        "residual_degrees_of_freedom": float(design.residual_degrees_of_freedom),
        "decomposition_method": design.decomposition_method,
        "coefficient_covariance": covariance.tolist(),
        "contrast_covariance": contrast_covariance.tolist(),
    }


def _validate_manifest_resource_digests(
    bundle_root: Any,
    manifest_payload: Mapping[str, object],
) -> dict[str, object]:
    files = manifest_payload.get("files")
    _require(
        isinstance(files, Sequence) and not isinstance(files, (str, bytes, bytearray)),
        "packaged scientific resource manifest has no files array",
    )
    roles: list[str] = []
    verified_files: list[str] = []
    for position, entry in enumerate(files):
        _require(
            isinstance(entry, Mapping),
            f"packaged scientific resource manifest files[{position}] is not an object",
        )
        relative_path = _required_manifest_text(
            entry.get("relative_path"),
            field_name=f"files[{position}].relative_path",
        )
        expected_sha256 = _required_manifest_text(
            entry.get("sha256"),
            field_name=f"files[{position}].sha256",
        )
        role = _required_manifest_text(
            entry.get("role", "unknown"),
            field_name=f"files[{position}].role",
        )
        resource = _join_resource_path(bundle_root, relative_path)
        _require(
            resource.is_file(),
            f"packaged scientific resource is missing: {relative_path}",
        )
        observed_sha256 = hashlib.sha256(resource.read_bytes()).hexdigest()
        if observed_sha256 != expected_sha256:
            raise VerificationError(
                "packaged scientific resource digest mismatch for "
                f"{relative_path}: expected {expected_sha256}, got {observed_sha256}"
            )
        roles.append(role)
        verified_files.append(relative_path)

    required_roles = {
        "kinase_substrate",
        "site_sequences",
        "motif_scores",
        "motif_sizes",
        "attribution",
    }
    missing_roles = sorted(required_roles - set(roles))
    _require(
        missing_roles == [],
        "packaged scientific resource manifest is missing required roles: "
        + ", ".join(missing_roles),
    )
    return {
        "manifest_file_count": len(verified_files),
        "verified_files": verified_files,
        "verified_roles": sorted(set(roles)),
    }


@dataclass(frozen=True, slots=True)
class _DerivedFixture:
    phospho: Any
    site_metadata: Any
    parent_state: Any
    lineage: Any


@dataclass(frozen=True, slots=True)
class _PublicFrameOwnerCase:
    name: str
    make_numeric_source: Callable[[], Any]
    construct_from_numeric: Callable[[Any], object]
    observe_numeric: Callable[[object], Any]
    make_object_source: Callable[[object], Any] | None = None
    construct_from_object: Callable[[Any], object] | None = None
    observe_object_payload: Callable[[object], object] | None = None


@dataclass(frozen=True, slots=True)
class _JsonImmutabilityCase:
    name: str
    construct: Callable[[Mapping[str, object]], object]
    observe: Callable[[object], object]


def _public_frame_owner_cases() -> tuple[_PublicFrameOwnerCase, ...]:
    return (
        _PublicFrameOwnerCase(
            name="AnalysisReadyPhosphoDataset.phospho/site_metadata",
            make_numeric_source=_small_phospho,
            construct_from_numeric=lambda frame: _analysis_ready_dataset_from_frames(
                phospho=frame,
                site_metadata=_sample_site_metadata(),
            ),
            observe_numeric=lambda owner: owner.to_dataframe(),
            make_object_source=lambda payload: _object_payload_frame(
                _sample_site_metadata(),
                payload,
            ),
            construct_from_object=lambda frame: _analysis_ready_dataset_from_frames(
                phospho=_small_phospho(),
                site_metadata=frame,
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.site_metadata_dataframe()
            ),
        ),
        _PublicFrameOwnerCase(
            name="PhosphositeImportResult.phospho_matrix_candidate/site_metadata_candidate",
            make_numeric_source=_small_phospho,
            construct_from_numeric=lambda frame: _import_result_from_frames(
                phospho=frame,
                site_metadata=_sample_site_metadata(),
            ),
            observe_numeric=lambda owner: owner.phospho_matrix_candidate,
            make_object_source=lambda payload: _object_payload_frame(
                _sample_site_metadata(),
                payload,
            ),
            construct_from_object=lambda frame: _import_result_from_frames(
                phospho=_small_phospho(),
                site_metadata=frame,
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.site_metadata_candidate
            ),
        ),
        _PublicFrameOwnerCase(
            name="DifferentialAnalysisResult.contrast_tables",
            make_numeric_source=_differential_result_table,
            construct_from_numeric=_differential_result_from_table,
            observe_numeric=lambda owner: owner.table_for("B_vs_A"),
            make_object_source=lambda payload: _object_payload_frame(
                _differential_result_table(),
                payload,
            ),
            construct_from_object=_differential_result_from_table,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.table_for("B_vs_A")
            ),
        ),
        _PublicFrameOwnerCase(
            name="KinaseActivityResult.activity_matrix/target_table",
            make_numeric_source=_activity_matrix,
            construct_from_numeric=_activity_result_from_matrix,
            observe_numeric=lambda owner: owner.activity_matrix,
            make_object_source=lambda payload: _object_payload_frame(
                _activity_target_table(),
                payload,
            ),
            construct_from_object=_activity_result_from_target_table,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.target_table
            ),
        ),
        _PublicFrameOwnerCase(
            name="KinaseScoringResult.profile_scores",
            make_numeric_source=_kinase_score_matrix,
            construct_from_numeric=_kinase_scoring_result_from_matrix,
            observe_numeric=lambda owner: owner.profile_scores,
        ),
        _PublicFrameOwnerCase(
            name="KinasePredictionResult.pred_mat/substrate_list",
            make_numeric_source=_prediction_matrix,
            construct_from_numeric=_kinase_prediction_result_from_matrix,
            observe_numeric=lambda owner: owner.pred_mat,
            make_object_source=lambda payload: _object_payload_frame(
                _substrate_list(),
                payload,
            ),
            construct_from_object=_kinase_prediction_result_from_substrate_list,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.substrate_list
            ),
        ),
        _PublicFrameOwnerCase(
            name="KinaseWorkflowResult.substrate_contributions",
            make_numeric_source=_kinase_substrate_contribution_table,
            construct_from_numeric=_kinase_workflow_result_from_contributions,
            observe_numeric=lambda owner: owner.substrate_contributions,
        ),
        _PublicFrameOwnerCase(
            name="SignalomeWorkflowResult.expanded_signalome",
            make_numeric_source=_signalome_expanded_table,
            construct_from_numeric=_signalome_result_from_expanded_table,
            observe_numeric=lambda owner: owner.to_dataframe(),
            make_object_source=lambda payload: _object_payload_frame(
                _signalome_expanded_table(),
                payload,
            ),
            construct_from_object=_signalome_result_from_expanded_table,
            observe_object_payload=lambda owner: _first_object_payload(
                owner.to_dataframe()
            ),
        ),
        _PublicFrameOwnerCase(
            name="ReferenceBundle.kinase_substrate_map",
            make_numeric_source=_reference_kinase_map,
            construct_from_numeric=lambda frame: _reference_bundle_from_kinase_map(
                frame
            ),
            observe_numeric=lambda owner: owner.kinase_substrate_map_dataframe(),
            make_object_source=lambda payload: _object_payload_frame(
                _reference_kinase_map(),
                payload,
            ),
            construct_from_object=lambda frame: _reference_bundle_from_kinase_map(
                frame
            ),
            observe_object_payload=lambda owner: _first_object_payload(
                owner.kinase_substrate_map_dataframe()
            ),
        ),
    )


def _json_immutability_cases() -> tuple[_JsonImmutabilityCase, ...]:
    return (
        _JsonImmutabilityCase(
            name="ResultCaveat.details",
            construct=lambda payload: _result_caveat_with_details(payload),
            observe=lambda owner: owner.details,
        ),
        _JsonImmutabilityCase(
            name="PhosphositeImportResult.diagnostics",
            construct=lambda payload: _import_result_from_frames(
                phospho=_small_phospho(),
                site_metadata=_sample_site_metadata(),
                diagnostics=payload,
            ),
            observe=lambda owner: owner.diagnostics,
        ),
        _JsonImmutabilityCase(
            name="EnrichmentWorkflowResult.diagnostics",
            construct=lambda payload: _enrichment_result_with_payload(
                "diagnostics",
                payload,
            ),
            observe=lambda owner: owner.diagnostics,
        ),
        _JsonImmutabilityCase(
            name="EnrichmentWorkflowResult.method_metadata",
            construct=lambda payload: _enrichment_result_with_payload(
                "method_metadata",
                payload,
            ),
            observe=lambda owner: owner.method_metadata,
        ),
        _JsonImmutabilityCase(
            name="EnrichmentWorkflowResult.background_summary",
            construct=lambda payload: _enrichment_result_with_payload(
                "background_summary",
                payload,
            ),
            observe=lambda owner: owner.background_summary,
        ),
        _JsonImmutabilityCase(
            name="EnrichmentWorkflowResult.set_collection_summary",
            construct=lambda payload: _enrichment_result_with_payload(
                "set_collection_summary",
                payload,
            ),
            observe=lambda owner: owner.set_collection_summary,
        ),
        _JsonImmutabilityCase(
            name="BatchCorrectionResult.diagnostics",
            construct=lambda payload: _batch_correction_result_with_diagnostics(
                payload
            ),
            observe=lambda owner: owner.diagnostics,
        ),
        _JsonImmutabilityCase(
            name="ProteinAwarePreparationReport.policy_parameters",
            construct=lambda payload: _protein_aware_report_with_policy_parameters(
                payload
            ),
            observe=lambda owner: owner.policy_parameters,
        ),
        _JsonImmutabilityCase(
            name="IntensityScaleEstablishmentProvenance.parameters",
            construct=lambda payload: _intensity_scale_provenance_with_parameters(
                payload
            ),
            observe=lambda owner: owner.parameters,
        ),
        _JsonImmutabilityCase(
            name="KinaseScoringResult.score_scale_metadata",
            construct=lambda payload: _kinase_scoring_result_with_metadata(payload),
            observe=lambda owner: owner.score_scale_metadata,
        ),
        _JsonImmutabilityCase(
            name="KinaseWorkflowAttritionProvenance.metrics",
            construct=lambda payload: _kinase_attrition_with_payload(
                "metrics",
                payload,
            ),
            observe=lambda owner: owner.metrics,
        ),
        _JsonImmutabilityCase(
            name="KinaseWorkflowAttritionProvenance.policy",
            construct=lambda payload: _kinase_attrition_with_payload(
                "policy",
                payload,
            ),
            observe=lambda owner: owner.policy,
        ),
        _JsonImmutabilityCase(
            name="KinaseWorkflowAttritionProvenance.policy_violations",
            construct=lambda payload: _kinase_attrition_with_payload(
                "policy_violations",
                payload,
            ),
            observe=lambda owner: owner.policy_violations[0],
        ),
    )


def _trusted_dataset(*, identity_details: Mapping[str, object]) -> Any:
    import pandas as pd

    from phospy.api import AnalysisReadyPhosphoDataset, Organism

    index = _site_index()
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
            index=index.copy(),
        ),
        site_metadata=_site_metadata(index, sites=("Y182", "S9")),
        organism=Organism.RAT,
        intensity_scale_state=_supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=_supported_linear_processing_state(has_total_matrix=False),
        trusted_construction_assertions=_trusted_assertions(
            identity_details=identity_details,
        ),
    )


def _small_phospho() -> Any:
    import pandas as pd

    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 1.0],
        },
        index=_site_index().copy(),
    )


def _sample_site_metadata() -> Any:
    return _site_metadata(_site_index(), sites=("Y182", "S9"))


def _analysis_ready_dataset_from_frames(*, phospho: Any, site_metadata: Any) -> object:
    from phospy.api import AnalysisReadyPhosphoDataset, Organism

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        return AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            intensity_scale_state=_supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=_supported_linear_processing_state(has_total_matrix=False),
        )


def _import_result_from_frames(
    *,
    phospho: Any,
    site_metadata: Any,
    diagnostics: Mapping[str, object] | None = None,
) -> object:
    from phospy.api import PhosphositeImportResult

    return PhosphositeImportResult(
        phospho_matrix_candidate=phospho,
        site_metadata_candidate=site_metadata,
        sample_column_mapping={"sample_a": "sample_a"},
        diagnostics=diagnostics,
    )


def _differential_result_table() -> Any:
    import numpy as np

    table = _sample_site_metadata()
    table.loc[:, "logFC"] = [1.0, -1.0]
    table.loc[:, "t"] = [2.0, -2.0]
    table.loc[:, "P.Value"] = [0.05, 0.10]
    table.loc[:, "adj.P.Val"] = [0.10, 0.10]
    table.loc[:, "result_status"] = ["tested", "tested"]
    table.loc[:, "result_status_reason"] = ["", ""]
    table.loc[:, "numeric_payload"] = np.asarray([1.0, 2.0])
    return table


def _differential_prior_diagnostics(index: Any) -> object:
    import numpy as np
    import pandas as pd

    from phospy.science.differential.models.diagnostics import (
        EmpiricalBayesPriorDiagnostics,
    )

    return EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=1.0,
        base_prior_degrees_of_freedom=10.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
    )


def _differential_result_from_table(table: Any) -> object:
    import numpy as np
    import pandas as pd

    from phospy.api.results import DifferentialAnalysisResult

    index = table.index.copy()
    return DifferentialAnalysisResult(
        residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="residual_variance",
        ),
        posterior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="posterior_residual_variance",
        ),
        prior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom_series_value=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
        prior_variance=1.0,
        prior_degrees_of_freedom=10.0,
        residual_degrees_of_freedom=4.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_differential_prior_diagnostics(index),
        mean_variance_trend_diagnostics=None,
        contrast_tables={"B_vs_A": table},
    )


def _kinase_score_matrix() -> Any:
    import pandas as pd

    return pd.DataFrame(
        {
            "MAP2K6": [0.8, 0.2],
            "AKT1": [0.2, 0.8],
        },
        index=_site_index().copy(),
    )


def _prediction_matrix() -> Any:
    import pandas as pd

    return pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=_site_index().copy(),
    )


def _activity_matrix() -> Any:
    import pandas as pd

    return pd.DataFrame(
        {"MAP2K6": [1.0, 2.0], "AKT1": [0.5, 1.5]},
        index=pd.Index(["sample_a", "sample_b"]),
    )


def _activity_target_table() -> Any:
    return _target_table_rows(score_values=(0.9, 0.8))


def _activity_result_from_matrix(matrix: Any) -> object:
    import pandas as pd

    from phospy.api.results import KinaseActivityResult

    return KinaseActivityResult(
        activity_matrix=matrix,
        thresholded_substrate_mean_activity=_activity_matrix(),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=_activity_target_table(),
    )


def _activity_result_from_target_table(table: Any) -> object:
    import pandas as pd

    from phospy.api.results import KinaseActivityResult

    return KinaseActivityResult(
        activity_matrix=_activity_matrix(),
        thresholded_substrate_mean_activity=_activity_matrix(),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=table,
    )


def _kinase_scoring_result_from_matrix(matrix: Any) -> object:
    from phospy.api.results import KinaseScoringResult

    return KinaseScoringResult(profile_scores=matrix)


def _kinase_scoring_result_with_metadata(payload: Mapping[str, object]) -> object:
    from phospy.api.results import KinaseScoringResult

    return KinaseScoringResult(
        profile_scores=_kinase_score_matrix(),
        score_scale_metadata=payload,
    )


def _kinase_prediction_result_from_matrix(matrix: Any) -> object:
    from phospy.api.results import KinasePredictionResult

    return KinasePredictionResult(pred_mat=matrix)


def _kinase_prediction_result_from_substrate_list(table: Any) -> object:
    from phospy.api.results import KinasePredictionResult

    return KinasePredictionResult(
        pred_mat=_prediction_matrix(),
        substrate_list=table,
    )


def _substrate_list() -> Any:
    return _target_table_rows(score_values=(0.9, 0.8), include_rank=True)


def _kinase_substrate_contribution_table() -> Any:
    import pandas as pd

    from phospy.tables.kinase import KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS

    return pd.DataFrame.from_records(
        [
            {
                "kinase": "MAP2K6",
                "substrate_site": "MAPK14;Y182;",
                "substrate_identifier": "MAPK14;Y182;",
                "value_used_in_scoring": 0.8,
                "score_component": "rank_weighted_fusion_scores",
                "score_source": "profile_only_motif_missing_or_constant",
                "reference_source_name": "fixture",
                "reference_source_version": "v1",
                "reference_bundle_id": "fixture_bundle",
                "reference_identifier_namespace": "display_id",
                "status": "included",
                "exclusion_reason": None,
                "ambiguous": False,
            }
        ],
        columns=pd.Index(KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS),
    )


def _kinase_workflow_result_from_contributions(table: Any) -> object:
    from phospy.api.results import (
        KinasePredictionResult,
        KinaseScoringResult,
        KinaseWorkflowResult,
    )

    return KinaseWorkflowResult(
        dataset=_analysis_ready_dataset_from_frames(
            phospho=_small_phospho(),
            site_metadata=_sample_site_metadata(),
        ),
        references=_reference_bundle_from_kinase_map(_reference_kinase_map()),
        scoring_result=KinaseScoringResult(profile_scores=_kinase_score_matrix()),
        prediction_result=KinasePredictionResult(pred_mat=_prediction_matrix()),
        substrate_contributions=table,
    )


def _signalome_expanded_table() -> Any:
    import pandas as pd

    from phospy.science.signalomes.constants import (
        DISPLAY_ID_COLUMN,
        SITE_KEY_COLUMN,
    )

    return pd.DataFrame(
        {
            SITE_KEY_COLUMN: _site_index().astype(str).tolist(),
            DISPLAY_ID_COLUMN: ["MAPK14;Y182;", "GSK3B;S9;"],
            "numeric_payload": [1.0, 2.0],
        },
        index=_site_index().copy(),
    )


def _signalome_result_from_expanded_table(table: Any) -> object:
    from phospy.api.results import SignalomeWorkflowResult
    from phospy.science.signalomes.models import (
        KinaseNetwork,
        SignalomeAssignments,
        SignalomeModules,
    )

    return SignalomeWorkflowResult(
        dataset=_analysis_ready_dataset_from_frames(
            phospho=_small_phospho(),
            site_metadata=_sample_site_metadata(),
        ),
        kinase_result=_kinase_workflow_result_from_contributions(
            _kinase_substrate_contribution_table()
        ),
        module_assignments=SignalomeAssignments._from_owned(
            table=_empty_signalome_assignments_table()
        ),
        signalome_modules=SignalomeModules._from_owned(
            table=_empty_signalome_modules_table()
        ),
        kinase_network=KinaseNetwork._from_owned(edges=_empty_kinase_network_edges()),
        expanded_signalome=table,
    )


def _empty_signalome_assignments_table() -> Any:
    import pandas as pd

    from phospy.science.signalomes.constants import (
        DISPLAY_ID_COLUMN,
        GENE_SYMBOL_COLUMN,
        ISOFORM_ID_COLUMN,
        MODULE_ID_COLUMN,
        MODULE_TOP_KINASE_CANDIDATES_COLUMN,
        MODULE_TOP_KINASE_COLUMN,
        MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
        MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
        MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
        PROTEIN_ACCESSION_COLUMN,
        PROTEIN_COLUMN,
        SITE_COLUMN,
        SITE_KEY_COLUMN,
        TOP_KINASE_CANDIDATES_COLUMN,
        TOP_KINASE_COLUMN,
        TOP_KINASE_IS_AMBIGUOUS_COLUMN,
        TOP_KINASE_SELECTION_POLICY_COLUMN,
        TOP_KINASE_TIE_COUNT_COLUMN,
        TOP_KINASE_WEIGHTS_COLUMN,
        TOP_SCORE_COLUMN,
    )

    columns = (
        SITE_KEY_COLUMN,
        DISPLAY_ID_COLUMN,
        GENE_SYMBOL_COLUMN,
        SITE_COLUMN,
        PROTEIN_COLUMN,
        PROTEIN_ACCESSION_COLUMN,
        ISOFORM_ID_COLUMN,
        MODULE_ID_COLUMN,
        TOP_KINASE_COLUMN,
        TOP_SCORE_COLUMN,
        TOP_KINASE_CANDIDATES_COLUMN,
        TOP_KINASE_WEIGHTS_COLUMN,
        TOP_KINASE_TIE_COUNT_COLUMN,
        TOP_KINASE_IS_AMBIGUOUS_COLUMN,
        TOP_KINASE_SELECTION_POLICY_COLUMN,
        MODULE_TOP_KINASE_COLUMN,
        MODULE_TOP_KINASE_CANDIDATES_COLUMN,
        MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
        MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
        MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
    )
    return pd.DataFrame(columns=columns, index=_site_index()[:0].copy())


def _empty_signalome_modules_table() -> Any:
    import pandas as pd

    return pd.DataFrame({"MAP2K6": pd.Series(dtype="float64")})


def _empty_kinase_network_edges() -> Any:
    import pandas as pd

    from phospy.science.signalomes.constants import (
        CORRELATION_COLUMN,
        SOURCE_KINASE_COLUMN,
        TARGET_KINASE_COLUMN,
    )

    return pd.DataFrame(
        columns=(SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN, CORRELATION_COLUMN)
    )


def _reference_kinase_map() -> Any:
    import pandas as pd

    return pd.DataFrame(
        {
            "kinase": ["MAP2K6", "AKT1"],
            "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            "numeric_payload": [1.0, 2.0],
        }
    )


def _reference_bundle_from_kinase_map(kinase_substrate_map: Any) -> object:
    from phospy.api import Organism, ReferenceBundle

    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=kinase_substrate_map,
        site_sequences=_reference_site_sequences(),
    )


def _reference_site_sequences() -> Any:
    import pandas as pd

    return pd.DataFrame(
        {
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ]
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )


def _target_table_rows(
    *,
    score_values: tuple[float, float],
    include_rank: bool = False,
) -> Any:
    import pandas as pd

    payload: dict[str, object] = {
        "kinase": ["MAP2K6", "AKT1"],
        "site_id": ["MAPK14;Y182;", "GSK3B;S9;"],
        "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
        "site_key": _site_index().astype(str).tolist(),
        "display_id": ["MAPK14;Y182;", "GSK3B;S9;"],
        "score": list(score_values),
    }
    if include_rank:
        payload["rank"] = [1, 1]
    return pd.DataFrame(payload)


def _result_caveat_with_details(payload: Mapping[str, object]) -> object:
    from phospy.api import ResultCaveat

    return ResultCaveat(
        code="public_boundary_json",
        severity="info",
        message="JSON state is protected.",
        details=payload,
    )


def _enrichment_result_with_payload(
    field_name: str,
    payload: Mapping[str, object],
) -> object:
    from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    from phospy.contracts.configs import EnrichmentConfig
    from phospy.contracts.results import EnrichmentWorkflowResult
    from phospy.science.enrichment.models import GeneSetCollection

    values: dict[str, object] = {
        "diagnostics": {},
        "method_metadata": {},
        "background_summary": {},
        "set_collection_summary": {},
    }
    values[field_name] = payload
    return EnrichmentWorkflowResult(
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={"mapk_pathway": ("AKT1", "MAPK1")},
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        ),
        config=EnrichmentConfig(),
        **values,
    )


def _batch_correction_result_with_diagnostics(
    payload: Mapping[str, object],
) -> object:
    import pandas as pd

    from phospy.science.datasets.preprocessing.batch_correction import (
        BATCH_CORRECTION_STATUS_APPLIED,
        BatchCorrectionDiagnostics,
        BatchCorrectionPolicy,
        BatchCorrectionReport,
        BatchCorrectionResult,
    )

    report = BatchCorrectionReport(
        status=BATCH_CORRECTION_STATUS_APPLIED,
        policy=BatchCorrectionPolicy(
            method="linear_residualize_batch",
            batch_column="batch",
            condition_column="condition",
        ),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=2,
            batch_levels=("b1", "b2"),
            condition_levels=("A", "B"),
            matrix_shape_before=(1, 1),
            matrix_shape_after=(1, 1),
        ),
    )
    return BatchCorrectionResult(
        corrected_matrix=pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["site_a"], name="site_id"),
        ),
        report=report,
        diagnostics=payload,
    )


def _protein_aware_report_with_policy_parameters(
    payload: Mapping[str, object],
) -> object:
    from phospy.science.datasets.preprocessing.protein_aware_alignment import (
        ProteinAwarePreparationEligibility,
        ProteinAwareSampleAlignmentDiagnostics,
    )
    from phospy.science.datasets.preprocessing.protein_aware_preparation import (
        ProteinAwarePreparationReport,
        ProteinAwareSiteEligibility,
    )
    from phospy.science.datasets.preprocessing.protein_mapping import (
        ProteinMappingStatus,
    )

    return ProteinAwarePreparationReport(
        site_eligibility=(
            ProteinAwareSiteEligibility(
                site_key="site_a",
                eligibility=(
                    ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY
                ),
                mapping_status=ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER,
                reasons=("compact verifier fixture",),
            ),
        ),
        sample_alignment=ProteinAwareSampleAlignmentDiagnostics(
            phospho_sample_columns=("sample_a",),
            total_protein_sample_columns=("sample_a",),
            exact_sample_order_match=True,
            sample_order_compatible=True,
            reordered_sample_columns=False,
            allow_reordered_samples=False,
            missing_total_protein_samples=(),
            extra_total_protein_samples=(),
        ),
        policy_parameters=payload,
    )


def _intensity_scale_provenance_with_parameters(
    payload: Mapping[str, object],
) -> object:
    from phospy.science.transformations.models import (
        IntensityScaleEstablishmentMode,
        IntensityScaleEstablishmentProvenance,
        IntensityScaleEstablishmentSource,
        IntensityScaleEvidenceLevel,
    )

    return IntensityScaleEstablishmentProvenance(
        scale="linear",
        mode=IntensityScaleEstablishmentMode.DECLARED,
        source=IntensityScaleEstablishmentSource.DECLARED_BY_USER,
        evidence_level=IntensityScaleEvidenceLevel.DECLARED_BY_USER,
        parameters=payload,
    )


def _kinase_attrition_with_payload(
    field_name: str,
    payload: Mapping[str, object],
) -> object:
    from phospy.api.results import KinaseWorkflowAttritionProvenance

    values: dict[str, object] = {
        "metrics": {"retained": 2},
        "policy": {"minimum_scored_fraction": 0.5},
        "policy_violations": (),
    }
    if field_name == "policy_violations":
        values[field_name] = (payload,)
    else:
        values[field_name] = payload
    return KinaseWorkflowAttritionProvenance(
        metrics=values["metrics"],
        policy=values["policy"],
        policy_outcome="warned",
        policy_violations=values["policy_violations"],
    )


def _expect_direct_constructor_rejects_stale_provenance(dataset: Any) -> str:
    return _expect_public_dataset_constructor_rejects_stale_provenance(dataset)


def _expect_public_dataset_constructor_rejects_stale_provenance(
    dataset: Any,
    *,
    constructor: Callable[..., object] | None = None,
) -> str:
    from phospy.api import AnalysisReadyPhosphoDataset
    from phospy.errors import DatasetValidationError

    phospho = dataset.phospho
    phospho.iloc[0, 0] = float(phospho.iloc[0, 0]) + 10.0

    def _construct_with_stale_provenance() -> Any:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            dataset_constructor = constructor or AnalysisReadyPhosphoDataset
            return dataset_constructor(
                phospho=phospho,
                site_metadata=dataset.site_metadata,
                sample_metadata=dataset.sample_metadata,
                total=dataset.total,
                comparisons=dataset.comparisons,
                imputation_observation_mask=(
                    dataset.imputation_observed_mask_dataframe()
                ),
                organism=dataset.organism,
                intensity_scale_state=dataset.intensity_scale_state,
                processing_state=dataset.processing_state,
                provenance=dataset.provenance,
                trusted_construction_assertions=(
                    dataset.trusted_construction_assertions
                ),
            )

    message = _expect_raises(
        DatasetValidationError,
        _construct_with_stale_provenance,
        contains="dataset.phospho",
    )
    _require(
        "expected exact digest" in message and "actual exact digest" in message,
        "stale direct-constructor provenance error did not report expected and "
        "actual table digests",
    )
    return message


def _trusted_assertions(*, identity_details: Mapping[str, object]) -> Any:
    from phospy.provenance import (
        TrustedDatasetConstructionAssertions,
        TrustedDatasetConstructionEvidence,
    )

    return TrustedDatasetConstructionAssertions(
        identity=TrustedDatasetConstructionEvidence.evidence(
            source="protein-scoped site_key export",
            details=identity_details,
        ),
        intensity_scale=TrustedDatasetConstructionEvidence.evidence(
            source="installed artifact fixture",
            policy="require_established_intensity_scale_state",
            details={"scale": "linear"},
        ),
        quantitative_meaning=TrustedDatasetConstructionEvidence.evidence(
            source="installed artifact fixture",
            policy="analysis-ready quantitative matrix",
        ),
        aligned_structure=TrustedDatasetConstructionEvidence.evidence(
            source="pre-export table alignment audit",
            policy="require_identical_site_indexes_and_sample_axes",
        ),
        localisation=TrustedDatasetConstructionEvidence.evidence(
            source="localisation_confidence column",
            policy="require_threshold",
            threshold=0.75,
        ),
        sequence=TrustedDatasetConstructionEvidence.evidence(
            source="site_sequence column",
        ),
        reference_context=TrustedDatasetConstructionEvidence.waiver(
            reason="compact installed-artifact fixture has no source manifest",
        ),
        asserted_by="distribution-artifact-verifier",
        assertion_source="installed artifact corrected-path verifier",
    )


def _differential_dataset() -> Any:
    import pandas as pd

    from phospy.api import AnalysisReadyPhosphoDataset, Organism

    sites = ("Y182", "S9", "T308")
    index = _site_index(
        protein_identifiers=("MAPK14", "GSK3B", "AKT1"),
        sites=sites,
    )
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0, 1.0],
                "A_2": [1.1, 2.1, 1.1],
                "B_1": [2.1, 2.0, 1.0],
                "B_2": [2.0, 2.2, 0.9],
            },
            index=index.copy(),
        ),
        site_metadata=_site_metadata(index, sites=sites),
        organism=Organism.RAT,
        intensity_scale_state=_supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=_supported_log2_processing_state(has_total_matrix=False),
        trusted_construction_assertions=_trusted_assertions(
            identity_details={"columns": ["site_key"], "fixture": "differential"},
        ),
    )


def _derived_fixture() -> _DerivedFixture:
    import pandas as pd

    from phospy.api import AnalysisReadyPhosphoDataset, Organism
    from phospy.provenance import (
        DerivedQuantitativeDataProvenance,
        DerivedSampleMapping,
    )
    from phospy.science.datasets.derived_quantitative import (
        CertifiedDerivedQuantitativeParentState,
    )

    sites = ("Y182", "S9")
    index = _site_index(sites=sites)
    site_metadata = _site_metadata(index, sites=sites)
    parent_dataset = AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=pd.DataFrame(
            {
                "bio_a_t1": [1.0, 3.0],
                "bio_a_t2": [3.0, 5.0],
                "bio_b_t1": [5.0, 7.0],
                "bio_b_t2": [7.0, 9.0],
            },
            index=index.copy(),
        ),
        site_metadata=site_metadata.copy(deep=True),
        organism=Organism.RAT,
        intensity_scale_state=_supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=_supported_linear_processing_state(has_total_matrix=False),
        trusted_construction_assertions=_trusted_assertions(
            identity_details={"columns": ["site_key"], "fixture": "derived_parent"},
        ),
    )
    parent_phospho = parent_dataset.phospho
    derived_phospho = pd.DataFrame(
        {
            "bio_a": parent_phospho.loc[:, ["bio_a_t1", "bio_a_t2"]].mean(axis=1),
            "bio_b": parent_phospho.loc[:, ["bio_b_t1", "bio_b_t2"]].mean(axis=1),
        },
        index=index.copy(),
    )
    parent_state = CertifiedDerivedQuantitativeParentState.from_dataset(parent_dataset)
    derived_fingerprints = _fingerprints_for_tables(
        phospho=derived_phospho,
        site_metadata=site_metadata,
    )
    lineage = DerivedQuantitativeDataProvenance(
        derivation_type="technical_replicate_aggregation",
        parent_dataset_type=parent_state.parent_dataset_type,
        derived_dataset_type="DerivedAnalysisReadyPhosphoDataset",
        parent_dataset_fingerprints=parent_state.parent_dataset_fingerprints,
        derived_dataset_fingerprints=derived_fingerprints,
        sample_mapping=(
            DerivedSampleMapping(
                output_sample_id="bio_a",
                input_sample_ids=("bio_a_t1", "bio_a_t2"),
                condition="A",
                biological_replicate_id="bio_a",
                technical_replicate_ids=("t1", "t2"),
            ),
            DerivedSampleMapping(
                output_sample_id="bio_b",
                input_sample_ids=("bio_b_t1", "bio_b_t2"),
                condition="B",
                biological_replicate_id="bio_b",
                technical_replicate_ids=("t1", "t2"),
            ),
        ),
        aggregation_method="mean",
        input_intensity_scale="linear",
        output_intensity_scale="linear",
        quantitative_meaning="phosphosite_abundance",
        missingness_policy={"policy": "complete_matrix"},
        matrices_transformed={
            "phospho": True,
            "sample_metadata": False,
            "total_protein": False,
            "imputation_observation_mask": False,
            "comparisons": False,
        },
        implementation="distribution-artifact-verifier",
        implementation_version="1",
    )
    return _DerivedFixture(
        phospho=derived_phospho,
        site_metadata=site_metadata,
        parent_state=parent_state,
        lineage=lineage,
    )


def _fingerprints_for_tables(*, phospho: Any, site_metadata: Any) -> tuple[Any, ...]:
    from phospy.provenance import fingerprint_optional_table_strict

    fingerprints = [
        fingerprint_optional_table_strict(phospho, name="dataset.phospho"),
        fingerprint_optional_table_strict(site_metadata, name="dataset.site_metadata"),
    ]
    return tuple(item for item in fingerprints if item is not None)


def _supported_linear_intensity_scale_state(*, has_total_matrix: bool) -> Any:
    return _supported_intensity_scale_state(
        scale="linear",
        has_total_matrix=has_total_matrix,
    )


def _supported_log2_intensity_scale_state(*, has_total_matrix: bool) -> Any:
    return _supported_intensity_scale_state(
        scale="log2",
        has_total_matrix=has_total_matrix,
    )


def _supported_linear_processing_state(*, has_total_matrix: bool) -> Any:
    return _supported_processing_state(
        scale="linear",
        has_total_matrix=has_total_matrix,
    )


def _supported_log2_processing_state(*, has_total_matrix: bool) -> Any:
    return _supported_processing_state(
        scale="log2",
        has_total_matrix=has_total_matrix,
    )


def _supported_processing_state(*, scale: str, has_total_matrix: bool) -> Any:
    from phospy.science.datasets.builders.preprocessing import (
        build_dataset_processing_state,
    )
    from phospy.science.datasets.preprocessing.models import PreprocessingPlan

    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=_supported_intensity_scale_state(
            scale=scale,
            has_total_matrix=has_total_matrix,
        ),
    )


def _supported_intensity_scale_state(*, scale: str, has_total_matrix: bool) -> Any:
    import pandas as pd

    from phospy.science.datasets.builders.transformation_resolver import (
        DatasetIntensityScaleResolver,
    )
    from phospy.science.transformations.models import (
        IntensityScaleKind,
        IntensityScaleState,
        MatrixIntensityScaleState,
    )
    from phospy.science.transformations.transformers import IdentityTransformer

    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["GENEA;S1;"], name="site_id"),
    )
    total = None
    if has_total_matrix:
        total = pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["GENEA"], name="protein_id"),
        )
    if scale == "linear":
        expected_kind = IntensityScaleKind.LINEAR
        phospho_state = MatrixIntensityScaleState.linear(established_by="trusted.input")
        total_state = (
            MatrixIntensityScaleState.linear(established_by="trusted.input")
            if has_total_matrix
            else None
        )
    elif scale == "log2":
        expected_kind = IntensityScaleKind.LOG2
        phospho_state = MatrixIntensityScaleState.log2(established_by="trusted.input")
        total_state = (
            MatrixIntensityScaleState.log2(established_by="trusted.input")
            if has_total_matrix
            else None
        )
    else:
        raise VerificationError(f"unsupported verifier intensity scale: {scale}")

    declared_state = IntensityScaleState(phospho=phospho_state, total=total_state)
    return (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=total,
            expected_scale_kind=expected_kind,
            declared_input_scale_state=declared_state,
            input_declaration_source="scripts.verify_distribution_artifact",
        )
        .intensity_scale_state
    )


def _site_index(
    *,
    protein_identifiers: tuple[str, ...] = ("MAPK14", "GSK3B"),
    sites: tuple[str, ...] = ("Y182", "S9"),
) -> Any:
    import pandas as pd

    return pd.Index(
        [
            _protein_site_key(protein_identifier=protein_identifier, site=site)
            for protein_identifier, site in zip(
                protein_identifiers,
                sites,
                strict=True,
            )
        ],
        name="site_key",
    )


def _site_metadata(index: Any, *, sites: tuple[str, ...]) -> Any:
    import pandas as pd

    genes = tuple(str(key.protein_identifier) for key in _decoded_index(index))
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(genes, sites, strict=True)
            ],
            **_site_key_context_columns(index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
            "localisation_confidence": [0.95] * len(index),
        },
        index=index.copy(),
    )


def _protein_site_key(*, protein_identifier: str, site: str) -> str:
    from phospy.api import Organism
    from phospy.science.sites.site_keys import (
        ProteinScopedPhosphositeKey,
        encode_site_key,
    )

    return encode_site_key(
        ProteinScopedPhosphositeKey(
            organism=Organism.RAT,
            protein_namespace="protein_id",
            protein_identifier=protein_identifier,
            residue=site[0],
            position=int(site[1:]),
        )
    )


def _site_key_context_columns(site_keys: Any) -> dict[str, list[str]]:
    decoded_keys = _decoded_index(site_keys)
    return {
        "organism": [key.organism.value for key in decoded_keys],
        "protein_namespace": [key.protein_namespace for key in decoded_keys],
        "protein_identifier": [key.protein_identifier for key in decoded_keys],
    }


def _decoded_index(index: Any) -> tuple[Any, ...]:
    from phospy.science.sites.site_keys import decode_site_key

    return tuple(
        decode_site_key(
            value,
            field_name="distribution_artifact_verifier.site_key",
            error_type=ValueError,
        )
        for value in index.astype(str).tolist()
    )


def _assert_numeric_frame_owner_case_isolated(case: _PublicFrameOwnerCase) -> None:
    source = case.make_numeric_source()
    owner = case.construct_from_numeric(source)
    before = _require_frame(case.observe_numeric(owner), field_name=case.name).copy(
        deep=True
    )

    _mutate_first_numeric_cell(source, 999.0)
    _assert_frame_equal(
        _require_frame(case.observe_numeric(owner), field_name=case.name),
        before,
        label=f"{case.name} retained caller-owned DataFrame input",
    )

    exported = _require_frame(case.observe_numeric(owner), field_name=case.name)
    _mutate_first_numeric_cell(exported, 777.0)
    _assert_frame_equal(
        _require_frame(case.observe_numeric(owner), field_name=case.name),
        before,
        label=f"{case.name} returned a mutable public DataFrame alias",
    )


def _assert_object_frame_owner_case_isolated(case: _PublicFrameOwnerCase) -> None:
    _require(
        case.make_object_source is not None
        and case.construct_from_object is not None
        and case.observe_object_payload is not None,
        f"{case.name} object-cell ownership probe is incomplete",
    )
    payload = _mutable_object_payload()
    source = case.make_object_source(payload)
    owner = case.construct_from_object(source)

    _mutate_object_payload(payload, "caller")
    observed_payload = case.observe_object_payload(owner)
    _require(
        _object_payload_state(observed_payload) == _PUBLIC_FRAME_OBJECT_PAYLOAD_STATE,
        f"{case.name} retained caller-owned nested object-dtype cell",
    )

    exported_payload = case.observe_object_payload(owner)
    _mutate_object_payload(exported_payload, "export")
    _require(
        _object_payload_state(case.observe_object_payload(owner))
        == _PUBLIC_FRAME_OBJECT_PAYLOAD_STATE,
        f"{case.name} returned a mutable public nested object-dtype cell alias",
    )


def _assert_json_owner_case_isolated(case: _JsonImmutabilityCase) -> None:
    source = _mutable_json_payload()
    owner = case.construct(source)
    before = _json_state(case.observe(owner))

    _mutate_json_payload(source, "caller")
    _require(
        _json_state(case.observe(owner)) == before,
        f"{case.name} retained caller-owned JSON mapping state",
    )

    exported_payload = _detached_json_payload(case.observe(owner))
    _mutate_json_payload(exported_payload, "export")
    _require(
        _json_state(case.observe(owner)) == before,
        f"{case.name} returned nested mutable JSON state",
    )


def _require_frame(value: object, *, field_name: str) -> Any:
    import pandas as pd

    _require(
        isinstance(value, pd.DataFrame), f"{field_name} did not expose a DataFrame"
    )
    return value


def _assert_frame_equal(left: Any, right: Any, *, label: str) -> None:
    import pandas as pd

    try:
        pd.testing.assert_frame_equal(left, right)
    except AssertionError as exc:
        raise VerificationError(f"{label}: {exc}") from exc


def _first_numeric_column(frame: Any) -> str:
    import pandas as pd

    for column in frame.columns:
        if pd.api.types.is_numeric_dtype(frame.loc[:, column]):
            return str(column)
    raise VerificationError("public DataFrame ownership probe has no numeric column")


def _mutate_first_numeric_cell(frame: Any, value: float) -> None:
    column = _first_numeric_column(frame)
    frame.loc[frame.index[0], column] = value


def _object_payload_frame(frame: Any, payload: object) -> Any:
    copied = frame.copy(deep=True)
    copied.loc[:, _PUBLIC_FRAME_OBJECT_PAYLOAD_COLUMN] = [
        payload,
        *[_mutable_object_payload() for _ in range(len(copied.index) - 1)],
    ]
    return copied


def _first_object_payload(frame: Any) -> object:
    _require(
        _PUBLIC_FRAME_OBJECT_PAYLOAD_COLUMN in frame.columns,
        "public DataFrame ownership probe lost the nested object payload column",
    )
    return frame.loc[frame.index[0], _PUBLIC_FRAME_OBJECT_PAYLOAD_COLUMN]


def _mutable_object_payload() -> dict[str, object]:
    import numpy as np

    return {
        "list": ["list-start"],
        "dict": {"inner": ["dict-start"]},
        "array": np.asarray([1.0, 2.0]),
        "set": {"set-start"},
        "nested": [
            {"array": np.asarray([3.0, 4.0])},
            {"set": {"nested-set-start"}},
            ["nested-list-start"],
        ],
    }


def _mutate_object_payload(payload: object, marker: str) -> None:
    import numpy as np

    _require(isinstance(payload, dict), "object payload must be a dict")
    list_value = payload["list"]
    dict_value = payload["dict"]
    array_value = payload["array"]
    set_value = payload["set"]
    nested_value = payload["nested"]
    _require(isinstance(list_value, list), "object payload list is invalid")
    _require(isinstance(dict_value, dict), "object payload dict is invalid")
    _require(isinstance(array_value, np.ndarray), "object payload array is invalid")
    _require(isinstance(set_value, set), "object payload set is invalid")
    _require(isinstance(nested_value, list), "object payload nested list is invalid")
    nested_array_mapping = nested_value[0]
    nested_set_mapping = nested_value[1]
    nested_list = nested_value[2]
    _require(isinstance(nested_array_mapping, dict), "nested array mapping is invalid")
    _require(isinstance(nested_set_mapping, dict), "nested set mapping is invalid")
    _require(isinstance(nested_list, list), "nested list is invalid")
    nested_array = nested_array_mapping["array"]
    nested_set = nested_set_mapping["set"]
    _require(isinstance(nested_array, np.ndarray), "nested array is invalid")
    _require(isinstance(nested_set, set), "nested set is invalid")

    list_value.append(f"{marker}-list")
    dict_inner = dict_value["inner"]
    _require(isinstance(dict_inner, list), "nested dict inner list is invalid")
    dict_inner.append(f"{marker}-dict")
    array_value[0] = 99.0
    set_value.add(f"{marker}-set")
    nested_array[0] = 88.0
    nested_set.add(f"{marker}-nested-set")
    nested_list.append(f"{marker}-nested-list")


def _object_payload_state(payload: object) -> dict[str, tuple[object, ...]]:
    import numpy as np

    _require(isinstance(payload, dict), "object payload must be a dict")
    nested_value = payload["nested"]
    _require(isinstance(nested_value, list), "object payload nested list is invalid")
    nested_array_mapping = nested_value[0]
    nested_set_mapping = nested_value[1]
    nested_list = nested_value[2]
    _require(isinstance(nested_array_mapping, dict), "nested array mapping is invalid")
    _require(isinstance(nested_set_mapping, dict), "nested set mapping is invalid")
    _require(isinstance(nested_list, list), "nested list is invalid")
    array_value = payload["array"]
    nested_array = nested_array_mapping["array"]
    set_value = payload["set"]
    nested_set = nested_set_mapping["set"]
    dict_value = payload["dict"]
    list_value = payload["list"]
    _require(isinstance(array_value, np.ndarray), "object payload array is invalid")
    _require(isinstance(nested_array, np.ndarray), "nested array is invalid")
    _require(isinstance(set_value, set), "object payload set is invalid")
    _require(isinstance(nested_set, set), "nested set is invalid")
    _require(isinstance(dict_value, dict), "object payload dict is invalid")
    _require(isinstance(list_value, list), "object payload list is invalid")
    dict_inner = dict_value["inner"]
    _require(isinstance(dict_inner, list), "nested dict inner list is invalid")
    return {
        "list": tuple(list_value),
        "dict": tuple(dict_inner),
        "array": tuple(float(value) for value in array_value.tolist()),
        "set": tuple(sorted(str(value) for value in set_value)),
        "nested_array": tuple(float(value) for value in nested_array.tolist()),
        "nested_set": tuple(sorted(str(value) for value in nested_set)),
        "nested_list": tuple(nested_list),
    }


def _mutable_json_payload() -> dict[str, object]:
    return {
        "nested": {
            "items": [1],
            "metadata": {"markers": ["start"]},
        },
        "rows": [{"values": [1, 2]}],
    }


def _mutate_json_payload(payload: object, marker: str) -> None:
    _require(isinstance(payload, dict), "JSON payload must be a mutable dict")
    nested = payload["nested"]
    rows = payload["rows"]
    _require(isinstance(nested, dict), "JSON payload nested value must be a dict")
    _require(isinstance(rows, list), "JSON payload rows value must be a list")
    items = nested["items"]
    metadata = nested["metadata"]
    _require(isinstance(items, list), "JSON payload nested items must be a list")
    _require(isinstance(metadata, dict), "JSON payload metadata must be a dict")
    markers = metadata["markers"]
    _require(isinstance(markers, list), "JSON payload metadata markers must be a list")
    row = rows[0]
    _require(isinstance(row, dict), "JSON payload row must be a dict")
    values = row["values"]
    _require(isinstance(values, list), "JSON payload row values must be a list")
    items.append(marker)
    markers.append(marker)
    values.append(marker)


def _detached_json_payload(value: object) -> dict[str, object]:
    copier = getattr(value, "copy", None)
    if callable(copier):
        payload = copier()
    else:
        payload = copy.deepcopy(value)
    _require(isinstance(payload, dict), "JSON payload export must be a dict")
    return payload


def _json_state(value: object) -> object:
    if isinstance(value, Mapping):
        return tuple(
            (str(key), _json_state(item))
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return tuple(_json_state(item) for item in value)
    return value


def _try_mutate_public_provenance_details(dataset: Any) -> None:
    try:
        public_parameters = dataset.provenance.workflow_parameters
        construction = _mapping(public_parameters["construction"])
        assertion_payload = _mapping(construction["trusted_construction_assertions"])
        identity_payload = _mapping(assertion_payload["identity"])
        mutable_public_details = _mapping(identity_payload["details"])
        columns = mutable_public_details.get("columns")
        if isinstance(columns, list):
            columns.append("public-only")
        mutable_public_details["schema"] = {"version": 42}
    except (AttributeError, KeyError, TypeError):
        return


def _expect_raises(
    expected_type: type[BaseException],
    callback: Callable[[], object],
    *,
    contains: str | None = None,
) -> str:
    try:
        callback()
    except expected_type as exc:
        message = str(exc)
        if contains is not None and contains not in message:
            raise VerificationError(
                f"expected {expected_type.__name__} message to contain "
                f"{contains!r}; got {message!r}"
            ) from exc
        return message
    except Exception as exc:
        raise VerificationError(
            f"expected {expected_type.__name__}, got {type(exc).__name__}: {exc}"
        ) from exc
    raise VerificationError(f"expected {expected_type.__name__} to be raised")


def _module_file(module: object, *, module_name: str) -> Path:
    module_file = getattr(module, "__file__", None)
    _require(
        isinstance(module_file, str) and bool(module_file.strip()),
        f"{module_name} has no resolvable __file__",
    )
    return Path(module_file).resolve()


def _distribution_init_paths(distribution: metadata.Distribution) -> tuple[Path, ...]:
    files = distribution.files
    _require(
        files is not None,
        "importlib.metadata returned no file inventory for phospy",
    )
    init_paths = tuple(
        Path(distribution.locate_file(file_path)).resolve()
        for file_path in files
        if file_path.as_posix().replace("\\", "/") == "phospy/__init__.py"
    )
    _require(
        init_paths != (),
        "importlib.metadata file inventory does not include phospy/__init__.py",
    )
    return init_paths


def _assert_loaded_phospy_modules_not_from_repository(
    context: VerificationContext,
) -> None:
    if context.package_root is None:
        return
    loaded_paths: dict[str, str] = {}
    for module_name, module in sorted(sys.modules.items()):
        if module_name != "phospy" and not module_name.startswith("phospy."):
            continue
        module_file = getattr(module, "__file__", None)
        if not isinstance(module_file, str) or not module_file.strip():
            continue
        module_path = Path(module_file).resolve()
        _reject_if_inside_repository(
            module_path,
            context.repository_root,
            description=f"loaded module {module_name}",
        )
        if not _is_relative_to(module_path, context.package_root):
            raise VerificationError(
                f"loaded module {module_name} did not originate from the installed "
                f"package root {context.package_root.as_posix()}: "
                f"{module_path.as_posix()}"
            )
        loaded_paths[module_name] = module_path.as_posix()
    _require(
        "phospy" in loaded_paths,
        "no loaded phospy package module was found after identity verification",
    )


def _reject_if_inside_repository(
    path: Path,
    repository_root: Path,
    *,
    description: str,
) -> None:
    if _is_relative_to(path, repository_root):
        raise VerificationError(
            f"{description} resolved inside the source checkout: {path.as_posix()}"
        )


def _summarize_check_statuses(
    check_results: Sequence[CheckResult],
) -> dict[str, str]:
    statuses: dict[str, str] = {
        key: "skipped" for keys in SUMMARY_CHECK_KEYS_BY_DETAIL.values() for key in keys
    }
    for result in check_results:
        for key in SUMMARY_CHECK_KEYS_BY_DETAIL.get(result.name, ()):
            statuses[key] = result.status
    return statuses


def _dependency_snapshot() -> dict[str, str]:
    snapshot: dict[str, str] = {}
    for distribution in metadata.distributions():
        try:
            name = distribution.metadata.get("Name")
            version = distribution.version
        except Exception:
            continue
        if isinstance(name, str) and name.strip() and isinstance(version, str):
            snapshot[name.strip().lower()] = version
    return dict(sorted(snapshot.items()))


def _payload_sha256(payload: Mapping[str, object]) -> str:
    serialized = json.dumps(
        _json_safe(payload),
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(serialized).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json_mapping(path: Path) -> dict[str, object]:
    try:
        with path.open("rb") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise VerificationError(
            f"build manifest is not valid JSON: {path.as_posix()}: {exc}"
        ) from exc
    _require(isinstance(payload, Mapping), "build manifest must be a JSON object")
    return dict(payload)


def _require_artifact_filename_matches_kind(path: Path, artifact_kind: str) -> None:
    if artifact_kind == "wheel":
        _require(
            path.name.endswith(".whl"),
            f"wheel artifact path must end in .whl: {path.as_posix()}",
        )
        return
    if artifact_kind == "sdist":
        _require(
            path.name.endswith(".tar.gz"),
            f"sdist artifact path must end in .tar.gz: {path.as_posix()}",
        )
        return
    raise VerificationError(f"unsupported artifact kind: {artifact_kind!r}")


def _required_text(value: object, *, field_name: str) -> str:
    _require(
        isinstance(value, str) and bool(value.strip()),
        f"{field_name} must be non-empty text",
    )
    return value.strip()


def _require_text_equals(
    value: object,
    expected: str,
    *,
    field_name: str,
) -> None:
    observed = _required_text(value, field_name=field_name)
    _require(
        observed == expected,
        f"{field_name} mismatch: expected {expected!r}, got {observed!r}",
    )


def _require_sha256_digest(value: str, *, field_name: str) -> None:
    _require(
        value.startswith("sha256:"),
        f"{field_name} must use sha256: prefix",
    )
    _require_sha256_hex(value.removeprefix("sha256:"), field_name=field_name)


def _require_sha256_hex(value: str, *, field_name: str) -> None:
    _require(
        len(value) == 64
        and all(character in "0123456789abcdef" for character in value),
        f"{field_name} must be 64 lowercase hexadecimal characters",
    )


def _require_isolated_interpreter() -> None:
    _require(
        bool(sys.flags.isolated),
        "artifact verifier must run with isolated Python: use python -I",
    )


def _require_external_working_directory(repository_root: Path) -> None:
    current_directory = Path.cwd().resolve()
    _require(
        not _is_relative_to(current_directory, repository_root),
        "artifact verifier working directory is inside the source checkout; "
        "run it from a temporary directory outside the repository",
    )


def _join_resource_path(root: Any, relative_path: str) -> Any:
    parsed = PurePosixPath(relative_path)
    _require(
        not parsed.is_absolute() and ".." not in parsed.parts,
        f"packaged scientific resource path is unsafe: {relative_path!r}",
    )
    resource = root
    for part in parsed.parts:
        resource = resource.joinpath(part)
    return resource


def _required_manifest_text(value: object, *, field_name: str) -> str:
    _require(
        isinstance(value, str) and bool(value.strip()),
        f"packaged scientific resource manifest {field_name} must be non-empty text",
    )
    return value.strip()


def _mapping(value: object) -> dict[str, object]:
    _require(isinstance(value, Mapping), "expected a JSON object payload")
    return dict(value)


def _assert_close(
    actual: float,
    expected: float,
    *,
    label: str,
    abs_tol: float = 1.0e-9,
) -> None:
    if not math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=abs_tol):
        raise VerificationError(
            f"{label} mismatch: expected {expected:.12g}, got {actual:.12g}"
        )


def _assert_sequence_close(
    actual: Sequence[float],
    expected: Sequence[float],
    *,
    label: str,
    abs_tol: float = 1.0e-9,
) -> None:
    _require(
        len(actual) == len(expected),
        f"{label} length mismatch: expected {len(expected)}, got {len(actual)}",
    )
    for position, (actual_value, expected_value) in enumerate(
        zip(actual, expected, strict=True)
    ):
        _assert_close(
            float(actual_value),
            float(expected_value),
            label=f"{label}[{position}]",
            abs_tol=abs_tol,
        )


def _same_resolved_path(left: Path, right: Path) -> bool:
    return left.resolve() == right.resolve()


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _json_safe(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bool | int | str) or value is None:
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        try:
            return _json_safe(value.item())
        except Exception:
            return str(value)
    return str(value)


_CHECKS: tuple[tuple[str, CheckFunction], ...] = (
    ("artifact-manifest-binding", _check_artifact_manifest_binding),
    ("installed-package-identity", _check_installed_package_identity),
    ("packaged-scientific-resources", _check_packaged_scientific_resources),
    ("public-boundary-integrity", _check_public_boundary_integrity),
    (
        "corrected-construction-and-provenance-path",
        _check_corrected_construction_and_provenance_path,
    ),
    (
        "corrected-derived-and-ownership-path",
        _check_corrected_derived_and_ownership_path,
    ),
    ("corrected-differential-path", _check_corrected_differential_path),
)


if __name__ == "__main__":
    raise SystemExit(main())
