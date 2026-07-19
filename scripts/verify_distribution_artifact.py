from __future__ import annotations

import argparse
import copy
import hashlib
import importlib
import json
import math
import platform
import sys
import time
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
    manifest_payload = json.loads(manifest_resource.read_text(encoding="utf-8"))
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
        **digest_report,
    }


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
    return {
        "trusted_assertion_fingerprint": (
            dataset.trusted_construction_assertions.assertion_fingerprint
        ),
        "provenance_hash": before_hash,
        "round_trip_verified": True,
        "organism_reference_contradiction_rejected": True,
    }


def _check_corrected_derived_and_ownership_path(
    context: VerificationContext,
) -> Mapping[str, object]:
    from phospy.api import Organism
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
    from phospy.provenance import fingerprint_optional_table

    fingerprints = [
        fingerprint_optional_table(phospho, name="dataset.phospho"),
        fingerprint_optional_table(site_metadata, name="dataset.site_metadata"),
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
