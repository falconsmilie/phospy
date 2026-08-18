"""Verify built PhosPy wheel and sdist by installing and executing them.

This script is intentionally independent of repository pytest fixtures. It
creates isolated virtual environments outside the checkout, installs exactly one
wheel and one sdist from ``dist/``, and runs an isolated Python probe from a
temporary working directory so source-tree imports cannot satisfy the check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from dataclasses import dataclass
from email.parser import Parser
from pathlib import Path
from typing import Any

SUPPORTED_PYTHON_VERSIONS = ("3.11", "3.12")
EXPECTED_REQUIRES_PYTHON = ">=3.11,<3.13"
EXPECTED_REQUIRES_PYTHON_SPECIFIERS = frozenset(EXPECTED_REQUIRES_PYTHON.split(","))
WHEEL_SUFFIX = ".whl"
SDIST_SUFFIX = ".tar.gz"

INSTALLED_PROBE_SOURCE = r"""
from __future__ import annotations

import hashlib
import importlib.resources as resources
import json
import pathlib
import sys
import tempfile

import pandas as pd

import phospy
import phospy.advanced as advanced
import phospy.api as api
import phospy.contracts.configs as contract_configs
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.advanced import DatasetMissingDataConfig, DatasetNormalisationConfig
from phospy.advanced import DatasetSiteMatrixConfig, KinaseReliabilityProfile
from phospy.advanced import PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
from phospy.advanced import KinaseScoringConfig
from phospy.advanced import ReferenceContextCompatibilityPolicy
from phospy.advanced import publish_dataset
from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy import KinaseWorkflow
from phospy.api import Contrast, DatasetBuildRequest, DatasetLocalisationConfig
from phospy.api import DatasetPreprocessingConfig, DifferentialAnalysisRequest
from phospy.api import ExperimentalDesign, KinaseWorkflowRequest, Organism
from phospy.api import ReferencePreset, SampleDesignRecord
from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot
from phospy.io.bundles.kinase import load_kinase_workflow_bundle
from phospy.io.bundles.kinase import save_kinase_workflow_bundle
from phospy.science.references.resources import load_bundled_kinase_substrate_map
from phospy.science.references.resources import load_bundled_motif_scores
from phospy.science.references.resources import load_bundled_motif_sizes
from phospy.science.references.resources import load_bundled_reference_manifest
from phospy.science.references.resources import load_bundled_site_sequences


REQUIRED_API_NAMES = (
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DatasetBuildRequest",
    "DatasetPreprocessingConfig",
    "DifferentialAnalysisRequest",
    "DifferentialAnalysisWorkflow",
    "ExperimentalDesign",
    "KinaseWorkflow",
    "KinaseWorkflowRequest",
    "Organism",
    "ReferencePreset",
    "SignalomeWorkflow",
)


def main() -> None:
    repo_root = pathlib.Path(sys.argv[1]).resolve()
    environment_root = pathlib.Path(sys.argv[2]).resolve()
    artifact_kind = sys.argv[3]

    package_file = pathlib.Path(phospy.__file__).resolve()
    _require_inside(package_file, environment_root, "phospy.__file__")
    _require_outside(package_file, repo_root, "phospy.__file__")

    public_surface_report = _verify_public_surface()
    resource_report = _verify_bundled_rat_reference_resources()
    scientific_report = _exercise_public_scientific_contracts()

    print(
        json.dumps(
            {
                "status": "ok",
                "artifact_kind": artifact_kind,
                "python": sys.version.split()[0],
                "executable": str(pathlib.Path(sys.executable).resolve()),
                "phospy_file": str(package_file),
                "public_surface_report": public_surface_report,
                "resource_report": resource_report,
                "scientific_report": scientific_report,
            },
            sort_keys=True,
        )
    )


def _verify_public_surface() -> dict[str, object]:
    for name in REQUIRED_API_NAMES:
        if name not in api.__all__:
            raise AssertionError(f"phospy.api.__all__ is missing {name!r}")
        getattr(api, name)
    for name in (
        "DatasetIntensityTransformConfig",
        "KinaseScoringConfig",
        "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION",
        "ReferenceContextCompatibilityPolicy",
        "publish_dataset",
    ):
        if name in api.__all__:
            raise AssertionError(f"advanced API name leaked into phospy.api: {name}")
        if name not in advanced.__all__:
            raise AssertionError(f"phospy.advanced.__all__ is missing {name!r}")
        getattr(advanced, name)
    if contract_configs.PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION != (
        "duplicate_correlation"
    ):
        raise AssertionError("contracts.configs duplicate-correlation policy changed")
    if PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION != "duplicate_correlation":
        raise AssertionError("advanced duplicate-correlation policy export changed")
    if "duplicate_correlation" not in contract_configs.SUPPORTED_PAIRED_DESIGN_POLICIES:
        raise AssertionError("supported paired-design policies omit duplicate_correlation")
    if phospy.AnalysisReadyDatasetBuilder is not api.AnalysisReadyDatasetBuilder:
        raise AssertionError("root AnalysisReadyDatasetBuilder is not api export")
    if phospy.DifferentialAnalysisWorkflow is not api.DifferentialAnalysisWorkflow:
        raise AssertionError("root DifferentialAnalysisWorkflow is not api export")
    if phospy.KinaseWorkflow is not api.KinaseWorkflow:
        raise AssertionError("root KinaseWorkflow is not api export")
    ticket_1_report = _verify_withdrawn_peptide_to_site_boundary()
    return {
        "required_api_names": list(REQUIRED_API_NAMES),
        "ticket_1_posthoc_peptide_to_site_boundary": ticket_1_report,
    }


def _verify_withdrawn_peptide_to_site_boundary() -> dict[str, object]:
    import phospy.science.differential as differential_public
    import phospy.science.differential.aggregation as aggregation_public
    from phospy.errors import PhosPyInputError
    from phospy.science.differential.aggregation.experimental import (
        PeptideToSiteAggregator,
    )

    withdrawn_public_exports = {
        "PeptideDifferentialEstimateTable",
        "PeptideToSiteAggregationConfig",
        "PeptideToSiteAggregationExecutor",
        "PeptideToSiteAggregationResult",
        "PeptideToSiteAggregator",
        "PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL",
        "PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT",
        "SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES",
    }
    for module_name, module in (
        ("phospy.science.differential", differential_public),
        ("phospy.science.differential.aggregation", aggregation_public),
    ):
        leaked = withdrawn_public_exports & set(module.__all__)
        if leaked:
            raise AssertionError(
                f"{module_name} still exports withdrawn post-hoc symbols: "
                f"{sorted(leaked)}"
            )

    if aggregation_public.PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS != (
        "unsupported_withdrawn_posthoc_estimate_combination_v1"
    ):
        raise AssertionError("peptide-to-site aggregation support status is not withdrawn")

    try:
        PeptideToSiteAggregator().run({"logFC": [-9.75]})
    except PhosPyInputError as exc:
        message = str(exc)
    else:
        raise AssertionError("withdrawn peptide-to-site compatibility shell executed")

    for required in (
        "withdrawn from public support",
        "coherent combined effect/inference",
        "executable peptide-to-site mapping semantics",
    ):
        if required not in message:
            raise AssertionError(
                "withdrawn peptide-to-site compatibility error is missing "
                f"{required!r}: {message}"
            )
    return {
        "status": "withdrawn_asserted",
        "support_status": aggregation_public.PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
        "checked_error_tokens": [
            "withdrawn from public support",
            "coherent combined effect/inference",
            "executable peptide-to-site mapping semantics",
        ],
    }


def _verify_bundled_rat_reference_resources() -> dict[str, object]:
    package_root = resources.files("phospy")
    bundle_root = (
        package_root.joinpath("data")
        .joinpath("reference_bundles")
        .joinpath("rat")
        .joinpath("l6_native")
    )
    manifest_resource = bundle_root.joinpath("manifest.json")
    if not manifest_resource.is_file():
        raise AssertionError("bundled rat l6_native manifest is missing")
    manifest_payload = json.loads(manifest_resource.read_text(encoding="utf-8"))
    if manifest_payload.get("reference_id") != "l6_native":
        raise AssertionError("bundled rat manifest reference_id changed")

    manifest = load_bundled_reference_manifest(Organism.RAT)
    if manifest.reference_id != "l6_native":
        raise AssertionError("typed bundled rat manifest reference_id changed")
    if manifest.reference_version != manifest_payload.get("reference_version"):
        raise AssertionError("typed manifest version does not match resource JSON")

    verified_files: dict[str, str] = {}
    for file_manifest in manifest_payload["files"]:
        relative_path = file_manifest["relative_path"]
        declared_hash = file_manifest["sha256"]
        resource = bundle_root.joinpath(relative_path)
        if not resource.is_file():
            raise AssertionError(f"manifest-declared resource is missing: {relative_path}")
        actual_hash = hashlib.sha256(resource.read_bytes()).hexdigest()
        if actual_hash != declared_hash:
            raise AssertionError(
                "manifest-declared resource hash mismatch: "
                f"{relative_path}; expected={declared_hash}; actual={actual_hash}"
            )
        verified_files[str(relative_path)] = actual_hash

    substrate_map = load_bundled_kinase_substrate_map(Organism.RAT)
    site_sequences = load_bundled_site_sequences(Organism.RAT)
    motif_scores = load_bundled_motif_scores(Organism.RAT)
    motif_sizes = load_bundled_motif_sizes(Organism.RAT)
    if substrate_map.empty or site_sequences.empty:
        raise AssertionError("bundled rat reference core tables must be non-empty")
    if motif_scores is None or motif_scores.empty:
        raise AssertionError("bundled rat motif_scores must be present and non-empty")
    if motif_sizes is None or motif_sizes.empty:
        raise AssertionError("bundled rat motif_sizes must be present and non-empty")

    return {
        "reference_id": manifest.reference_id,
        "reference_version": manifest.reference_version,
        "verified_resource_count": len(verified_files),
        "verified_resource_sha256": verified_files,
        "substrate_rows": int(substrate_map.shape[0]),
        "site_sequence_rows": int(site_sequences.shape[0]),
        "motif_score_shape": [int(motif_scores.shape[0]), int(motif_scores.shape[1])],
        "motif_size_count": int(motif_sizes.shape[0]),
    }


def _exercise_public_scientific_contracts() -> dict[str, object]:
    imputed_dataset = _build_dataset(include_missing=True)
    if int(imputed_dataset.phospho.isna().sum().sum()) != 0:
        raise AssertionError("row-median preprocessing left missing values")
    if imputed_dataset.processing_state.missing_data.complete_matrix is not True:
        raise AssertionError("missing-data processing state is not complete")
    if imputed_dataset.processing_state.missing_data.imputed is not True:
        raise AssertionError("missing-data processing state did not record imputation")

    dataset = _build_dataset(include_missing=False)
    if dataset.phospho.index.name != "site_key":
        raise AssertionError("dataset builder did not establish site_key identity")
    if dataset.processing_state.missing_data.imputed is not False:
        raise AssertionError("non-missing dataset incorrectly records imputation")

    differential_result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_design(),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
        )
    )
    differential_table = differential_result.table_for("B_vs_A")
    required_columns = {
        "site_key",
        "display_id",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
    }
    if not required_columns.issubset(differential_table.columns):
        raise AssertionError("differential result table is missing public columns")
    adjusted = differential_table["adj.P.Val"].dropna()
    if not adjusted.between(0.0, 1.0, inclusive="both").all():
        raise AssertionError("differential adjusted p-values are outside [0, 1]")
    duplicate_report = _exercise_duplicate_correlation_differential_contract(dataset)

    kinase_request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            reliability_profile=KinaseReliabilityProfile.CUSTOM,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        activity_config=None,
        site_sequence_conflict_policy="prefer_reference",
    )
    kinase_result = KinaseWorkflow().run(kinase_request)
    if kinase_result.scoring_result.profile_scores.empty:
        raise AssertionError("kinase profile scores are empty")
    if kinase_result.prediction_result.pred_mat.empty:
        raise AssertionError("kinase prediction matrix is empty")
    if kinase_result.references.manifest is None:
        raise AssertionError("kinase bundled reference manifest was not retained")
    if kinase_result.references.manifest.reference_id != "l6_native":
        raise AssertionError("kinase workflow did not use the rat l6_native manifest")
    bundle_report = _exercise_result_bundle_round_trip(
        kinase_request=kinase_request,
        kinase_result=kinase_result,
    )
    publisher_report = _exercise_advanced_table_publisher(dataset)

    return {
        "imputed_dataset_shape": [
            int(imputed_dataset.phospho.shape[0]),
            int(imputed_dataset.phospho.shape[1]),
        ],
        "dataset_shape": [int(dataset.phospho.shape[0]), int(dataset.phospho.shape[1])],
        "differential_rows": int(differential_table.shape[0]),
        "duplicate_correlation": duplicate_report,
        "kinase_profile_shape": [
            int(kinase_result.scoring_result.profile_scores.shape[0]),
            int(kinase_result.scoring_result.profile_scores.shape[1]),
        ],
        "kinase_prediction_shape": [
            int(kinase_result.prediction_result.pred_mat.shape[0]),
            int(kinase_result.prediction_result.pred_mat.shape[1]),
        ],
        "bundle_round_trip": bundle_report,
        "advanced_table_publisher": publisher_report,
    }


def _exercise_duplicate_correlation_differential_contract(dataset) -> dict[str, object]:
    duplicate_result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_paired_design(),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
            config=advanced.DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
            ),
        )
    )
    if type(duplicate_result).__name__ != "DifferentialAnalysisResult":
        raise AssertionError("duplicate_correlation returned unexpected result type")
    duplicate_table = duplicate_result.table_for("B_vs_A")
    if {"logFC", "t", "P.Value", "adj.P.Val"}.difference(duplicate_table.columns):
        raise AssertionError("duplicate_correlation result table lost public columns")
    if duplicate_result.policy_provenance is None:
        raise AssertionError("duplicate_correlation result lacks policy provenance")
    duplicate = duplicate_result.policy_provenance.duplicate_correlation
    if duplicate is None:
        raise AssertionError("duplicate_correlation provenance was not attached")
    if duplicate.normalised_paired_design_policy != "duplicate_correlation":
        raise AssertionError("duplicate_correlation provenance recorded wrong policy")
    if duplicate.block_treatment != "consensus_correlation":
        raise AssertionError("duplicate_correlation block treatment changed")
    if duplicate.covariance_structure != "compound_symmetry":
        raise AssertionError("duplicate_correlation covariance structure changed")
    if duplicate.gls_fit_status != "fit":
        raise AssertionError("duplicate_correlation GLS did not report fit status")
    if duplicate.consensus.consensus_correlation is None:
        raise AssertionError("duplicate_correlation consensus was not recorded")
    if duplicate_result.policy_provenance.design.block_column_names:
        raise AssertionError("duplicate_correlation design contains block columns")
    if not any(
        caveat.code == "differential_duplicate_correlation_consensus"
        for caveat in duplicate_result.caveats
    ):
        raise AssertionError("duplicate_correlation consensus caveat is missing")
    return {
        "result_type": type(duplicate_result).__name__,
        "rows": int(duplicate_table.shape[0]),
        "policy": duplicate.normalised_paired_design_policy,
        "block_treatment": duplicate.block_treatment,
        "covariance_structure": duplicate.covariance_structure,
        "sample_count": int(duplicate.sample_count),
        "block_count": int(duplicate.block_count),
        "consensus_correlation": float(duplicate.consensus.consensus_correlation),
    }


def _exercise_advanced_table_publisher(dataset) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="phospy-installed-publisher-") as tmp_dir:
        output_root = pathlib.Path(tmp_dir) / "published_tables"
        written = publish_dataset(dataset, output_root, output_format="csv")
        expected_manifest = output_root / "dataset" / "manifest.json"
        if written["dataset.manifest"] != expected_manifest:
            raise AssertionError("publisher returned an unexpected manifest path")
        if not expected_manifest.is_file():
            raise AssertionError("publisher did not write the dataset manifest")
        stale_path = output_root / "dataset" / "stale.txt"
        stale_path.write_text("stale", encoding="utf-8")
        try:
            publish_dataset(dataset, output_root, output_format="csv")
        except api.PhosPyInputError as exc:
            if str(output_root) not in str(exc):
                raise AssertionError(
                    "publisher existing-root error did not identify destination"
                ) from exc
        else:
            raise AssertionError("publisher replaced existing output without overwrite=True")

        overwritten = publish_dataset(
            dataset,
            output_root,
            output_format="csv",
            overwrite=True,
        )
        if stale_path.exists():
            raise AssertionError("publisher overwrite left a stale file behind")
        if list(pathlib.Path(tmp_dir).glob(".published_tables.tmp-*")):
            raise AssertionError("publisher left a staging directory after success")
        if list(pathlib.Path(tmp_dir).glob(".published_tables.previous-*")):
            raise AssertionError("publisher left a backup directory after success")
        if any(
            ".tmp-" in part or ".previous-" in part
            for path in overwritten.values()
            for part in path.parts
        ):
            raise AssertionError("publisher returned staging or backup paths")
        manifest_payload = json.loads(expected_manifest.read_text(encoding="utf-8"))
        return {
            "manifest_output_format": manifest_payload["output_format"],
            "written_key_count": len(overwritten),
        }


def _exercise_result_bundle_round_trip(
    *,
    kinase_request: KinaseWorkflowRequest,
    kinase_result,
) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="phospy-installed-result-bundle-") as tmp_dir:
        bundle_root = pathlib.Path(tmp_dir) / "kinase_bundle"
        config_snapshot = KinaseWorkflowConfigSnapshot.from_request(kinase_request)
        written = save_kinase_workflow_bundle(
            kinase_result,
            bundle_root,
            config_snapshot=config_snapshot,
            output_format="csv",
        )
        if written["manifest"] != bundle_root / "manifest.json":
            raise AssertionError("bundle writer returned an unexpected manifest path")
        loaded = load_kinase_workflow_bundle(bundle_root)
        if loaded.config_snapshot != config_snapshot:
            raise AssertionError("installed artifact bundle config snapshot changed")
        if loaded.result.provenance != kinase_result.provenance:
            raise AssertionError("installed artifact bundle provenance changed")
        if loaded.result.prediction_result.pred_mat.shape != kinase_result.prediction_result.pred_mat.shape:
            raise AssertionError("installed artifact bundle prediction shape changed")
        manifest_payload = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
        phospho_entry = manifest_payload["dataset"]["tables"]["phospho"]
        if not isinstance(phospho_entry, dict):
            raise AssertionError("installed artifact bundle table entry is not content-addressed")
        for required_key in ("path", "sha256", "byte_size", "logical_type", "shape"):
            if required_key not in phospho_entry:
                raise AssertionError(f"bundle table entry is missing {required_key!r}")
        return {
            "manifest_version": int(loaded.manifest_version),
            "phospho_entry_logical_type": phospho_entry["logical_type"],
            "phospho_entry_shape": phospho_entry["shape"],
        }


def _build_dataset(*, include_missing: bool):
    index = ["TSC2;S939;", "GSK3B;S9;", "MAPK14;Y182;", "GSK3A;S21;"]
    phospho = pd.DataFrame(
        {
            "A_1": [10.0, 9.0, 7.5, 8.0],
            "A_2": [10.4, 9.1, 7.7, 8.2],
            "B_1": [12.0, 9.3, 7.4, 9.1],
            "B_2": [12.2, 9.4, 7.6, 9.0],
        },
        index=index,
    )
    if include_missing:
        phospho.iloc[1, 1] = None
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3B", "MAPK14", "GSK3A"],
            "site": ["S939", "S9", "Y182", "S21"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "_______MSGRPRTTSFAESCKPVQQPSAFG",
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            ],
            "display_id": index,
            "organism": ["rat"] * 4,
            "protein_namespace": ["protein_id"] * 4,
            "protein_identifier": ["TSC2", "GSK3B", "MAPK14", "GSK3A"],
            "localisation_confidence": [0.95] * 4,
            "protein_id": ["TSC2", "GSK3B", "MAPK14", "GSK3A"],
        },
        index=index,
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                localisation=DatasetLocalisationConfig(
                    mode="require_threshold",
                    confidence_column="localisation_confidence",
                    min_confidence=0.75,
                ),
                intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
                normalisation=DatasetNormalisationConfig(policy="median_center"),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                    input_scale="linear",
                ),
                site_matrix=DatasetSiteMatrixConfig(policy="as_input"),
            ),
        )
    )


def _design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_2",
            ),
        )
    )


def _paired_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_2",
                block_id="pair_2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_2",
                block_id="pair_2",
            ),
        )
    )


def _require_inside(path: pathlib.Path, root: pathlib.Path, label: str) -> None:
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise AssertionError(f"{label} is outside expected root: {path} not in {root}") from exc


def _require_outside(path: pathlib.Path, root: pathlib.Path, label: str) -> None:
    try:
        path.relative_to(root)
    except ValueError:
        return
    raise AssertionError(f"{label} unexpectedly resolves inside checkout: {path}")


main()
"""


class InstalledDistributionVerificationError(RuntimeError):
    """Built distribution installation or installed-package execution failed."""


@dataclass(frozen=True)
class DistributionArtifacts:
    wheel: Path
    sdist: Path


@dataclass(frozen=True)
class InstalledDistributionReport:
    artifact_kind: str
    artifact_path: Path
    artifact_sha256: str
    environment_root: Path
    run_directory: Path
    phospy_file: Path
    python_version: str
    requires_python: str
    resource_count: int
    ticket_1_boundary_status: str
    constraint_path: Path | None
    constraint_sha256: str | None


def find_distribution_artifacts(dist_dir: str | Path) -> DistributionArtifacts:
    """Return exactly one wheel and one sdist from a distribution directory."""

    root = Path(dist_dir)
    wheels = sorted(path for path in root.glob(f"*{WHEEL_SUFFIX}") if path.is_file())
    sdists = sorted(path for path in root.glob(f"*{SDIST_SUFFIX}") if path.is_file())
    if len(wheels) != 1:
        raise InstalledDistributionVerificationError(
            f"expected exactly one wheel in {root}, found {len(wheels)}"
        )
    if len(sdists) != 1:
        raise InstalledDistributionVerificationError(
            f"expected exactly one sdist in {root}, found {len(sdists)}"
        )
    return DistributionArtifacts(wheel=wheels[0].resolve(), sdist=sdists[0].resolve())


def verify_installed_distributions(
    *,
    dist_dir: str | Path,
    repo_root: str | Path,
    python_executable: str | Path = sys.executable,
    constraint: str | Path | None = None,
    use_system_site_packages: bool = False,
    install_dependencies: bool = True,
    build_isolation: bool = True,
    install_packaging_tools: bool = True,
    ignore_requires_python: bool = False,
) -> tuple[InstalledDistributionReport, ...]:
    """Install and execute the built wheel and sdist in isolated environments."""

    artifacts = find_distribution_artifacts(dist_dir)
    resolved_repo_root = Path(repo_root).resolve()
    resolved_python = str(python_executable)
    resolved_constraint = None if constraint is None else Path(constraint).resolve()

    reports: list[InstalledDistributionReport] = []
    with tempfile.TemporaryDirectory(prefix="phospy-installed-dist-") as tmp_dir:
        work_root = Path(tmp_dir).resolve()
        _require_path_outside(
            work_root, resolved_repo_root, label="verification tempdir"
        )
        for artifact_kind, artifact_path in (
            ("wheel", artifacts.wheel),
            ("sdist", artifacts.sdist),
        ):
            reports.append(
                _verify_one_artifact(
                    artifact_kind=artifact_kind,
                    artifact_path=artifact_path,
                    work_root=work_root,
                    repo_root=resolved_repo_root,
                    python_executable=resolved_python,
                    constraint=resolved_constraint,
                    use_system_site_packages=use_system_site_packages,
                    install_dependencies=install_dependencies,
                    build_isolation=build_isolation,
                    install_packaging_tools=install_packaging_tools,
                    ignore_requires_python=ignore_requires_python,
                )
            )
    return tuple(reports)


def _verify_one_artifact(
    *,
    artifact_kind: str,
    artifact_path: Path,
    work_root: Path,
    repo_root: Path,
    python_executable: str,
    constraint: Path | None,
    use_system_site_packages: bool = False,
    install_dependencies: bool = True,
    build_isolation: bool = True,
    install_packaging_tools: bool = True,
    ignore_requires_python: bool = False,
) -> InstalledDistributionReport:
    requires_python = _validate_artifact_requires_python(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
    )
    environment_root = work_root / f"{artifact_kind}-venv"
    run_directory = work_root / f"{artifact_kind}-run"
    run_directory.mkdir()

    venv_command = [python_executable, "-m", "venv"]
    if use_system_site_packages:
        venv_command.append("--system-site-packages")
    venv_command.append(str(environment_root))
    _run(
        venv_command,
        cwd=work_root,
        repo_root=repo_root,
        context=f"create {artifact_kind} verification environment",
    )
    environment_python = _environment_python(environment_root)
    if not environment_python.is_file():
        raise InstalledDistributionVerificationError(
            f"virtual environment Python is missing: {environment_python}"
        )

    if install_packaging_tools:
        _run_pip(
            environment_python,
            "install packaging tools",
            repo_root=repo_root,
            cwd=work_root,
            constraint=constraint,
            arguments=["install", "--upgrade", "pip", "setuptools", "wheel"],
        )
    _run_pip(
        environment_python,
        f"install {artifact_kind}",
        repo_root=repo_root,
        cwd=work_root,
        constraint=constraint,
        arguments=[
            "install",
            *(("--no-deps",) if not install_dependencies else ()),
            *(("--no-build-isolation",) if not build_isolation else ()),
            *(("--ignore-requires-python",) if ignore_requires_python else ()),
            str(artifact_path),
        ],
    )
    _run_pip(
        environment_python,
        f"check {artifact_kind} environment",
        repo_root=repo_root,
        cwd=work_root,
        constraint=None,
        arguments=["check"],
    )
    result = _run(
        [
            str(environment_python),
            "-I",
            "-c",
            INSTALLED_PROBE_SOURCE,
            str(repo_root),
            str(environment_root.resolve()),
            artifact_kind,
        ],
        cwd=run_directory,
        repo_root=repo_root,
        context=f"execute installed {artifact_kind} probe",
    )
    payload = _probe_payload(result.stdout, artifact_kind=artifact_kind)
    phospy_file = Path(_required_string(payload, "phospy_file")).resolve()
    _require_path_inside(
        phospy_file, environment_root.resolve(), label="phospy.__file__"
    )
    _require_path_outside(phospy_file, repo_root, label="phospy.__file__")
    resource_report = payload.get("resource_report")
    if not isinstance(resource_report, dict):
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe did not report resource results"
        )
    resource_count = resource_report.get("verified_resource_count")
    if not isinstance(resource_count, int) or resource_count < 1:
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe verified no resources"
        )
    ticket_1_boundary_status = _ticket_1_boundary_status(
        payload,
        artifact_kind=artifact_kind,
    )
    return InstalledDistributionReport(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
        artifact_sha256=_sha256_file(artifact_path),
        environment_root=environment_root.resolve(),
        run_directory=run_directory.resolve(),
        phospy_file=phospy_file,
        python_version=_required_string(payload, "python"),
        requires_python=requires_python,
        resource_count=int(resource_count),
        ticket_1_boundary_status=ticket_1_boundary_status,
        constraint_path=None if constraint is None else constraint.resolve(),
        constraint_sha256=None if constraint is None else _sha256_file(constraint),
    )


def _validate_artifact_requires_python(
    *,
    artifact_kind: str,
    artifact_path: Path,
) -> str:
    requires_python = _artifact_requires_python(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
    )
    if (
        _requires_python_specifiers(requires_python)
        != EXPECTED_REQUIRES_PYTHON_SPECIFIERS
    ):
        raise InstalledDistributionVerificationError(
            f"{artifact_kind} metadata Requires-Python mismatch for {artifact_path}: "
            f"expected {EXPECTED_REQUIRES_PYTHON!r}, found {requires_python!r}"
        )
    return requires_python


def _artifact_requires_python(*, artifact_kind: str, artifact_path: Path) -> str:
    if artifact_kind == "wheel":
        metadata_text = _wheel_metadata_text(artifact_path)
    elif artifact_kind == "sdist":
        metadata_text = _sdist_metadata_text(artifact_path)
    else:
        raise InstalledDistributionVerificationError(
            f"unsupported artifact kind for metadata verification: {artifact_kind!r}"
        )
    value = Parser().parsestr(metadata_text).get("Requires-Python")
    if value is None or not value.strip():
        raise InstalledDistributionVerificationError(
            f"{artifact_kind} metadata is missing Requires-Python: {artifact_path}"
        )
    return value.strip()


def _requires_python_specifiers(value: str) -> frozenset[str]:
    return frozenset(part.strip() for part in value.split(",") if part.strip())


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _wheel_metadata_text(artifact_path: Path) -> str:
    try:
        with zipfile.ZipFile(artifact_path) as wheel:
            metadata_names = sorted(
                name
                for name in wheel.namelist()
                if name.endswith(".dist-info/METADATA")
            )
            if len(metadata_names) != 1:
                raise InstalledDistributionVerificationError(
                    f"expected exactly one wheel METADATA file in {artifact_path}, "
                    f"found {len(metadata_names)}"
                )
            return wheel.read(metadata_names[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise InstalledDistributionVerificationError(
            f"wheel is not a readable zip archive: {artifact_path}"
        ) from exc
    except UnicodeDecodeError as exc:
        raise InstalledDistributionVerificationError(
            f"wheel METADATA is not UTF-8 decodable: {artifact_path}"
        ) from exc


def _sdist_metadata_text(artifact_path: Path) -> str:
    try:
        with tarfile.open(artifact_path, "r:gz") as sdist:
            metadata_members = sorted(
                (
                    member
                    for member in sdist.getmembers()
                    if member.isfile() and member.name.endswith("/PKG-INFO")
                ),
                key=lambda member: member.name,
            )
            top_level_members = [
                member for member in metadata_members if member.name.count("/") == 1
            ]
            candidates = top_level_members or metadata_members
            if len(candidates) != 1:
                raise InstalledDistributionVerificationError(
                    f"expected exactly one sdist PKG-INFO file in {artifact_path}, "
                    f"found {len(candidates)}"
                )
            extracted = sdist.extractfile(candidates[0])
            if extracted is None:
                raise InstalledDistributionVerificationError(
                    f"sdist PKG-INFO could not be read: {candidates[0].name}"
                )
            return extracted.read().decode("utf-8")
    except tarfile.TarError as exc:
        raise InstalledDistributionVerificationError(
            f"sdist is not a readable gzipped tar archive: {artifact_path}"
        ) from exc
    except UnicodeDecodeError as exc:
        raise InstalledDistributionVerificationError(
            f"sdist PKG-INFO is not UTF-8 decodable: {artifact_path}"
        ) from exc


def _ticket_1_boundary_status(
    payload: dict[str, Any],
    *,
    artifact_kind: str,
) -> str:
    public_surface_report = payload.get("public_surface_report")
    if not isinstance(public_surface_report, dict):
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe did not report public-surface results"
        )
    ticket_1_report = public_surface_report.get(
        "ticket_1_posthoc_peptide_to_site_boundary"
    )
    if not isinstance(ticket_1_report, dict):
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe did not report Ticket 1 boundary results"
        )
    status = ticket_1_report.get("status")
    if status != "withdrawn_asserted":
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe reported unsupported Ticket 1 "
            f"boundary status: {ticket_1_report!r}"
        )
    return str(status)


def _run_pip(
    environment_python: Path,
    context: str,
    *,
    repo_root: Path,
    cwd: Path,
    constraint: Path | None,
    arguments: list[str],
) -> subprocess.CompletedProcess[str]:
    command = [
        str(environment_python),
        "-m",
        "pip",
        "--disable-pip-version-check",
        arguments[0],
    ]
    environment_overrides = None
    if constraint is not None and arguments[0] == "install":
        command.extend(["-c", str(constraint)])
        environment_overrides = {"PIP_BUILD_CONSTRAINT": str(constraint)}
    command.extend(arguments[1:])
    return _run(
        command,
        cwd=cwd,
        repo_root=repo_root,
        context=context,
        environment_overrides=environment_overrides,
    )


def _run(
    command: list[str],
    *,
    cwd: Path,
    repo_root: Path,
    context: str,
    environment_overrides: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    _require_path_outside(cwd.resolve(), repo_root, label=f"{context} cwd")
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    env.pop("PIP_CONSTRAINT", None)
    env.pop("PIP_BUILD_CONSTRAINT", None)
    if environment_overrides is not None:
        env.update(environment_overrides)
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise InstalledDistributionVerificationError(
            f"{context} failed with exit code {result.returncode}\n"
            f"command: {_format_command(command)}\n\n"
            f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
        )
    return result


def _probe_payload(stdout: str, *, artifact_kind: str) -> dict[str, Any]:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines:
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe produced no JSON output"
        )
    try:
        payload = json.loads(lines[-1])
    except json.JSONDecodeError as exc:
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe produced invalid JSON: {lines[-1]!r}"
        ) from exc
    if not isinstance(payload, dict):
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe JSON is not an object"
        )
    if payload.get("status") != "ok":
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} probe did not report ok status: {payload!r}"
        )
    if payload.get("artifact_kind") != artifact_kind:
        raise InstalledDistributionVerificationError(
            f"installed probe reported wrong artifact kind: {payload!r}"
        )
    return payload


def _required_string(payload: dict[str, Any], key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise InstalledDistributionVerificationError(
            f"installed probe payload field {key!r} must be a non-empty string"
        )
    return value


def _environment_python(environment_root: Path) -> Path:
    if os.name == "nt":
        return environment_root / "Scripts" / "python.exe"
    return environment_root / "bin" / "python"


def _require_path_inside(path: Path, root: Path, *, label: str) -> None:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise InstalledDistributionVerificationError(
            f"{label} is outside expected root: {path.resolve()} not in {root.resolve()}"
        ) from exc


def _require_path_outside(path: Path, root: Path, *, label: str) -> None:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return
    raise InstalledDistributionVerificationError(
        f"{label} unexpectedly resolves inside checkout: {path.resolve()}"
    )


def _format_command(command: list[str]) -> str:
    if os.name == "nt":
        return subprocess.list2cmdline(command)
    return shlex.join(command)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Install exactly one built PhosPy wheel and one built PhosPy sdist "
            "outside the checkout, then run isolated installed-package checks."
        )
    )
    parser.add_argument(
        "--dist-dir",
        default="dist",
        help="Directory containing exactly one .whl and one .tar.gz artifact.",
    )
    parser.add_argument(
        "--repo-root",
        default=Path(__file__).resolve().parents[1],
        type=Path,
        help="Repository root to exclude from installed import origins.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used to create verification environments.",
    )
    parser.add_argument(
        "--constraint",
        type=Path,
        default=None,
        help="Optional pip constraint file used for artifact installations.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        reports = verify_installed_distributions(
            dist_dir=args.dist_dir,
            repo_root=args.repo_root,
            python_executable=args.python,
            constraint=args.constraint,
        )
    except InstalledDistributionVerificationError as exc:
        print(f"installed distribution verification failed: {exc}", file=sys.stderr)
        return 1

    print(
        json.dumps(
            {
                "status": "ok",
                "supported_python_versions": SUPPORTED_PYTHON_VERSIONS,
                "dependency_constraints": (
                    None
                    if args.constraint is None
                    else {
                        "path": str(Path(args.constraint).resolve()),
                        "sha256": _sha256_file(Path(args.constraint).resolve()),
                    }
                ),
                "verified": [
                    {
                        "artifact_kind": report.artifact_kind,
                        "artifact_path": str(report.artifact_path),
                        "artifact_sha256": report.artifact_sha256,
                        "python": report.python_version,
                        "requires_python": report.requires_python,
                        "phospy_file": str(report.phospy_file),
                        "resource_count": report.resource_count,
                        "ticket_1_boundary_status": report.ticket_1_boundary_status,
                        "constraint_path": (
                            None
                            if report.constraint_path is None
                            else str(report.constraint_path)
                        ),
                        "constraint_sha256": report.constraint_sha256,
                    }
                    for report in reports
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
