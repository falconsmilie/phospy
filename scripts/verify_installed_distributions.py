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
from dataclasses import dataclass, field
from email.parser import Parser
from pathlib import Path
from typing import Any

SUPPORTED_PYTHON_VERSIONS = ("3.11", "3.12")
EXPECTED_REQUIRES_PYTHON = ">=3.11,<3.13"
EXPECTED_REQUIRES_PYTHON_SPECIFIERS = frozenset(EXPECTED_REQUIRES_PYTHON.split(","))
WHEEL_SUFFIX = ".whl"
SDIST_SUFFIX = ".tar.gz"

RETIRED_CAPABILITY_ROUTES = (
    "phospy.contracts.kinase_reference_projection",
    "phospy.science.prediction.motif_scoring",
    "phospy.science.prediction.motif_scoring.frequency_matrices",
    "phospy.science.prediction.motif_scoring.kinase_library",
    "phospy.science.prediction.motif_scoring.library_validation",
    "phospy.science.prediction.motif_scoring.models",
    "phospy.science.prediction.motif_scoring.scaling",
    "phospy.science.prediction.motif_scoring.scoring",
    "phospy.science.prediction.scoring",
    "phospy.science.references.kinase_library",
    "phospy.science.references.kinase_library_loading",
    "phospy.science.references.kinase_library_models",
    "phospy.science.references.validation.kinase_library",
    "phospy.science.tables.kinase",
    "phospy.validation.references.kinase_library",
    "phospy.workflows.kinase.contributions",
    "phospy.workflows.kinase.kinase_library_scoring",
    "phospy.workflows.kinase.reference_projection",
    "phospy.workflows.kinase.scoring_mode_contracts",
    "phospy.workflows.kinase.scoring_runner",
    "phospy.workflows.kinase.scoring_transforms",
    "phospy.workflows.kinase",
    "phospy.contracts.results.kinase",
    "phospy.contracts.results.kinase_boundary",
)

RETIRED_HORIZONTAL_ROUTES = (
    "phospy.contracts",
    "phospy.contracts.configs",
    "phospy.contracts.requests",
    "phospy.contracts.results",
    "phospy.contracts.results.base",
    "phospy.science",
    "phospy.science.configs",
    "phospy.science.datasets",
    "phospy.science.references",
    "phospy.science.references.validation",
    "phospy.science.tables",
    "phospy.workflows",
    "phospy.validation.configs",
    "phospy.validation.configs.preprocessing",
    "phospy.validation.datasets",
    "phospy.validation.datasets.analysis_ready",
    "phospy.validation.datasets.fragpipe",
    "phospy.validation.datasets.group_coverage_filter",
    "phospy.validation.datasets.importers",
    "phospy.validation.datasets.maxquant",
    "phospy.validation.references",
    "phospy.validation.workflows",
)

RETIRED_ROUTE_IMPORT_PROBE_SOURCE = r"""
from __future__ import annotations

import importlib
import json
import sys

route = sys.argv[1]
try:
    importlib.import_module(route)
except ModuleNotFoundError as error:
    missing_name = error.name
    requested_route_or_absent_ancestor = isinstance(missing_name, str) and (
        missing_name == route or route.startswith(missing_name + ".")
    )
    print(
        json.dumps(
            {
                "route": route,
                "status": (
                    "absent"
                    if requested_route_or_absent_ancestor
                    else "unrelated_missing_dependency"
                ),
                "exception_type": "ModuleNotFoundError",
                "module_not_found_error_name": missing_name,
                "requested_route_or_absent_ancestor": (
                    requested_route_or_absent_ancestor
                ),
            },
            sort_keys=True,
        )
    )
    raise SystemExit(0 if requested_route_or_absent_ancestor else 2)
else:
    print(json.dumps({"route": route, "status": "present"}, sort_keys=True))
    raise SystemExit(1)
"""

INSTALLED_PROBE_SOURCE = (
    r"""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata as metadata
import importlib.util
import importlib.resources as resources
import inspect
import json
import pathlib
import platform
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import fields
from enum import Enum
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

import phospy
import phospy.advanced as advanced
import phospy.api as api
import phospy._api_inventory as api_inventory
import phospy.evidence.config as evidence_config
import phospy.kinase.activity.config as kinase_activity_config
import phospy.kinase.prediction.config as kinase_prediction_config
import phospy.kinase.scoring.config as kinase_scoring_config
import phospy.preprocessing.batch_correction.contracts.configuration as batch_config
import phospy.preprocessing.batch_correction.contracts.missingness as missingness_config
import phospy.preprocessing.contracts.configuration as preprocessing_config
import phospy.references.config as references_config
import phospy.signalome.config as signalome_config
import phospy.kinase.policies as kinase_policies
import phospy.kinase.caveats as kinase_caveats
import phospy.differential.policies as differential_policies
from phospy.advanced import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1,
)
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.advanced import DatasetMissingDataConfig, DatasetNormalisationConfig
from phospy.advanced import DatasetProteinAwarePreparationConfig
from phospy.advanced import DatasetSiteMatrixConfig, KinaseReliabilityProfile
from phospy.advanced import DifferentialProteinAwareModelConfig
from phospy.advanced import PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
from phospy.advanced import KinaseScoringConfig
from phospy.advanced import LocalisationRequirement
from phospy.advanced import ReferenceContextCompatibilityPolicy
from phospy.advanced import SignalomeClusteringConfig, SignalomeConfig
from phospy.advanced import SignalomeOutputConfig, SignalomeScientificConfig
from phospy.advanced import SignalomeValidationConfig
from phospy.advanced import (
    ControlSiteSet,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    ObservationMask,
    OriginallyMissingCellTracking,
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryResult,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
    SpsRuvBatchCorrectionConfig,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.advanced import publish_dataset, publish_kinase_workflow
from phospy.advanced import publish_signalome_workflow
from phospy import AnalysisReadyDatasetBuilder
from phospy import KinaseWorkflow
from phospy import SignalomeWorkflow
from phospy.api import Contrast, DatasetBuildRequest, DatasetLocalisationConfig
from phospy.api import DatasetPreprocessingConfig, DifferentialAnalysisRequest
from phospy.api import DifferentialAnalysisWorkflow
from phospy.api import ExperimentalDesign, KinaseWorkflowRequest, Organism
from phospy.api import ReferencePreset, SampleDesignRecord
from phospy.kinase.result_boundary import (
    KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES,
    KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
    KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
    KINASE_SCORING_LIMITATION_CAVEAT_CODE,
    SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    DownstreamScoreSelectionPolicy,
    DownstreamScoreSource,
    KinaseScoreMatrix,
)
from phospy.kinase.prediction.config import (
    KinaseAdaptivePolicy as OwnerKinaseAdaptivePolicy,
)
from phospy.kinase.prediction.config import (
    KinasePredictionConfig as OwnerKinasePredictionConfig,
)
from phospy.kinase.prediction.config import (
    KinasePredictionMode as OwnerKinasePredictionMode,
)
from phospy.kinase.prediction.results import (
    KinasePredictionResult as OwnerKinasePredictionResult,
)
from phospy.kinase.prediction.tables import KinasePredictionMatrix
from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot
from phospy.io.importers.contracts import (
    PhosphositeImporter,
    PhosphositeImportRequest,
)
from phospy.signalome.requests import SignalomeWorkflowRequest
from phospy.io.bundles.kinase import load_kinase_workflow_bundle
from phospy.io.bundles.kinase import save_kinase_workflow_bundle
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot
from phospy.io.bundles.signalome import load_signalome_workflow_bundle
from phospy.io.bundles.signalome import save_signalome_workflow_bundle
from phospy.references.resources import load_bundled_kinase_substrate_map
from phospy.references.resources import load_bundled_motif_scores
from phospy.references.resources import load_bundled_motif_sizes
from phospy.references.resources import load_bundled_reference_manifest
from phospy.references.resources import load_bundled_site_sequences
from phospy.references.builder import ReferenceBundleBuilder as OwnerReferenceBundleBuilder
from phospy.references.provenance.serialization import (
    reference_from_payload,
    reference_to_payload,
)
from phospy.preprocessing.provenance.models import PreprocessingStageProvenance
from phospy.provenance.models.environment import EnvironmentProvenance
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.serialization.tables import table_fingerprint_to_payload
from phospy.run_lineage.models import RunProvenance
from phospy.run_lineage.serialization import (
    from_payload as run_provenance_from_payload,
    to_payload as run_provenance_to_payload,
)
from phospy.preprocessing.batch_correction.execution import (
    RuvIIIReplicateStructure,
    run_ruv_iii,
)
from phospy.preprocessing.batch_correction.environment import (
    BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES,
    collect_batch_correction_environment_provenance,
)
from phospy.preprocessing.batch_correction.execution.linear import (
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.preprocessing.batch_correction.execution.ruv_iii_executor import (
    _prepare_matrix as prepare_ruv_iii_missingness,
)
from phospy.preprocessing.control_sites.contracts import (
    ControlSiteSet as OwnerControlSiteSet,
)
from phospy.preprocessing.control_sites.discovery_contracts import (
    SpsDiscoveryRequest as OwnerSpsDiscoveryRequest,
    SpsReferenceDataset as OwnerSpsReferenceDataset,
)
from phospy.preprocessing.control_sites.discovery_results import (
    SpsDiscoveryResult as OwnerSpsDiscoveryResult,
)
from phospy.preprocessing.control_sites.workflow import (
    SpsDiscoveryWorkflow as OwnerSpsDiscoveryWorkflow,
)
from phospy.provenance.serialization.batch_correction import (
    batch_correction_provenance_from_payload,
    batch_correction_provenance_to_payload,
)
from phospy.preprocessing.batch_correction import (
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
)
from phospy.preprocessing.provenance.serialization import (
    stage_from_payload,
    stage_to_payload,
)


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
    retired_route_report = _verify_retired_routes()
    prediction_identity_report = _verify_prediction_identities_and_import_orders()
    reference_builder_report = _exercise_reference_bundle_builder()
    import_order_report = _verify_references_evidence_import_orders(
        repo_root=repo_root,
        environment_root=environment_root,
    )
    import_order_report.extend(
        _verify_enrichment_import_orders(
            repo_root=repo_root,
            environment_root=environment_root,
        )
    )
    import_order_report.extend(
        _verify_control_site_import_orders(
            repo_root=repo_root,
            environment_root=environment_root,
        )
    )
    import_order_report.extend(
        _verify_batch_correction_boundary_import_orders(
            repo_root=repo_root,
            environment_root=environment_root,
        )
    )
    import_order_report.extend(
        _verify_protein_aware_boundary_import_orders(
            repo_root=repo_root,
            environment_root=environment_root,
        )
    )
    resource_report = _verify_bundled_rat_reference_resources()
    scientific_report = _exercise_public_scientific_contracts()
    enrichment_report = _exercise_enrichment_capability()

    print(
        json.dumps(
            {
                "status": "ok",
                "artifact_kind": artifact_kind,
                "python": sys.version.split()[0],
                "executable": str(pathlib.Path(sys.executable).resolve()),
                "phospy_file": str(package_file),
                "runtime_version": str(phospy.__version__),
                "environment_versions": {
                    "python": sys.version.split()[0],
                    "phospy": str(phospy.__version__),
                    "numpy": metadata.version("numpy"),
                    "pandas": metadata.version("pandas"),
                    "scipy": metadata.version("scipy"),
                    "scikit-learn": metadata.version("scikit-learn"),
                    "pip": metadata.version("pip"),
                },
                "runtime_platform": {
                    "system": platform.system(),
                    "release": platform.release(),
                    "machine": platform.machine(),
                    "platform": platform.platform(),
                },
                "dependency_identities": sorted(
                    {
                        str(distribution.metadata["Name"]): distribution.version
                        for distribution in metadata.distributions()
                        if distribution.metadata["Name"]
                    }.items()
                ),
                "public_surface_report": public_surface_report,
                "retired_route_report": retired_route_report,
                "prediction_identity_report": prediction_identity_report,
                "reference_builder_report": reference_builder_report,
                "import_order_report": import_order_report,
                "resource_report": resource_report,
                "scientific_report": scientific_report,
                "enrichment_report": enrichment_report,
            },
            sort_keys=True,
        )
    )


def _verify_public_surface() -> dict[str, object]:
    expected_counts = {"phospy": 3, "phospy.api": 66, "phospy.advanced": 124}
    observed_counts = {
        "phospy": len(phospy.__all__),
        "phospy.api": len(api.__all__),
        "phospy.advanced": len(advanced.__all__),
    }
    if observed_counts != expected_counts:
        raise AssertionError(
            f"aggregate export counts changed: {observed_counts!r}"
        )
    expected_inventories = {
        "phospy": api_inventory.ROOT_PUBLIC_API,
        "phospy.api": api_inventory.ORDINARY_PUBLIC_API,
        "phospy.advanced": api_inventory.ADVANCED_PUBLIC_API,
    }
    observed_inventories = {
        "phospy": tuple(phospy.__all__),
        "phospy.api": tuple(api.__all__),
        "phospy.advanced": tuple(advanced.__all__),
    }
    expected_inventory_sha256 = {
        "phospy": "e1ba016281f8a1d46bf4c2cabe61525c4d3e2b7d3fed53ea833daeeec9c02120",
        "phospy.api": "75988b66501a85463a19df9a9c16a6048906123b19f727be1e79828ef96783be",
        "phospy.advanced": "9944ed70bbfa63c18a9e49768b791b641e41dd8fb84531dcc226aa481482421d",
    }
    inventory_sha256 = {
        name: hashlib.sha256(
            json.dumps(list(names), separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        for name, names in observed_inventories.items()
    }
    if observed_inventories != expected_inventories:
        raise AssertionError("ordered aggregate exports do not match inventories")
    if inventory_sha256 != expected_inventory_sha256:
        raise AssertionError(
            f"ordered aggregate export baseline changed: {inventory_sha256!r}"
        )
    for name in REQUIRED_API_NAMES:
        if name not in api.__all__:
            raise AssertionError(f"phospy.api.__all__ is missing {name!r}")
        getattr(api, name)
    for name in (
        "DatasetIntensityTransformConfig",
        "DatasetProteinAwarePreparationConfig",
        "DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1",
        "DifferentialProteinAwareModelConfig",
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
    if differential_policies.PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION != (
        "duplicate_correlation"
    ):
        raise AssertionError("Differential duplicate-correlation policy changed")
    if PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION != "duplicate_correlation":
        raise AssertionError("advanced duplicate-correlation policy export changed")
    if (
        "duplicate_correlation"
        not in differential_policies.SUPPORTED_PAIRED_DESIGN_POLICIES
    ):
        raise AssertionError("supported paired-design policies omit duplicate_correlation")
    if phospy.AnalysisReadyDatasetBuilder is not api.AnalysisReadyDatasetBuilder:
        raise AssertionError("root AnalysisReadyDatasetBuilder is not api export")
    for removed_root_name in (
        "AnalysisReadyPhosphoDataset",
        "DifferentialAnalysisWorkflow",
    ):
        if hasattr(phospy, removed_root_name):
            raise AssertionError(f"removed root name survived: {removed_root_name}")
    if phospy.KinaseWorkflow is not api.KinaseWorkflow:
        raise AssertionError("root KinaseWorkflow is not api export")
    if api.ReferenceBundleBuilder is not OwnerReferenceBundleBuilder:
        raise AssertionError("ReferenceBundleBuilder is not the References owner")
    _verify_direct_aggregate_identities()
    foreign_request_definitions = {
        "PhosphositeImportRequest": PhosphositeImportRequest,
        "PhosphositeImporter": PhosphositeImporter,
    }
    for name, definition in foreign_request_definitions.items():
        if definition.__module__ != "phospy.io.importers.contracts":
            raise AssertionError(
                f"importer request definition owner changed: {name}"
            )
        if getattr(api, name) is not definition:
            raise AssertionError(
                f"foreign request API identity changed: {name}"
            )
    if SignalomeWorkflowRequest.__module__ != "phospy.signalome.requests":
        raise AssertionError("Signalome workflow request owner changed")
    if api.SignalomeWorkflowRequest is not SignalomeWorkflowRequest:
        raise AssertionError("Signalome workflow request API identity changed")
    if KinaseWorkflowRequest.__module__ != "phospy.kinase.requests":
        raise AssertionError("Kinase workflow request owner changed")
    for removed_module in (
        "phospy._deprecations",
        "phospy.api._compat",
        "phospy.api.requests",
        "phospy.api.results",
        "phospy.api.workflows",
        "phospy.advanced.configs",
        "phospy.advanced.results",
        "phospy.contracts.configs.enrichment",
        "phospy.contracts.enrichment_identifier_sets",
        "phospy.contracts.results.enrichment",
        "phospy.science.configs.enrichment",
        "phospy.science.enrichment.collections",
        "phospy.science.enrichment.models",
        "phospy.science.enrichment.ora",
        "phospy.validation.workflows.enrichment",
        "phospy.workflows.enrichment.caveats",
        "phospy.workflows.enrichment.executor",
        "phospy.workflows.enrichment.interpreter",
        "phospy.workflows.enrichment.models",
        "phospy.workflows.enrichment.ora_execution",
        "phospy.workflows.enrichment.provenance",
        "phospy.workflows.enrichment.public",
        "phospy.workflows.enrichment.result_assembly",
        "phospy.workflows.enrichment.set_filtering",
        "phospy.workflows.enrichment.validator",
        "phospy.workflows.enrichment.workflow",
        "phospy.tables.activity",
        "phospy.contracts.configs.preprocessing.control_sites",
        "phospy.science.configs.preprocessing.control_sites",
        "phospy.science.batch_correction",
        "phospy.science.datasets.preprocessing.control_sites",
        "phospy.science.datasets.preprocessing.protein_aware_models",
        "phospy.contracts.results.preprocessing",
        "phospy.science.batch_correction.sps_discovery",
        "phospy.validation.workflows.batch_correction.control_sites",
        "phospy.validation.workflows.batch_correction.control_site_workflow",
        "phospy.validation.workflows.batch_correction.sps_discovery",
        "phospy.workflows.batch_correction.sps_discovery",
    ):
        try:
            removed_spec = importlib.util.find_spec(removed_module)
        except ModuleNotFoundError:
            removed_spec = None
        if removed_spec is not None:
            raise AssertionError(f"removed compatibility module survived: {removed_module}")
    protein_aware_report = _exercise_protein_aware_configuration_contract()
    ticket_1_report = _verify_withdrawn_peptide_to_site_boundary()
    return {
        "required_api_names": list(REQUIRED_API_NAMES),
        "aggregate_export_counts": observed_counts,
        "aggregate_ordered_names": {
            name: list(names) for name, names in observed_inventories.items()
        },
        "aggregate_inventory_sha256": inventory_sha256,
        "direct_identity_count": sum(observed_counts.values()),
        "foreign_request_definitions": {
            name: definition.__module__
            for name, definition in foreign_request_definitions.items()
        },
        "protein_aware_configuration": protein_aware_report,
        "ticket_1_posthoc_peptide_to_site_boundary": ticket_1_report,
    }


def _verify_retired_routes() -> dict[str, object]:
    retired_modules = __RETIRED_CAPABILITY_ROUTES__ + __RETIRED_HORIZONTAL_ROUTES__
    source = __RETIRED_ROUTE_IMPORT_PROBE_SOURCE__
    route_results = []
    for module_name in retired_modules:
        completed = subprocess.run(
            (sys.executable, "-I", "-c", source, module_name),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise AssertionError(
                "retired compatibility route absence was not established: "
                f"{module_name}; returncode={completed.returncode}; "
                f"stdout={completed.stdout!r}; stderr={completed.stderr!r}"
            )
        try:
            route_result = json.loads(completed.stdout)
        except json.JSONDecodeError as error:
            raise AssertionError(
                f"retired route probe emitted invalid JSON for {module_name}: "
                f"{completed.stdout!r}"
            ) from error
        if not isinstance(route_result, dict):
            raise AssertionError(
                f"retired route probe result is not an object: {route_result!r}"
            )
        if route_result != {
            "exception_type": "ModuleNotFoundError",
            "module_not_found_error_name": route_result.get(
                "module_not_found_error_name"
            ),
            "requested_route_or_absent_ancestor": True,
            "route": module_name,
            "status": "absent",
        }:
            raise AssertionError(
                f"retired route probe did not prove route absence: {route_result!r}"
            )
        missing_name = route_result["module_not_found_error_name"]
        if not isinstance(missing_name, str) or not (
            missing_name == module_name
            or module_name.startswith(missing_name + ".")
        ):
            raise AssertionError(
                f"retired route probe reported unrelated missing module: "
                f"{route_result!r}"
            )
        route_results.append(route_result)
    return {
        "status": "ok",
        "fresh_process_failure": True,
        "retired_module_count": len(route_results),
        "retired_modules": list(retired_modules),
        "route_results": route_results,
    }


def _verify_prediction_identities_and_import_orders() -> dict[str, object]:
    identities = {
        "api.KinasePredictionResult": (
            api.KinasePredictionResult,
            OwnerKinasePredictionResult,
        ),
        "advanced.KinasePredictionConfig": (
            advanced.KinasePredictionConfig,
            OwnerKinasePredictionConfig,
        ),
        "advanced.KinasePredictionMode": (
            advanced.KinasePredictionMode,
            OwnerKinasePredictionMode,
        ),
        "advanced.KinaseAdaptivePolicy": (
            advanced.KinaseAdaptivePolicy,
            OwnerKinaseAdaptivePolicy,
        ),
    }
    changed = [name for name, (facade, owner) in identities.items() if facade is not owner]
    if changed:
        raise AssertionError(f"installed prediction identities changed: {changed!r}")

    code = r'''\
import importlib
import json
import sys

retired = json.loads(sys.argv[1])
order = sys.argv[2]

def supported():
    import phospy.api as api
    import phospy.advanced as advanced
    from phospy.kinase.prediction.config import KinaseAdaptivePolicy
    from phospy.kinase.prediction.config import KinasePredictionConfig
    from phospy.kinase.prediction.config import KinasePredictionMode
    from phospy.kinase.prediction.results import KinasePredictionResult
    assert api.KinasePredictionResult is KinasePredictionResult
    assert advanced.KinasePredictionConfig is KinasePredictionConfig
    assert advanced.KinasePredictionMode is KinasePredictionMode
    assert advanced.KinaseAdaptivePolicy is KinaseAdaptivePolicy

def absent(route):
    try:
        importlib.import_module(route)
    except ModuleNotFoundError as error:
        if not isinstance(error.name, str) or not (
            error.name == route or route.startswith(error.name + ".")
        ):
            raise
    else:
        raise AssertionError(f"retired route imported: {route}")

if order == "supported-first":
    supported()
    for route in retired:
        absent(route)
elif order == "retired-first":
    for route in retired:
        absent(route)
    supported()
else:
    for route in retired:
        supported()
        absent(route)
print(json.dumps({"order": order, "status": "ok"}, sort_keys=True))
'''
    retired = list(__RETIRED_CAPABILITY_ROUTES__ + __RETIRED_HORIZONTAL_ROUTES__)
    results = []
    for order in ("supported-first", "retired-first", "alternating"):
        completed = subprocess.run(
            (
                sys.executable,
                "-I",
                "-c",
                code,
                json.dumps(retired),
                order,
            ),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode:
            raise AssertionError(
                f"installed prediction import-order check failed: {order}; "
                f"stdout={completed.stdout!r}; stderr={completed.stderr!r}"
            )
        result = json.loads(completed.stdout)
        if result != {"order": order, "status": "ok"}:
            raise AssertionError(f"malformed prediction import-order result: {result!r}")
        results.append(result)
    return {
        "status": "ok",
        "exact_identity_count": len(identities),
        "identities": sorted(identities),
        "import_orders": results,
    }


def _verify_direct_aggregate_identities() -> None:
    value_owner_overrides = {
        "DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1": (
            "phospy.differential.policies"
        ),
        "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION": (
            "phospy.differential.policies"
        ),
        "PAIRED_DESIGN_POLICY_FIXED_BLOCK": "phospy.differential.policies",
        "PAIRED_DESIGN_POLICY_REJECT": "phospy.differential.policies",
        "DifferentialImputedValuePolicy": "phospy.differential.policies",
        "DifferentialReliabilityProfile": "phospy.differential.policies",
        "MultipleTestingMethod": "phospy.differential.policies",
        "PairedDesignPolicy": "phospy.differential.policies",
        "TechnicalReplicatePolicy": "phospy.differential.policies",
        "EnrichmentIdentifierKind": "phospy.enrichment.contracts.vocabulary",
        "EnrichmentMethod": "phospy.enrichment.contracts.vocabulary",
        "EnrichmentOutsideBackgroundPolicy": (
            "phospy.enrichment.contracts.vocabulary"
        ),
        "MultipleTestingCorrection": "phospy.statistics.multiple_testing",
        "DatasetBatchCorrectionMethod": (
            "phospy.preprocessing.batch_correction.contracts.configuration"
        ),
        "DatasetPreprocessingBatchCorrectionConfig": (
            "phospy.preprocessing.batch_correction.contracts.configuration"
        ),
        "SpsRuvBatchCorrectionMethod": (
            "phospy.preprocessing.batch_correction.contracts.configuration"
        ),
        "SignalomeAssignmentPolicy": "phospy.signalome.config",
        "SignalomeCandidateScoringPolicy": "phospy.signalome.config",
        "SignalomeKinaseNetworkPolicy": "phospy.signalome.config",
        "SignalomeMode": "phospy.signalome.config",
        "SignalomeScorePreconditioningPolicy": "phospy.signalome.config",
    }
    for module in (phospy, api, advanced):
        for name in module.__all__:
            value = getattr(module, name)
            config_owner = next(
                (
                    owner
                    for owner in (
                        evidence_config,
                        kinase_activity_config,
                        kinase_prediction_config,
                        kinase_scoring_config,
                        batch_config,
                        missingness_config,
                        preprocessing_config,
                        references_config,
                        signalome_config,
                    )
                    if hasattr(owner, name) and value is getattr(owner, name)
                ),
                None,
            )
            if config_owner is None:
                owner_module_name = value_owner_overrides.get(name)
                owner_qualname = name
                if owner_module_name is None:
                    owner_module_name = getattr(value, "__module__", None)
                    owner_qualname = getattr(value, "__qualname__", name)
                if owner_module_name is None:
                    raise AssertionError(f"aggregate export has no owner module: {name}")
                owner_value = importlib.import_module(owner_module_name)
                for part in owner_qualname.split("."):
                    owner_value = getattr(owner_value, part)
            else:
                owner_value = getattr(config_owner, name)
            if value is not owner_value:
                raise AssertionError(
                    f"aggregate export is not its direct owner identity: {module.__name__}.{name}"
                )


def _exercise_reference_bundle_builder() -> dict[str, object]:
    if api.ReferenceBundleBuilder is not OwnerReferenceBundleBuilder:
        raise AssertionError("public ReferenceBundleBuilder owner identity changed")
    if tuple(inspect.signature(api.ReferenceBundleBuilder).parameters) != ():
        raise AssertionError("ReferenceBundleBuilder constructor signature changed")
    builder = api.ReferenceBundleBuilder()
    if not isinstance(builder, OwnerReferenceBundleBuilder):
        raise AssertionError("zero-argument ReferenceBundleBuilder construction failed")

    with tempfile.TemporaryDirectory(prefix="phospy-installed-reference-builder-") as tmp_dir:
        source_root = pathlib.Path(tmp_dir)
        kinase_path = source_root / "kinase_substrates.csv"
        sequence_path = source_root / "site_sequences.tsv"
        pd.DataFrame(
            {
                "kinase": ["akt1", "map2k1"],
                "site_id": ["MAPK1;S123;", "MAPK1;Y185;"],
                "organism": ["human", "human"],
            }
        ).to_csv(kinase_path, index=False)
        pd.DataFrame(
            {
                "site_id": ["MAPK1;S123;", "MAPK1;Y185;"],
                "site_sequence": ["A" * 15 + "S" + "A" * 15, "A" * 15 + "Y" + "A" * 15],
                "gene_symbol": ["MAPK1", "MAPK1"],
                "protein_accession": ["P28482", "P28482"],
                "organism": ["human", "human"],
            }
        ).to_csv(sequence_path, index=False, sep="\t")

        bundle = builder.run(
            api.ReferenceBundleBuildRequest(
                organism=api.Organism.HUMAN,
                kinase_substrate_path=kinase_path,
                site_sequence_path=sequence_path,
                source_name="installed synthetic local reference",
                source_version="installed-probe-v1",
                retrieved_at="2026-10-01",
                license="synthetic test license",
                redistribution_status="synthetic local probe; not for redistribution",
                identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
                bundle_id="installed_synthetic_human_v1",
            )
        )
        if bundle.kinase_substrate_map.to_dict(orient="records") != [
            {"kinase": "AKT1", "substrate_site": "MAPK1;S123;"},
            {"kinase": "MAP2K1", "substrate_site": "MAPK1;Y185;"},
        ]:
            raise AssertionError("installed reference builder rows changed")
        expected_site_sequences = {
            "MAPK1;S123;": {
                "site_sequence": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "gene_symbol": "MAPK1",
                "protein_accession": "P28482",
                "display_id": "MAPK1;S123;",
            },
            "MAPK1;Y185;": {
                "site_sequence": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "gene_symbol": "MAPK1",
                "protein_accession": "P28482",
                "display_id": "MAPK1;Y185;",
            },
        }
        if bundle.site_sequences.to_dict(orient="index") != expected_site_sequences:
            raise AssertionError("installed reference builder sequence rows changed")
        if bundle.manifest is None or bundle.provenance is None:
            raise AssertionError("installed reference builder omitted metadata")
        manifest = bundle.manifest
        expected_source_sha256 = {
            "kinase_substrate": hashlib.sha256(kinase_path.read_bytes()).hexdigest(),
            "site_sequences": hashlib.sha256(sequence_path.read_bytes()).hexdigest(),
        }
        reference_version_payload = (
            f"kinase_substrate:{expected_source_sha256['kinase_substrate']}\n"
            f"site_sequences:{expected_source_sha256['site_sequences']}\n"
        )
        expected_reference_version = "local-snapshot-sha256-" + hashlib.sha256(
            reference_version_payload.encode("ascii")
        ).hexdigest()
        observed_manifest_metadata = {
            "reference_id": manifest.reference_id,
            "display_name": manifest.display_name,
            "organism": manifest.organism,
            "taxonomy_id": manifest.taxonomy_id,
            "protein_namespace": manifest.protein_namespace,
            "reference_version": manifest.reference_version,
            "table_sha256": manifest.table_sha256,
            "source_name": manifest.source_name,
            "source_version": manifest.source_version,
            "retrieved_at": manifest.retrieved_at.isoformat(),
            "license_name": manifest.license_name,
            "redistribution_status": manifest.redistribution_status.value,
            "redistribution_notes": manifest.redistribution_notes,
            "generated_by": manifest.generated_by,
            "generated_at_utc": manifest.generated_at_utc,
            "manifest_schema_version": manifest.manifest_schema_version,
            "sequence_context_policy": manifest.sequence_context_policy,
            "sequence_window_length": manifest.sequence_window_length,
            "sequence_center_index": manifest.sequence_center_index,
            "allowed_sequence_alphabet": manifest.allowed_sequence_alphabet,
            "organism_common_name": manifest.organism_common_name,
            "supports": manifest.supports,
            "limitations": manifest.limitations,
        }
        expected_manifest_metadata = {
            "reference_id": "installed_synthetic_human_v1",
            "display_name": "installed synthetic local reference",
            "organism": "human",
            "taxonomy_id": 9606,
            "protein_namespace": "display_id (GENE_SYMBOL;RESIDUE;)",
            "reference_version": expected_reference_version,
            "table_sha256": expected_source_sha256["kinase_substrate"],
            "source_name": "installed synthetic local reference",
            "source_version": "installed-probe-v1",
            "retrieved_at": "2026-10-01",
            "license_name": "synthetic test license",
            "redistribution_status": "unresolved",
            "redistribution_notes": (
                "synthetic local probe; not for redistribution"
            ),
            "generated_by": "ReferenceBundleBuilder",
            "generated_at_utc": "2026-10-01T00:00:00Z",
            "manifest_schema_version": "1.1",
            "sequence_context_policy": "centered phosphosite sequence window",
            "sequence_window_length": 31,
            "sequence_center_index": 15,
            "allowed_sequence_alphabet": "ACDEFGHIKLMNPQRSTVWY",
            "organism_common_name": "human",
            "supports": ("kinase_workflow", "site_sequence_derivation"),
            "limitations": (
                "caller-supplied local source files; redistribution governed by request metadata",
            ),
        }
        if observed_manifest_metadata != expected_manifest_metadata:
            raise AssertionError("installed reference builder manifest metadata changed")
        source_files = bundle.manifest.source_files
        observed_source_sha256 = {
            role: str(source_files[role]["sha256"])
            for role in expected_source_sha256
        }
        if observed_source_sha256 != expected_source_sha256:
            raise AssertionError("installed reference builder source hashes changed")
        observed_source_file_metadata = {
            role: {
                "role": source_files[role]["role"],
                "format": source_files[role]["format"],
                "row_count": source_files[role]["row_count"],
                "column_names": source_files[role]["column_names"],
            }
            for role in expected_source_sha256
        }
        expected_source_file_metadata = {
            "kinase_substrate": {
                "role": "kinase_substrate",
                "format": "csv",
                "row_count": 2,
                "column_names": ["kinase", "site_id", "organism"],
            },
            "site_sequences": {
                "role": "site_sequences",
                "format": "tsv",
                "row_count": 2,
                "column_names": [
                    "site_id",
                    "site_sequence",
                    "gene_symbol",
                    "protein_accession",
                    "organism",
                ],
            },
        }
        if observed_source_file_metadata != expected_source_file_metadata:
            raise AssertionError("installed reference builder source metadata changed")
        provenance = bundle.provenance
        observed_provenance_metadata = {
            "source_type": provenance.source_type,
            "organism": provenance.organism.value,
            "bundle_id": provenance.bundle_id,
            "source_name": provenance.source_name,
            "source_version": provenance.source_version,
            "retrieved_at": provenance.retrieved_at,
            "identifier_namespace": provenance.identifier_namespace,
            "sequence_window": dict(provenance.sequence_window or {}),
        }
        expected_provenance_metadata = {
            "source_type": "local",
            "organism": "human",
            "bundle_id": "installed_synthetic_human_v1",
            "source_name": "installed synthetic local reference",
            "source_version": "installed-probe-v1",
            "retrieved_at": "2026-10-01",
            "identifier_namespace": "display_id (GENE_SYMBOL;RESIDUE;)",
            "sequence_window": {
                "upstream_residues": 15,
                "downstream_residues": 15,
                "central_residue_required": True,
            },
        }
        if observed_provenance_metadata != expected_provenance_metadata:
            raise AssertionError("installed reference builder provenance metadata changed")
        reference_payload = reference_to_payload(provenance)
        restored_reference = reference_from_payload(reference_payload)
        if reference_to_payload(restored_reference) != reference_payload:
            raise AssertionError("installed direct Reference provenance round trip changed")
        if restored_reference.__class__ is not provenance.__class__:
            raise AssertionError("installed direct Reference provenance identity changed")
        table_fingerprints = {
            item.name: {
                "rows": item.rows,
                "columns": item.columns,
                "index_name": item.index_name,
                "column_names": item.column_names,
                "dtypes": item.dtypes,
                "exact_hash_algorithm": item.exact_hash_algorithm,
                "exact_hash_value": item.exact_hash_value,
                "tolerance_hash_algorithm": item.tolerance_hash_algorithm,
                "tolerance_hash_value": item.tolerance_hash_value,
            }
            for item in bundle.provenance.table_fingerprints
        }
        expected_table_fingerprints = {
            "references.kinase_substrate_map": {
                "rows": 2,
                "columns": 2,
                "index_name": None,
                "column_names": ("kinase", "substrate_site"),
                "dtypes": ("string", "string"),
                "exact_hash_algorithm": "sha256-stable-json-v1",
                "exact_hash_value": (
                    "a490fb6499d5c9deb87a7ce28284db5e916be44527efc9657827274520e696c6"
                ),
                "tolerance_hash_algorithm": "sha256-float-round-8dp-v1",
                "tolerance_hash_value": (
                    "a490fb6499d5c9deb87a7ce28284db5e916be44527efc9657827274520e696c6"
                ),
            },
            "references.site_sequences": {
                "rows": 2,
                "columns": 4,
                "index_name": "site_id",
                "column_names": (
                    "site_sequence",
                    "gene_symbol",
                    "protein_accession",
                    "display_id",
                ),
                "dtypes": ("string", "string", "string", "string"),
                "exact_hash_algorithm": "sha256-stable-json-v1",
                "exact_hash_value": (
                    "d6ca2a89d15e7a9882e23e7e2ff55325b84e818125b7ae903a218f5fd88e021a"
                ),
                "tolerance_hash_algorithm": "sha256-float-round-8dp-v1",
                "tolerance_hash_value": (
                    "d6ca2a89d15e7a9882e23e7e2ff55325b84e818125b7ae903a218f5fd88e021a"
                ),
            },
        }
        if table_fingerprints != expected_table_fingerprints:
            raise AssertionError("installed reference builder table fingerprints changed")

        return {
            "constructor_signature": str(inspect.signature(api.ReferenceBundleBuilder)),
            "owner_identity": True,
            "kinase_substrate_rows": int(bundle.kinase_substrate_map.shape[0]),
            "site_sequence_rows": int(bundle.site_sequences.shape[0]),
            "reference_id": manifest.reference_id,
            "source_name": manifest.source_name,
            "manifest_metadata": observed_manifest_metadata,
            "provenance_metadata": observed_provenance_metadata,
            "provenance_payload_sha256": _probe_payload_sha256(reference_payload),
            "source_sha256": observed_source_sha256,
            "source_file_metadata": observed_source_file_metadata,
            "table_fingerprints": table_fingerprints,
        }


def _probe_json_value(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _probe_json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [_probe_json_value(item) for item in value]
    if isinstance(value, Enum):
        return _probe_json_value(value.value)
    if value is None or isinstance(value, str | int | float | bool):
        return value
    raise TypeError(f"unsupported enrichment probe value: {type(value).__name__}")


def _probe_payload_sha256(value: object) -> str:
    payload = json.dumps(
        _probe_json_value(value),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _probe_frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "index": [str(value) for value in frame.index],
        "columns": [str(value) for value in frame.columns],
        "values": frame.to_numpy().tolist(),
    }


def _fixed_enrichment_environment() -> EnvironmentProvenance:
    return EnvironmentProvenance(
        package_name="phospy",
        package_version="1.7.6",
        python_version="3.11.15",
        dependency_versions={
            "numpy": "2.3.5",
            "pandas": "2.3.3",
            "scipy": "1.16.3",
        },
        platform={
            "implementation": "CPython",
            "system": "Windows",
            "machine": "AMD64",
        },
        blas_lapack={"provider": "controlled-enrichment-regression"},
        thread_environment={
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
        },
        timezone="UTC",
        locale={"preferred_encoding": "UTF-8"},
        constraints_fingerprint={
            "algorithm": "sha256",
            "value": "enrichment-regression-constraints",
        },
    )


def _exercise_enrichment_capability() -> dict[str, object]:
    from phospy.enrichment.contracts.collections import (
        GeneSetCollection as OwnerGeneSetCollection,
    )
    from phospy.enrichment.contracts.collections import (
        PtmSetCollection as OwnerPtmSetCollection,
    )
    from phospy.enrichment.contracts.config import EnrichmentConfig as OwnerConfig
    from phospy.enrichment.contracts.requests import (
        EnrichmentWorkflowRequest as OwnerRequest,
    )
    from phospy.enrichment.contracts.results import (
        EnrichmentWorkflowResult as OwnerResult,
    )
    from phospy.enrichment.workflow import EnrichmentWorkflow as OwnerWorkflow

    expected_request_parameters = (
        "identifier_column",
        "identifier_kind",
        "set_collection",
        "background_universe",
        "config",
        "input_table",
        "selected_identifiers",
        "selected_identifier_provenance",
        "background_identifier_provenance",
    )
    if tuple(inspect.signature(api.EnrichmentWorkflow).parameters) != ():
        raise AssertionError("installed EnrichmentWorkflow constructor changed")
    if tuple(inspect.signature(api.EnrichmentWorkflowRequest).parameters) != (
        expected_request_parameters
    ):
        raise AssertionError("installed EnrichmentWorkflowRequest field order changed")
    identities = {
        "EnrichmentConfig": api.EnrichmentConfig is OwnerConfig,
        "EnrichmentWorkflowRequest": api.EnrichmentWorkflowRequest is OwnerRequest,
        "EnrichmentWorkflowResult": api.EnrichmentWorkflowResult is OwnerResult,
        "GeneSetCollection": api.GeneSetCollection is OwnerGeneSetCollection,
        "PtmSetCollection": api.PtmSetCollection is OwnerPtmSetCollection,
        "EnrichmentWorkflow": api.EnrichmentWorkflow is OwnerWorkflow,
    }
    if not all(identities.values()):
        raise AssertionError(f"installed Enrichment owner identities changed: {identities}")

    def identifier_provenance(source_label: str, count: int) -> object:
        return api.EnrichmentIdentifierSetProvenance(
            source_type="raw_identifier_list",
            source_label=source_label,
            identifier_count=count,
            upstream_workflow_id="enrichment-regression-control",
            upstream_result_id=f"{source_label}-v1",
        )

    gene_request = api.EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind="gene_symbol",
        set_collection=api.GeneSetCollection(
            sets={
                "SIGNALING": ("AKT1", "MAPK1", "MTOR"),
                "CELL_CYCLE": ("CDK1", "CDK2", "MAPK1"),
                "NO_BACKGROUND_MEMBERS": ("OUTSIDE_A", "OUTSIDE_B"),
                "TOO_LARGE": ("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
            },
            identifier_kind="gene_symbol",
            term_names={
                "SIGNALING": "Signaling",
                "CELL_CYCLE": "Cell cycle",
                "NO_BACKGROUND_MEMBERS": "No background members",
                "TOO_LARGE": "Too large",
            },
            source_name="synthetic_gene_sets",
            source_version="regression-v1",
        ),
        selected_identifiers=("AKT1", "UNKNOWN", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=api.EnrichmentConfig(
            min_set_size=1,
            max_set_size=3,
            selected_outside_background_policy="drop",
        ),
        selected_identifier_provenance=identifier_provenance("gene-selected", 3),
        background_identifier_provenance=identifier_provenance(
            "gene-background", 5
        ),
    )
    ptm_request = api.EnrichmentWorkflowRequest(
        identifier_column="site_key",
        identifier_kind="site_key",
        set_collection=api.PtmSetCollection(
            sets={
                "MOTIF_A": ("rat|P12345|S10", "rat|P12345|T20"),
                "MOTIF_B": ("rat|P12345|Y30", "rat|P99999|S5"),
                "MOTIF_EMPTY": ("rat|OUTSIDE|S1",),
            },
            identifier_kind="site_key",
            source_name="synthetic_ptm_sets",
            source_version="regression-v1",
        ),
        input_table=pd.DataFrame(
            {"site_key": ["rat|P12345|S10", " rat|P12345|T20 "]}
        ),
        background_universe=(
            "rat|P12345|S10",
            "rat|P12345|T20",
            "rat|P12345|Y30",
            "rat|P99999|S5",
        ),
        selected_identifier_provenance=identifier_provenance("ptm-selected", 2),
        background_identifier_provenance=identifier_provenance("ptm-background", 4),
    )
    with patch(
        "phospy.enrichment.provenance.collect_environment_provenance",
        return_value=_fixed_enrichment_environment(),
    ):
        workflow = api.EnrichmentWorkflow()
        results = {
            "gene": workflow.run(gene_request),
            "ptm": workflow.run(ptm_request),
        }
    expected_hashes = {
        "gene": {
            "records": "47c2a18f7d21985b9dddc8e48502d47766021b5a74a636fea8656bdece75836d",
            "diagnostics": "a7cbbe4e0bf00f2fc4b9b4050fa52973acba4ad00afd3d7df9fb42b0b33c1e45",
            "caveats": "dec14174e1fa041249436e401ca06bd49255ddbd05aec5cef50690468f1169c6",
            "selected_identifier_provenance": "24c8f4124e6645f1e0e3ff8c0bf757db6a9261108404b859720fbb548f916da0",
            "background_identifier_provenance": "dc2f9a9d8ecb161756fdf0d31dc98f78d469a676003d20e4a501cd2debedb014",
            "workflow_parameters": "d75f8f9a4aafd9b3625e39f02b626d472891569757d584acf73b2da91a70e48d",
            "input_tables": "76f506e6c9e661f0724e97aa022f83349f3b5f74a823d67177d5bb3f06640b4e",
            "output_tables": "ccc7cd07bd3a37df872906a6664cf8e820a5e24609c875baa6adaba59445b6d9",
            "run_provenance": "1a20b51c3ae39867e2d005c12dafde31077dc89f88e290b5a4bba40b23e0727d",
        },
        "ptm": {
            "records": "b45bbea8c8c1bcea1c23b93cd2c4e12e51d38e6f88d696ae52ab796e32cae5f6",
            "diagnostics": "06db27c78a4f4c5b5c32e669eb3e0aaf872ed81f0cd839bbb0d6cd74fb8ef53b",
            "caveats": "a09d3603bc0000479d8e93e35c5ad833b761ca04804b515dde19629913e3c265",
            "selected_identifier_provenance": "d72e0d5ab24c9ab01460b3155c471cf7e91f4ae583b86586c484d6556c15e607",
            "background_identifier_provenance": "cd069d4d03cc663a0ff6c38590e0072a4d30ecec0cdac5e6154b3132b4f7f71e",
            "workflow_parameters": "2bd9faf6c6a2bcf7390888ec0279208bc78fba5350aed8c8685ffa489ff31cc4",
            "input_tables": "025c4fd1fc741ed06e2e8d44c54eef7f65741fd405fe87a69a197d796e1ee677",
            "output_tables": "1e68f999766e008ce26433937343a064f0a88e9ce5f49f14958cf53a6540788f",
            "run_provenance": "d34915364a6b40ef0d80a59bddd1acaeace302c314f15a33246dbc4a7a2cee56",
        },
    }
    reports: dict[str, object] = {}
    for name, result in results.items():
        if result.provenance is None:
            raise AssertionError(f"installed {name} Enrichment omitted provenance")
        record_payloads = [
            {field.name: getattr(record, field.name) for field in fields(record)}
            for record in result.records
        ]
        caveat_payloads = [caveat.to_payload() for caveat in result.caveats]
        selected_payload = result.selected_identifier_provenance.to_payload()
        background_payload = result.background_identifier_provenance.to_payload()
        input_table_payloads = [
            table_fingerprint_to_payload(item)
            for item in result.provenance.input_tables
        ]
        output_table_payloads = [
            table_fingerprint_to_payload(item)
            for item in result.provenance.output_tables
        ]
        run_provenance_payload = run_provenance_to_payload(result.provenance)
        observed_hashes = {
            "records": _probe_payload_sha256(record_payloads),
            "diagnostics": _probe_payload_sha256(result.diagnostics),
            "caveats": _probe_payload_sha256(caveat_payloads),
            "selected_identifier_provenance": _probe_payload_sha256(
                selected_payload
            ),
            "background_identifier_provenance": _probe_payload_sha256(
                background_payload
            ),
            "workflow_parameters": _probe_payload_sha256(
                result.provenance.workflow_parameters
            ),
            "input_tables": _probe_payload_sha256(input_table_payloads),
            "output_tables": _probe_payload_sha256(output_table_payloads),
            "run_provenance": _probe_payload_sha256(run_provenance_payload),
        }
        if observed_hashes != expected_hashes[name]:
            raise AssertionError(
                f"installed {name} Enrichment complete RunProvenance payload changed: "
                f"{observed_hashes!r}"
            )
        parameters = result.provenance.workflow_parameters
        if parameters["correction_owner"] != "ora_engine":
            raise AssertionError("installed Enrichment correction owner changed")
        if parameters["number_of_tests"] != len(result.records):
            raise AssertionError("installed Enrichment tested family changed")
        reports[name] = {
            "records": _probe_json_value(record_payloads),
            "caveats": _probe_json_value(caveat_payloads),
            "identifier_set_provenance": {
                "selected": _probe_json_value(selected_payload),
                "background": _probe_json_value(background_payload),
            },
            "hashes": observed_hashes,
            "input_table_fingerprints": _probe_json_value(input_table_payloads),
            "output_table_fingerprints": _probe_json_value(output_table_payloads),
            "run_provenance_sha256": observed_hashes["run_provenance"],
            "correction_owner": parameters["correction_owner"],
            "number_of_tests": parameters["number_of_tests"],
        }
    return {
        "constructor_signature": str(inspect.signature(api.EnrichmentWorkflow)),
        "request_field_order": list(expected_request_parameters),
        "owner_identities": identities,
        "cases": reports,
    }


def _verify_references_evidence_import_orders(
    *,
    repo_root: pathlib.Path,
    environment_root: pathlib.Path,
) -> list[dict[str, object]]:
    owner_modules = {
        "references": (
            "phospy.references.builder",
            "phospy.references.models",
            "phospy.references.tables",
            "phospy.references.sequences.repository",
        ),
        "evidence": (
            "phospy.evidence.localisation",
            "phospy.evidence.sequence_resolution",
            "phospy.evidence.quantification_depth",
        ),
    }
    reports: list[dict[str, object]] = []
    for owners in (("references", "evidence"), ("evidence", "references")):
        modules = tuple(
            module_name for owner in owners for module_name in owner_modules[owner]
        )
        code = "\n".join(
            [
                "import importlib",
                "import inspect",
                "import pathlib",
                "import sys",
                "import pandas as pd",
                "repo_root = pathlib.Path(sys.argv[1]).resolve()",
                "environment_root = pathlib.Path(sys.argv[2]).resolve()",
                f"module_names = {modules!r}",
                "loaded = [importlib.import_module(name) for name in module_names]",
                "for module in loaded:",
                "    origin = pathlib.Path(module.__file__).resolve()",
                "    origin.relative_to(environment_root)",
                "    try:",
                "        origin.relative_to(repo_root)",
                "    except ValueError:",
                "        pass",
                "    else:",
                "        raise AssertionError(f'checkout import origin: {origin}')",
                "from phospy.references.builder import ReferenceBundleBuilder",
                "from phospy.references.tables import KinaseSubstrateReference",
                "from phospy.evidence.quantification_depth import validate_quantification_depth_series",
                "assert tuple(inspect.signature(ReferenceBundleBuilder).parameters) == ()",
                "assert isinstance(ReferenceBundleBuilder(), ReferenceBundleBuilder)",
                "table = KinaseSubstrateReference(frame=pd.DataFrame({"
                "'kinase': ['akt1'], 'substrate_site': ['MAPK1;S123;']}))",
                "assert table.frame.loc[0, 'kinase'] == 'AKT1'",
                "depth = validate_quantification_depth_series("
                "pd.Series([3.0], index=['site']), field_name='depth')",
                "assert depth.tolist() == [3.0]",
            ]
        )
        result = subprocess.run(
            [sys.executable, "-I", "-c", code, str(repo_root), str(environment_root)],
            cwd=pathlib.Path(tempfile.gettempdir()),
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            raise AssertionError(
                f"installed References/Evidence import order {owners!r} failed: "
                f"{result.stderr}"
            )
        reports.append({"order": list(owners), "status": "ok", "modules": list(modules)})
    return reports


def _verify_enrichment_import_orders(
    *,
    repo_root: pathlib.Path,
    environment_root: pathlib.Path,
) -> list[dict[str, object]]:
    orders = (
        (
            "phospy.io.readers.enrichment_sets",
            "phospy.enrichment.contracts.collections",
            "phospy.enrichment.contracts.requests",
        ),
        (
            "phospy.enrichment.contracts.requests",
            "phospy.enrichment.contracts.collections",
            "phospy.io.readers.enrichment_sets",
        ),
        (
            "phospy.api",
            "phospy.advanced",
            "phospy.enrichment.workflow",
            "phospy.enrichment.contracts.results",
        ),
        (
            "phospy.enrichment.contracts.results",
            "phospy.enrichment.workflow",
            "phospy.advanced",
            "phospy.api",
        ),
        (
            "phospy.provenance",
            "phospy.enrichment.provenance",
            "phospy.enrichment.contracts.results",
        ),
        (
            "phospy.enrichment.contracts.results",
            "phospy.enrichment.provenance",
            "phospy.provenance",
        ),
    )
    reports: list[dict[str, object]] = []
    for modules in orders:
        code = "\n".join(
            [
                "import importlib",
                "import pathlib",
                "import sys",
                "repo_root = pathlib.Path(sys.argv[1]).resolve()",
                "environment_root = pathlib.Path(sys.argv[2]).resolve()",
                f"module_names = {modules!r}",
                "loaded = [importlib.import_module(name) for name in module_names]",
                "for module in loaded:",
                "    origin = pathlib.Path(module.__file__).resolve()",
                "    origin.relative_to(environment_root)",
                "    try:",
                "        origin.relative_to(repo_root)",
                "    except ValueError:",
                "        pass",
                "    else:",
                "        raise AssertionError(f'checkout import origin: {origin}')",
                "import phospy.api as api",
                "from phospy.enrichment.contracts.requests import EnrichmentWorkflowRequest",
                "from phospy.enrichment.workflow import EnrichmentWorkflow",
                "assert api.EnrichmentWorkflowRequest is EnrichmentWorkflowRequest",
                "assert api.EnrichmentWorkflow is EnrichmentWorkflow",
                "workflow = api.EnrichmentWorkflow()",
                "gene = workflow.run(api.EnrichmentWorkflowRequest(",
                "    identifier_column='gene_symbol',",
                "    identifier_kind='gene_symbol',",
                "    set_collection=api.GeneSetCollection(",
                "        sets={'A': ('A',), 'B': ('B',)},",
                "        identifier_kind='gene_symbol',",
                "    ),",
                "    selected_identifiers=('A',),",
                "    background_universe=('A', 'B'),",
                "))",
                "ptm = workflow.run(api.EnrichmentWorkflowRequest(",
                "    identifier_column='site_key',",
                "    identifier_kind='site_key',",
                "    set_collection=api.PtmSetCollection(",
                "        sets={'M': ('rat|P1|S1',), 'N': ('rat|P1|T2',)},",
                "        identifier_kind='site_key',",
                "    ),",
                "    selected_identifiers=('rat|P1|S1',),",
                "    background_universe=('rat|P1|S1', 'rat|P1|T2'),",
                "))",
                "assert [record.term_id for record in gene.records] == ['A', 'B']",
                "assert [record.term_id for record in ptm.records] == ['M', 'N']",
                "assert gene.provenance.workflow_parameters['correction_owner'] == 'ora_engine'",
                "assert ptm.provenance.workflow_parameters['correction_owner'] == 'ora_engine'",
            ]
        )
        completed = subprocess.run(
            [sys.executable, "-I", "-c", code, str(repo_root), str(environment_root)],
            cwd=pathlib.Path(tempfile.gettempdir()),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"installed Enrichment import order {modules!r} failed: "
                f"{completed.stderr}"
            )
        reports.append(
            {
                "capability": "enrichment",
                "order": list(modules),
                "status": "ok",
                "gene_and_ptm_execution": "ok",
            }
        )
    return reports


def _verify_control_site_import_orders(
    *,
    repo_root: pathlib.Path,
    environment_root: pathlib.Path,
) -> list[dict[str, object]]:
    orders = (
        (
            "phospy.preprocessing.control_sites.discovery_contracts",
            "phospy.preprocessing.control_sites.discovery_provenance",
            "phospy.preprocessing.control_sites.discovery_results",
            "phospy.preprocessing.control_sites.workflow",
            "phospy.preprocessing.batch_correction.contracts.workflow",
            "phospy.preprocessing.batch_correction.planning",
        ),
        (
            "phospy.preprocessing.batch_correction.planning",
            "phospy.preprocessing.batch_correction.contracts.workflow",
            "phospy.preprocessing.control_sites.workflow",
            "phospy.preprocessing.control_sites.discovery_results",
            "phospy.preprocessing.control_sites.discovery_provenance",
            "phospy.preprocessing.control_sites.discovery_contracts",
        ),
    )
    reports: list[dict[str, object]] = []
    for modules in orders:
        code = "\n".join(
            [
                "import importlib",
                "import pathlib",
                "import sys",
                "import pandas as pd",
                "repo_root = pathlib.Path(sys.argv[1]).resolve()",
                "environment_root = pathlib.Path(sys.argv[2]).resolve()",
                f"module_names = {modules!r}",
                "loaded = [importlib.import_module(name) for name in module_names]",
                "for module in loaded:",
                "    origin = pathlib.Path(module.__file__).resolve()",
                "    origin.relative_to(environment_root)",
                "    try:",
                "        origin.relative_to(repo_root)",
                "    except ValueError:",
                "        pass",
                "    else:",
                "        raise AssertionError(f'checkout import origin: {origin}')",
                "from phospy.advanced import (ControlSiteSet, CorrectionMissingnessPolicy,",
                "    ObservationMask, OriginallyMissingCellTracking, SpsDiscoveryConfig,",
                "    SpsDiscoveryRequest, SpsDiscoveryResult, SpsDiscoveryWorkflow,",
                "    SpsReferenceDataset, SpsRuvBatchCorrectionConfig)",
                "from phospy.preprocessing.batch_correction import (",
                "    BatchCorrectionWorkflow, BatchCorrectionWorkflowRequest)",
                "sites = pd.Index((",
                "    'phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=P1|residue=S|position=1',",
                "    'phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=P2|residue=S|position=2',",
                "    'phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=P3|residue=S|position=3',",
                "), name='site_key')",
                "conditions = {'a1': 'a', 'a2': 'a', 'b1': 'b', 'b2': 'b'}",
                "references = tuple(SpsReferenceDataset.from_condition_relative_log2(",
                "    dataset_id=dataset_id, intensities=pd.DataFrame(values, index=sites),",
                "    condition_by_sample=conditions,",
                "    log2_scale_established_by=f'{dataset_id} scale',",
                "    baseline_centering_established_by=f'{dataset_id} baseline',",
                "    organism='rat', baseline_context='condition a',",
                "    reference_context='installed import-order fixture',",
                "    source_name=dataset_id, source_version='1',",
                "    source_uri=f'https://example.test/{dataset_id}')",
                "    for dataset_id, values in ((",
                "        'r1', {'a1': (0., 0., 0.), 'a2': (0., 0., 0.),",
                "               'b1': (.1, .2, 2.), 'b2': (.1, .2, 2.)}), (",
                "        'r2', {'a1': (0., 0., 0.), 'a2': (0., 0., 0.),",
                "               'b1': (.2, .3, 3.), 'b2': (.2, .3, 3.)})))",
                "result = SpsDiscoveryWorkflow().run(SpsDiscoveryRequest(",
                "    reference_datasets=references,",
                "    config=SpsDiscoveryConfig(top_n=2, minimum_shared_sites=2)))",
                "restored = SpsDiscoveryResult.from_payload(result.to_payload())",
                "controls = ControlSiteSet.from_payload(",
                "    restored.control_site_set.to_payload())",
                "assert controls.to_payload() == result.control_site_set.to_payload()",
                "samples = ('s1', 's2', 's3', 's4')",
                "target = pd.DataFrame({'s1': (1., 2., 5.), 's2': (2., 4., 5.5),",
                "    's3': (3., 1., 6.), 's4': (4., 3., 7.)}, index=sites)",
                "mask = ObservationMask(feature_ids=tuple(sites), sample_ids=samples)",
                "missingness = CorrectionMissingnessPolicy(",
                "    originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,",
                "    observation_mask=mask)",
                "public_config = SpsRuvBatchCorrectionConfig(",
                "    control_site_set=controls, batch_column='batch',",
                "    condition_columns=('condition',), missingness_policy=missingness,",
                "    n_unwanted_factors=1)",
                "corrected = BatchCorrectionWorkflow().run(BatchCorrectionWorkflowRequest(",
                "    phospho=target, config=public_config.to_internal_request(),",
                "    sample_metadata=pd.DataFrame({'batch': ('b1', 'b1', 'b2', 'b2'),",
                "        'condition': ('a', 'b', 'a', 'b')}, index=samples),",
                "    control_site_set=controls, missingness_policy=missingness))",
                "assert corrected.corrected_matrix.shape == target.shape",
                "assert corrected.provenance.selected_site_key_rows == result.selected_site_keys",
                "assert corrected.provenance.control_site_source['sps_discovery_identity'] == result.discovery_identity",
                "assert corrected.corrected_preprocessing_output.output_observation_mask is not None",
                "assert corrected.corrected_preprocessing_output.corrected_cell_status is not None",
            ]
        )
        completed = subprocess.run(
            [sys.executable, "-I", "-c", code, str(repo_root), str(environment_root)],
            cwd=pathlib.Path(tempfile.gettempdir()),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"installed Control-site/SPS import order {modules!r} failed: "
                f"{completed.stderr}"
            )
        reports.append(
            {
                "capability": "control_sites_sps",
                "order": list(modules),
                "status": "ok",
                "discovery_restore_and_correction": "ok",
            }
        )
    return reports


def _verify_protein_aware_boundary_import_orders(
    *,
    repo_root: pathlib.Path,
    environment_root: pathlib.Path,
) -> list[dict[str, object]]:
    reports: list[dict[str, object]] = []
    for dataset_first in (False, True):
        early_imports = ""
        stage_first_assertions = ""
        if dataset_first:
            early_imports = '''
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.datasets.internal_view import DatasetInternalView
from phospy.datasets.protein_aware import ProteinAwarePreparationReport, ProteinAwarePreparationResult
from phospy.run_lineage.models import RunProvenance
from phospy.datasets.construction.protein_aware import compose_protein_aware_preparation_result
'''
        else:
            stage_first_assertions = '''
assert not any(name == 'phospy.run_lineage' or name.startswith('phospy.run_lineage.') for name in sys.modules)
assert not any(name == 'phospy.datasets' or name.startswith('phospy.datasets.') for name in sys.modules)
assert not any(name == 'phospy.science.datasets' or name.startswith('phospy.science.datasets.') for name in sys.modules)
'''
        code = early_imports + r'''
import hashlib
import importlib
import importlib.util
import json
import pathlib
import sys
import pandas as pd

repo_root = pathlib.Path(sys.argv[1]).resolve()
environment_root = pathlib.Path(sys.argv[2]).resolve()
from phospy.preprocessing.contracts.configuration.total_protein import DatasetProteinAwarePreparationConfig
from phospy.preprocessing.protein_aware.mapping import ProteinMappingConfig
from phospy.preprocessing.results.protein_aware import ProteinAwarePreparationStageOutput
from phospy.preprocessing.stages.protein_aware import ProteinAwarePreparationStage
from phospy.quantitative.models import IntensityScaleState, MatrixIntensityScaleState

phospho = pd.DataFrame(
    {'sample_a': [1.0], 'sample_b': [2.0]},
    index=pd.Index(['MAPK14;Y182;'], name='site_id'),
)
site_metadata = pd.DataFrame(
    {
        'protein_id': ['P53778'],
        'gene_symbol': ['MAPK14'],
        'site': ['Y182'],
        'site_sequence': ['LDFGLARHTDDEMTGYVATRWYRAPEIMLNW'],
        'protein_identifier': ['MAPK14'],
        'protein_namespace': ['protein_id'],
        'localisation_confidence': [0.95],
    },
    index=phospho.index.copy(),
)
total = pd.DataFrame(
    {'sample_a': [10.0], 'sample_b': [11.0]},
    index=pd.Index(['P53778'], name='protein_id'),
)
state = IntensityScaleState(
    phospho=MatrixIntensityScaleState.log2(established_by='installed.phospho'),
    total=MatrixIntensityScaleState.log2(established_by='installed.total'),
)
stage_output = ProteinAwarePreparationStage().run(
    phospho=phospho,
    site_metadata=site_metadata,
    total=total,
    transformation_state=state,
    config=DatasetProteinAwarePreparationConfig(policy='prepare_model_inputs'),
    mapping_config=ProteinMappingConfig(protein_identifier_columns=('protein_id',)),
)
assert stage_output is not None
assert stage_output.__class__ is ProteinAwarePreparationStageOutput
''' + stage_first_assertions + r'''
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.datasets.internal_view import DatasetInternalView
from phospy.datasets.protein_aware import ProteinAwarePreparationReport, ProteinAwarePreparationResult
from phospy.run_lineage.models import RunProvenance
from phospy.datasets.construction.protein_aware import compose_protein_aware_preparation_result

result = compose_protein_aware_preparation_result(stage_output, provenance=None)
assert result.__class__ is ProteinAwarePreparationResult
assert result.report.__class__ is ProteinAwarePreparationReport
assert result.report.provenance is None
assert ProteinAwarePreparationStageOutput.__module__ == 'phospy.preprocessing.results.protein_aware'
assert AnalysisReadyPhosphoDataset.__module__ == 'phospy.datasets.models'
assert ProteinAwarePreparationReport.__module__ == 'phospy.datasets.protein_aware'
assert ProteinAwarePreparationResult.__module__ == 'phospy.datasets.protein_aware'
assert RunProvenance.__module__ == 'phospy.run_lineage.models'
assert DatasetInternalView.__module__ == 'phospy.datasets.internal_view'
payload = result.to_payload()
assert list(payload) == [
    'matched_pairs',
    'protein_covariate_matrix',
    'protein_covariate_matrix_shape',
    'report',
    'preparation_policy',
    'protein_mapping_policy',
    'provenance',
]
assert payload['report']['schema_version'] == 1
assert payload['provenance'] is None
import phospy
import phospy.api as api
import phospy.advanced as advanced
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    Organism,
)
dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale='linear',
        preprocessing_config=DatasetPreprocessingConfig(
            localisation=DatasetLocalisationConfig(
                mode='require_threshold',
                confidence_column='localisation_confidence',
                min_confidence=0.75,
            )
        ),
    )
)
assert dataset.__class__ is AnalysisReadyPhosphoDataset
dataset_view = DatasetInternalView(dataset)
assert dataset_view.phospho.equals(dataset.phospho)
for aggregate in (phospy, api, advanced):
    assert not hasattr(aggregate, 'ProteinAwarePreparationStageOutput')
for removed in (
    'phospy.contracts.dataset_build',
    'phospy.science.datasets.builders',
    'phospy.science.datasets.builders.protein_aware_binding',
    'phospy.science.datasets.preprocessing.protein_aware_models',
    'phospy.science.datasets.construction.protein_aware_binding',
    'phospy.science.datasets.internal_view',
    'phospy.science.datasets.protein_aware',
    'phospy.contracts.results.preprocessing',
):
    try:
        spec = importlib.util.find_spec(removed)
    except ModuleNotFoundError:
        spec = None
    assert spec is None
for module_name in (
    'phospy.preprocessing.stages.protein_aware',
    'phospy.preprocessing.results.protein_aware',
    'phospy.datasets.models',
    'phospy.datasets.construction.protein_aware',
    'phospy.datasets.internal_view',
    'phospy.datasets.protein_aware',
    'phospy.run_lineage.models',
):
    origin = pathlib.Path(importlib.import_module(module_name).__file__).resolve()
    origin.relative_to(environment_root)
    try:
        origin.relative_to(repo_root)
    except ValueError:
        pass
    else:
        raise AssertionError(f'checkout import origin: {origin}')
encoded = json.dumps(payload, sort_keys=True, separators=(',', ':')).encode('utf-8')
print(json.dumps({'payload_sha256': hashlib.sha256(encoded).hexdigest()}))
'''
        completed = subprocess.run(
            [sys.executable, "-I", "-c", code, str(repo_root), str(environment_root)],
            cwd=pathlib.Path(tempfile.gettempdir()),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            order_name = "dataset_first" if dataset_first else "stage_first"
            raise AssertionError(
                f"installed protein-aware import order {order_name} failed: "
                f"{completed.stderr}"
            )
        order_name = "dataset_first" if dataset_first else "stage_first"
        child_report = json.loads(completed.stdout.strip())
        reports.append(
            {
                "capability": "protein_aware_dataset_boundary",
                "order": order_name,
                "status": "ok",
                "payload_sha256": child_report["payload_sha256"],
            }
        )
    return reports


def _verify_batch_correction_boundary_import_orders(
    *,
    repo_root: pathlib.Path,
    environment_root: pathlib.Path,
) -> list[dict[str, object]]:
    orders = (
        (
            "phospy.provenance.serialization.batch_correction",
            "phospy.preprocessing.batch_correction.provenance.assembly",
        ),
        (
            "phospy.preprocessing.batch_correction.provenance.assembly",
            "phospy.provenance.serialization.batch_correction",
        ),
        (
            "phospy.preprocessing.batch_correction.adapters.preprocessing",
            "phospy.preprocessing.stages.batch_correction",
        ),
        (
            "phospy.preprocessing.stages.batch_correction",
            "phospy.preprocessing.batch_correction.adapters.preprocessing",
        ),
    )
    reports: list[dict[str, object]] = []
    for modules in orders:
        code = "\n".join(
            [
                "import importlib",
                "import pathlib",
                "import sys",
                "repo_root = pathlib.Path(sys.argv[1]).resolve()",
                "environment_root = pathlib.Path(sys.argv[2]).resolve()",
                f"module_names = {modules!r}",
                "loaded = [importlib.import_module(name) for name in module_names]",
                "for module in loaded:",
                "    origin = pathlib.Path(module.__file__).resolve()",
                "    origin.relative_to(environment_root)",
                "    try:",
                "        origin.relative_to(repo_root)",
                "    except ValueError:",
                "        pass",
                "    else:",
                "        raise AssertionError(f'checkout import origin: {origin}')",
                "if 'serialization' in module_names[0] or 'serialization' in module_names[1]:",
                "    from phospy.provenance.models import BatchCorrectionProvenance as aggregate",
                "    from phospy.provenance.models.batch_correction import BatchCorrectionProvenance as owner",
                "    import phospy.preprocessing.provenance.models as preprocessing_provenance",
                "    assert aggregate is owner",
                "    assert preprocessing_provenance._BatchCorrectionProvenance is owner",
                "else:",
                "    from phospy.preprocessing.batch_correction.adapters.preprocessing import SpsRuvStyleBatchCorrectionWorkflowRunner",
                "    from phospy.preprocessing.stages.batch_correction import BatchCorrectionStage",
                "    runner = SpsRuvStyleBatchCorrectionWorkflowRunner()",
                "    stage = BatchCorrectionStage(sps_ruv_runner=runner)",
                "    assert callable(runner.run)",
                "    assert stage._executor._sps_ruv_runner is runner",
            ]
        )
        completed = subprocess.run(
            [sys.executable, "-I", "-c", code, str(repo_root), str(environment_root)],
            cwd=pathlib.Path(tempfile.gettempdir()),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"installed batch-correction import order {modules!r} failed: "
                f"{completed.stderr}"
            )
        reports.append(
            {
                "capability": "batch_correction_boundary",
                "order": list(modules),
                "status": "ok",
            }
        )
    return reports


def _exercise_protein_aware_configuration_contract() -> dict[str, object]:
    preparation = DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")
    if preparation.policy != "prepare_model_inputs":
        raise AssertionError("protein-aware preparation policy changed")

    model = DifferentialProteinAwareModelConfig()
    if (
        model.method
        != DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
    ):
        raise AssertionError("protein-aware default method changed")

    differential_config = advanced.DifferentialAnalysisConfig(
        protein_aware_model=model,
    )
    if differential_config.protein_aware_model is not model:
        raise AssertionError("protein-aware differential config did not retain model")

    preprocessing_config = DatasetPreprocessingConfig(
        protein_aware_preparation=preparation,
    )
    if preprocessing_config.protein_aware_preparation is not preparation:
        raise AssertionError(
            "dataset preprocessing config did not retain protein-aware preparation"
        )

    return {
        "status": "constructed",
        "preparation_policy": str(preparation.policy),
        "method": str(model.method),
    }


def _verify_withdrawn_peptide_to_site_boundary() -> dict[str, object]:
    import phospy.differential as differential_public
    import phospy.differential.aggregation as aggregation_public
    from phospy.errors import PhosPyInputError
    from phospy.differential.aggregation.experimental import (
        PeptideToSiteAggregator,
    )
    from phospy.differential.aggregation.models import (
        PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
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
        ("phospy.differential", differential_public),
        ("phospy.differential.aggregation", aggregation_public),
    ):
        leaked = withdrawn_public_exports & set(getattr(module, "__all__", ()))
        if leaked:
            raise AssertionError(
                f"{module_name} still exports withdrawn post-hoc symbols: "
                f"{sorted(leaked)}"
            )

    if PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS != (
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
        "support_status": PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
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

    manifest = load_bundled_reference_manifest(Organism.RAT)
    if manifest.reference_id != "l6_native":
        raise AssertionError("typed bundled rat manifest reference_id changed")
    if manifest.reference_version != manifest_payload.get("reference_version"):
        raise AssertionError("typed manifest version does not match resource JSON")

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
    if dataset.preprocessing_report is None:
        raise AssertionError("default dataset build omitted its preprocessing report")
    no_op_correction = dataset.preprocessing_report.batch_correction
    if no_op_correction is None or no_op_correction.status != "disabled":
        raise AssertionError("default dataset build lost no-op correction execution")
    if no_op_correction.method != "none":
        raise AssertionError("default dataset build reported the wrong correction method")
    _assert_frame_snapshot_isolation(
        lambda: dataset.phospho,
        label="built dataset",
    )
    _assert_frame_snapshot_isolation(
        lambda: imputed_dataset.phospho,
        label="imputed dataset",
    )

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
    boundary_types = (
        DownstreamScoreSelectionPolicy,
        DownstreamScoreSource,
        KinaseScoreMatrix,
    )
    if any(
        boundary_type.__module__ != "phospy.kinase.result_boundary"
        for boundary_type in boundary_types
    ):
        raise AssertionError("installed Kinase boundary type ownership changed")
    if KinasePredictionMatrix.__module__ != "phospy.kinase.prediction.tables":
        raise AssertionError("installed prediction matrix ownership changed")
    if kinase_result.scoring_result.downstream_score_selection_policy is not (
        SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY
    ):
        raise AssertionError("installed Kinase selection policy identity changed")
    if KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES != {
        "kinase_non_default_reference_source",
        "kinase_reference_score_fallback",
        "kinase_non_phosr_equivalent_scoring",
    }:
        raise AssertionError("installed Kinase caveat grouping changed")
    if any(
        (
            hasattr(kinase_policies, "DownstreamScoreSource"),
            hasattr(
                kinase_caveats,
                "KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE",
            ),
            hasattr(
                kinase_caveats,
                "KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE",
            ),
            hasattr(kinase_caveats, "KINASE_SCORING_LIMITATION_CAVEAT_CODE"),
        )
    ):
        raise AssertionError("installed Kinase boundary retained a former route")
    expected_caveat_codes = {
        KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE:
            "kinase_non_default_reference_source",
        KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE:
            "kinase_reference_score_fallback",
        KINASE_SCORING_LIMITATION_CAVEAT_CODE:
            "kinase_non_phosr_equivalent_scoring",
    }
    if any(key != value for key, value in expected_caveat_codes.items()):
        raise AssertionError("installed Kinase boundary caveat identity changed")
    bundle_report = _exercise_result_bundle_round_trip(
        kinase_request=kinase_request,
        kinase_result=kinase_result,
    )
    signalome_report = _exercise_signalome_round_trip(kinase_result)
    publisher_report = _exercise_advanced_table_publisher(dataset)
    sps_ruv_report = _exercise_sps_ruv_contracts()

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
        "kinase_boundary": {
            "module": "phospy.kinase.result_boundary",
            "selection_policy": (
                kinase_result.scoring_result.downstream_score_selection_policy.name
            ),
            "caveat_code_count": len(KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES),
        },
        "bundle_round_trip": bundle_report,
        "signalome_round_trip": signalome_report,
        "advanced_table_publisher": publisher_report,
        "sps_ruv": sps_ruv_report,
    }


def _exercise_sps_ruv_contracts() -> dict[str, object]:
    owner_identities = {
        "ControlSiteSet": ControlSiteSet is OwnerControlSiteSet,
        "SpsDiscoveryRequest": SpsDiscoveryRequest is OwnerSpsDiscoveryRequest,
        "SpsDiscoveryResult": SpsDiscoveryResult is OwnerSpsDiscoveryResult,
        "SpsDiscoveryWorkflow": SpsDiscoveryWorkflow is OwnerSpsDiscoveryWorkflow,
        "SpsReferenceDataset": SpsReferenceDataset is OwnerSpsReferenceDataset,
    }
    if not all(owner_identities.values()):
        raise AssertionError(
            f"installed SPS/control public owner identities changed: {owner_identities!r}"
        )
    signature_parameters = {
        "ControlSiteSet": tuple(inspect.signature(ControlSiteSet).parameters),
        "ControlSiteSet.from_payload": tuple(
            inspect.signature(ControlSiteSet.from_payload).parameters
        ),
        "SpsDiscoveryRequest": tuple(
            inspect.signature(SpsDiscoveryRequest).parameters
        ),
        "SpsDiscoveryResult": tuple(
            inspect.signature(SpsDiscoveryResult).parameters
        ),
        "SpsDiscoveryResult.from_payload": tuple(
            inspect.signature(SpsDiscoveryResult.from_payload).parameters
        ),
        "SpsDiscoveryWorkflow": tuple(
            inspect.signature(SpsDiscoveryWorkflow).parameters
        ),
        "SpsDiscoveryWorkflow.run": tuple(
            inspect.signature(SpsDiscoveryWorkflow.run).parameters
        ),
        "SpsReferenceDataset": tuple(
            inspect.signature(SpsReferenceDataset).parameters
        ),
    }
    expected_signature_parameters = {
        "ControlSiteSet": (
            "annotations",
            "source_metadata",
            "unannotated_site_status",
        ),
        "ControlSiteSet.from_payload": ("payload",),
        "SpsDiscoveryRequest": ("reference_datasets", "config"),
        "SpsDiscoveryResult": (
            "selected_site_keys",
            "site_ranking",
            "provenance",
        ),
        "SpsDiscoveryResult.from_payload": ("payload",),
        "SpsDiscoveryWorkflow": (),
        "SpsDiscoveryWorkflow.run": ("self", "request"),
        "SpsReferenceDataset": (
            "dataset_id",
            "intensities",
            "condition_by_sample",
            "intensity_scale_state",
            "source_name",
            "source_version",
            "source_uri",
            "organism",
            "baseline_context",
            "reference_context",
        ),
    }
    if signature_parameters != expected_signature_parameters:
        raise AssertionError(
            "installed SPS/control signatures changed: "
            f"{signature_parameters!r}"
        )

    site_keys = pd.Index(
        (
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0001|residue=S|position=1",
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0002|residue=S|position=2",
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            "protein_identifier=P0003|residue=S|position=3",
        ),
        name="site_key",
    )
    condition_by_sample = {
        "a_1": "a",
        "a_2": "a",
        "b_1": "b",
        "b_2": "b",
    }
    references = tuple(
        SpsReferenceDataset.from_condition_relative_log2(
            dataset_id=dataset_id,
            intensities=pd.DataFrame(values, index=site_keys),
            condition_by_sample=condition_by_sample,
            log2_scale_established_by=f"{dataset_id} governed log2 transform",
            baseline_centering_established_by=(
                f"{dataset_id} governed condition-relative centering"
            ),
            organism="rat",
            baseline_context="condition a",
            reference_context="installed-distribution SPS smoke fixture",
            source_name=f"governed-{dataset_id}",
            source_version="1",
            source_uri=f"https://example.test/sps/{dataset_id}",
        )
        for dataset_id, values in (
            (
                "reference_a",
                {
                    "a_1": (0.0, 0.0, 0.0),
                    "a_2": (0.0, 0.0, 0.0),
                    "b_1": (0.1, 0.2, 2.0),
                    "b_2": (0.1, 0.2, 2.0),
                },
            ),
            (
                "reference_b",
                {
                    "a_1": (0.0, 0.0, 0.0),
                    "a_2": (0.0, 0.0, 0.0),
                    "b_1": (0.2, 0.3, 3.0),
                    "b_2": (0.2, 0.3, 3.0),
                },
            ),
        )
    )
    discovery = SpsDiscoveryWorkflow().run(
        SpsDiscoveryRequest(
            reference_datasets=references,
            config=SpsDiscoveryConfig(top_n=2, minimum_shared_sites=2),
        )
    )
    discovery_payload = discovery.to_payload()
    current_restored = SpsDiscoveryResult.from_payload(discovery_payload)
    controls = ControlSiteSet.from_payload(
        current_restored.control_site_set.to_payload()
    )
    if len(controls.annotations) != 2:
        raise AssertionError("installed SPS discovery did not select two controls")
    if current_restored.to_payload() != discovery_payload:
        raise AssertionError("installed current SPS result restoration changed payload")

    legacy_payload = json.loads(json.dumps(discovery_payload))
    legacy_payload["identity_schema"] = "phospy-sps-discovery-result-v1"
    legacy_payload["provenance"]["algorithm_version"] = "1.0.0"
    legacy_payload["provenance"].pop("cross_count_ordering_policy")
    for record in legacy_payload["site_ranking"]:
        for statistic in record["dataset_statistics"]:
            statistic.pop("rank_quantile")
    legacy_identity_payload = {
        key: value
        for key, value in legacy_payload.items()
        if key != "discovery_identity"
    }
    legacy_digest = hash_json_payload(legacy_identity_payload)
    legacy_payload["discovery_identity"] = f"sha256-stable-json-v1:{legacy_digest}"
    legacy_restored = SpsDiscoveryResult.from_payload(legacy_payload)
    if legacy_restored.provenance.algorithm_version != "1.0.0":
        raise AssertionError("installed legacy SPS algorithm identity was rewritten")
    if legacy_restored.discovery_identity == discovery.discovery_identity:
        raise AssertionError("installed legacy SPS interpretation lost its identity")
    legacy_controls_payload = legacy_restored.control_site_set.to_payload()
    legacy_controls = ControlSiteSet.from_payload(legacy_controls_payload)
    if legacy_controls.to_payload() != legacy_controls_payload:
        raise AssertionError("installed legacy SPS controls changed during restoration")
    if legacy_controls.source_metadata.source_version != "1.0.0":
        raise AssertionError("installed legacy SPS control source version was rewritten")
    if (
        legacy_controls.source_metadata.sps_discovery_identity
        != legacy_restored.discovery_identity
    ):
        raise AssertionError("installed legacy SPS control lineage identity changed")

    samples = ("sample_1", "sample_2", "sample_3", "sample_4")
    phospho = pd.DataFrame(
        {
            "sample_1": (1.0, 2.0, 5.0),
            "sample_2": (2.0, 4.0, 5.5),
            "sample_3": (3.0, 1.0, 6.0),
            "sample_4": (4.0, 3.0, 7.0),
        },
        index=site_keys,
    )
    phospho_before_linear = phospho.copy(deep=True)
    linear = LinearResidualizeBatchCorrectionEngine().run(
        phospho=phospho,
        batch_labels=("batch_1", "batch_1", "batch_2", "batch_2"),
        condition_labels=("a", "b", "a", "b"),
        config=DatasetBatchCorrectionConfig(method="linear_residualize_batch"),
    )
    if linear.report.status != "applied":
        raise AssertionError("installed linear residualisation was not applied")
    if linear.corrected_matrix.shape != phospho.shape:
        raise AssertionError("installed linear residualisation changed matrix shape")
    if not phospho.equals(phospho_before_linear):
        raise AssertionError("installed linear residualisation mutated its input")
    correction_environment = collect_batch_correction_environment_provenance(
        use_cache=False
    )
    if tuple(correction_environment.dependency_versions) != (
        BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
    ):
        raise AssertionError("installed correction dependency order changed")
    expected_runtime_dependencies = {
        dependency: metadata.version(dependency)
        for dependency in BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
    }
    if correction_environment.python_version != sys.version.split()[0]:
        raise AssertionError(
            "installed correction environment recorded the wrong Python version"
        )
    if dict(correction_environment.dependency_versions) != (
        expected_runtime_dependencies
    ):
        raise AssertionError(
            "installed correction environment recorded the wrong dependency versions"
        )
    if correction_environment.package_version != str(phospy.__version__):
        raise AssertionError(
            "installed correction environment recorded the wrong PhosPy version"
        )
    missing_input = phospho.copy(deep=True)
    governed_cell = (str(missing_input.index[2]), "sample_1")
    missing_input.loc[governed_cell[0], governed_cell[1]] = float("nan")
    governed_mask = ObservationMask(
        feature_ids=tuple(str(value) for value in missing_input.index),
        sample_ids=tuple(str(value) for value in missing_input.columns),
        originally_missing_cells=(governed_cell,),
    )
    temporary_policy = TemporaryImputationPolicy(
        allowed=True,
        method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
        method_parameters=(("min_observed_values", 2),),
    )
    prepared_missingness = prepare_ruv_iii_missingness(
        phospho=missing_input,
        plan=SimpleNamespace(
            observation_mask=governed_mask,
            temporary_imputation_policy=temporary_policy,
        ),
    )
    if not prepared_missingness.temporary_completion_applied:
        raise AssertionError("installed RUV-III omitted temporary completion")
    if pd.isna(prepared_missingness.working.loc[governed_cell[0], governed_cell[1]]):
        raise AssertionError("installed RUV-III did not complete the governed cell")
    restored_missingness = prepared_missingness.working.mask(
        prepared_missingness.actual_missing,
        float("nan"),
    )
    if not pd.isna(restored_missingness.loc[governed_cell[0], governed_cell[1]]):
        raise AssertionError("installed RUV-III did not restore the governed cell")
    if not pd.isna(missing_input.loc[governed_cell[0], governed_cell[1]]):
        raise AssertionError("installed RUV-III preparation mutated its input")
    observation_mask = ObservationMask(
        feature_ids=tuple(site_keys),
        sample_ids=samples,
    )
    missingness_policy = CorrectionMissingnessPolicy(
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        observation_mask=observation_mask,
    )
    public_correction_config = SpsRuvBatchCorrectionConfig(
        control_site_set=controls,
        batch_column="batch",
        condition_columns=("condition",),
        missingness_policy=missingness_policy,
        n_unwanted_factors=1,
    )
    corrected = BatchCorrectionWorkflow().run(
        BatchCorrectionWorkflowRequest(
            phospho=phospho,
            config=public_correction_config.to_internal_request(),
            sample_metadata=pd.DataFrame(
                {
                    "batch": ("batch_1", "batch_1", "batch_2", "batch_2"),
                    "condition": ("a", "b", "a", "b"),
                },
                index=samples,
            ),
            control_site_set=controls,
            missingness_policy=missingness_policy,
        )
    )
    if corrected.corrected_matrix.shape != phospho.shape:
        raise AssertionError("installed correction workflow changed matrix shape")
    correction_output = corrected.corrected_preprocessing_output
    if correction_output.output_observation_mask is None:
        raise AssertionError("installed correction omitted its observation mask")
    if correction_output.corrected_cell_status is None:
        raise AssertionError("installed correction omitted its cell classification")
    correction_provenance = batch_correction_provenance_to_payload(
        corrected.provenance
    )
    restored_correction = batch_correction_provenance_from_payload(
        correction_provenance
    )
    if batch_correction_provenance_to_payload(restored_correction) != correction_provenance:
        raise AssertionError("installed direct correction provenance round trip changed")
    if restored_correction.__class__ is not corrected.provenance.__class__:
        raise AssertionError("installed direct correction provenance identity changed")
    if corrected.provenance.python_version != sys.version.split()[0]:
        raise AssertionError("installed correction recorded the wrong Python version")
    if dict(corrected.provenance.dependency_versions) != expected_runtime_dependencies:
        raise AssertionError("installed correction recorded the wrong dependency versions")
    if corrected.provenance.phospy_version != str(phospy.__version__):
        raise AssertionError("installed correction recorded the wrong PhosPy version")
    malformed_correction = dict(correction_provenance)
    malformed_correction["schema_version"] = 999
    try:
        batch_correction_provenance_from_payload(malformed_correction)
    except api.PhosPyInputError as exc:
        if type(exc).__name__ != "PhosPyInputError":
            raise AssertionError(
                "installed malformed correction payload raised the wrong error type"
            ) from exc
    else:
        raise AssertionError("installed malformed correction payload was accepted")
    if corrected.provenance.selected_site_key_rows != discovery.selected_site_keys:
        raise AssertionError("installed correction lost restored selected controls")
    correction_sps_identity = corrected.provenance.control_site_source.get(
        "sps_discovery_identity"
    )
    if correction_sps_identity != discovery.discovery_identity:
        raise AssertionError("installed correction lost SPS discovery lineage")
    input_fingerprint = corrected.provenance.input_matrix_fingerprint
    output_fingerprint = corrected.provenance.output_matrix_fingerprint
    if output_fingerprint is None:
        raise AssertionError("installed correction omitted its output fingerprint")
    correction_stage = PreprocessingStageProvenance(
        stage="batch_correction",
        operation=corrected.provenance.requested_method,
        parameters={"method": corrected.provenance.requested_method},
        input_shape=(input_fingerprint.rows, input_fingerprint.columns),
        output_shape=(output_fingerprint.rows, output_fingerprint.columns),
        input_hash=input_fingerprint.exact_hash_value,
        output_hash=output_fingerprint.exact_hash_value,
        phospho_input_hash=input_fingerprint.exact_hash_value,
        phospho_output_hash=output_fingerprint.exact_hash_value,
        dropped_row_ids=(),
        dropped_row_count=0,
        batch_correction_provenance=corrected.provenance,
    )
    correction_stage_payload = stage_to_payload(correction_stage)
    restored_correction_stage = stage_from_payload(correction_stage_payload)
    if stage_to_payload(restored_correction_stage) != correction_stage_payload:
        raise AssertionError("installed direct preprocessing-stage round trip changed")
    if restored_correction_stage.__class__ is not correction_stage.__class__:
        raise AssertionError("installed preprocessing-stage provenance identity changed")
    correction_run = RunProvenance(
        environment=correction_environment,
        input_tables=(input_fingerprint,),
        preprocessing_stages=(correction_stage,),
        reference=None,
        workflow_name="installed_batch_correction_probe",
        workflow_parameters={"method": corrected.provenance.requested_method},
        random_state=None,
        random_seed_policy=None,
        output_tables=(output_fingerprint,),
    )
    nested_correction_payload = run_provenance_to_payload(correction_run)
    restored_run = run_provenance_from_payload(nested_correction_payload)
    if run_provenance_to_payload(restored_run) != nested_correction_payload:
        raise AssertionError("installed nested correction provenance round trip changed")
    restored_nested_correction = (
        restored_run.preprocessing_stages[0].batch_correction_provenance
    )
    if restored_nested_correction is None:
        raise AssertionError("installed nested correction provenance was omitted")
    if restored_nested_correction.__class__ is not corrected.provenance.__class__:
        raise AssertionError("installed nested correction provenance identity changed")
    correction_payload = {
        "corrected_matrix": _probe_frame_payload(corrected.corrected_matrix),
        "observation_mask": _probe_frame_payload(
            correction_output.output_observation_mask
        ),
        "cell_classification": _probe_frame_payload(
            correction_output.corrected_cell_status
        ),
        "diagnostics": corrected.diagnostics,
        "warnings": list(corrected.warnings),
        "caveats": [item.to_payload() for item in corrected.caveats],
        "provenance": correction_provenance,
        "fingerprints": {
            "input_matrix": correction_provenance["input_matrix_fingerprint"],
            "output_matrix": correction_provenance["output_matrix_fingerprint"],
            "observation_masks": correction_provenance["observation_masks"],
        },
        "sps_discovery_identity": correction_sps_identity,
    }
    expected_corrected_values = (
        (1.9999999999999996, 2.9999999999999996, 1.9999999999999996, 3.0),
        (1.5000000000000004, 3.5000000000000004, 1.5, 3.5),
        (5.625, 6.125, 5.375, 6.375),
    )
    import numpy as np

    np.testing.assert_allclose(
        corrected.corrected_matrix.to_numpy(dtype=float),
        np.asarray(expected_corrected_values, dtype=float),
        rtol=1e-10,
        atol=1e-12,
        equal_nan=True,
        err_msg="installed correction numerical matrix changed",
    )
    invariant_matrix = dict(correction_payload["corrected_matrix"])
    invariant_matrix["values"] = "compared-with-rtol-1e-10-atol-1e-12"
    invariant_correction_payload = {
        "corrected_matrix": invariant_matrix,
        "observation_mask": correction_payload["observation_mask"],
        "cell_classification": correction_payload["cell_classification"],
        "warnings": correction_payload["warnings"],
        "caveats": correction_payload["caveats"],
        "provenance_contract": {
            "schema_version": correction_provenance["schema_version"],
            "requested_method": correction_provenance["requested_method"],
            "preprocessing_stage_order": correction_provenance[
                "preprocessing_stage_order"
            ],
            "selected_site_key_rows": correction_provenance[
                "selected_site_key_rows"
            ],
            "python_version": "<verified-runtime-python>",
            "dependency_versions": {
                dependency: "<verified-runtime-version>"
                for dependency in BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
            },
            "output_fingerprint_present": (
                correction_provenance["output_matrix_fingerprint"] is not None
            ),
            "observation_mask_fingerprint_count": len(
                correction_provenance["observation_masks"]
            ),
            "sps_discovery_identity_preserved": (
                correction_sps_identity == discovery.discovery_identity
            ),
        },
    }
    sps_payload_hashes = {
        "current_result": _probe_payload_sha256(discovery_payload),
        "current_restored_result": _probe_payload_sha256(
            current_restored.to_payload()
        ),
        "legacy_input": _probe_payload_sha256(legacy_payload),
        "legacy_restored_result": _probe_payload_sha256(
            legacy_restored.to_payload()
        ),
        "restored_controls": _probe_payload_sha256(controls.to_payload()),
        "legacy_restored_controls": _probe_payload_sha256(
            legacy_controls.to_payload()
        ),
        "complete_correction": _probe_payload_sha256(invariant_correction_payload),
    }
    expected_sps_payload_hashes = {
        "current_result": "aa90a52d702fe74d797a6fb174234344c081a8f659e6863ce7748b324361f5bc",
        "current_restored_result": "aa90a52d702fe74d797a6fb174234344c081a8f659e6863ce7748b324361f5bc",
        "legacy_input": "47d33a0af8c9775efd189b887266e180638ee3f9f73fb27a7bb43832d6cf24c2",
        "legacy_restored_result": "847627653fbe197013eb3a180d8f7d4b683e7d098e8fecf52e35af75d0b55ff9",
        "restored_controls": "f8c7ca82e341f48e8cdfdd6c4ce10113bd4c2761ef6b3803c51b008d316a8d11",
        "legacy_restored_controls": "6f43b82b9ebaec9779fe1f9cf4fb02751369c715e797212d630a278a5e15fe23",
        "complete_correction": "aff9ff52945a7f00cb4470455d16eb6ec55792cca39721e02aacb3431c54be55",
    }
    if sps_payload_hashes != expected_sps_payload_hashes:
        raise AssertionError(
            f"complete SPS discovery/control/correction payload changed: "
            f"{sps_payload_hashes!r}"
        )

    replicate_structure = RuvIIIReplicateStructure.from_assignments(
        sample_order=samples,
        replicate_by_sample={
            "sample_1": "replicate_1",
            "sample_2": "replicate_1",
            "sample_3": "replicate_2",
            "sample_4": "replicate_2",
        },
    )
    ruv_iii_corrected = run_ruv_iii(
        phospho,
        control_site_keys=discovery.selected_site_keys,
        replicate_structure=replicate_structure,
        k=1,
    )
    if ruv_iii_corrected.corrected_matrix.shape != phospho.shape:
        raise AssertionError("installed RUV-III changed the phospho matrix shape")
    if ruv_iii_corrected.diagnostics.method != "ruv_iii":
        raise AssertionError("installed RUV-III reported an unexpected method")
    return {
        "sps_selected_control_count": len(controls.annotations),
        "sps_discovery_identity": discovery.discovery_identity,
        "sps_legacy_discovery_identity": legacy_payload["discovery_identity"],
        "sps_payload_hashes": sps_payload_hashes,
        "sps_owner_identities": owner_identities,
        "sps_signature_parameters": {
            key: list(value) for key, value in signature_parameters.items()
        },
        "correction_method": corrected.provenance.requested_method,
        "correction_output_shape": list(corrected.corrected_matrix.shape),
        "correction_sps_discovery_identity": correction_sps_identity,
        "correction_fingerprints": correction_payload["fingerprints"],
        "direct_correction_round_trip": "identity_preserved",
        "nested_correction_round_trip": "identity_preserved",
        "malformed_correction_payload": "PhosPyInputError",
        "runtime_correction_environment": {
            "python": corrected.provenance.python_version,
            "dependencies": dict(corrected.provenance.dependency_versions),
        },
        "linear_method": linear.report.method,
        "linear_output_shape": list(linear.corrected_matrix.shape),
        "correction_environment_dependencies": list(
            correction_environment.dependency_versions
        ),
        "ruv_iii_temporary_completion": "restored_governed_cell",
        "ruv_iii_method": ruv_iii_corrected.diagnostics.method,
        "ruv_iii_output_shape": list(ruv_iii_corrected.corrected_matrix.shape),
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
        processing_state_payload = manifest_payload["processing_state"]
        if "batch_correction" in processing_state_payload:
            raise AssertionError(
                "published processing-state payload unexpectedly contains batch_correction"
            )
        return {
            "manifest_output_format": manifest_payload["output_format"],
            "processing_state_batch_correction_omitted": True,
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
        if loaded.result.dataset.processing_state != kinase_result.dataset.processing_state:
            raise AssertionError("installed artifact bundle dataset processing state changed")
        pd.testing.assert_frame_equal(
            loaded.result.dataset.phospho,
            kinase_result.dataset.phospho,
            check_dtype=False,
        )
        _assert_frame_snapshot_isolation(
            lambda: loaded.result.dataset.phospho,
            label="reconstructed kinase dataset",
        )
        _assert_frame_snapshot_isolation(
            lambda: loaded.result.prediction_result.pred_mat,
            label="reconstructed kinase prediction",
        )
        manifest_payload = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
        phospho_entry = manifest_payload["dataset"]["tables"]["phospho"]
        if not isinstance(phospho_entry, dict):
            raise AssertionError("installed artifact bundle table entry is not content-addressed")
        for required_key in ("path", "sha256", "byte_size", "logical_type", "shape"):
            if required_key not in phospho_entry:
                raise AssertionError(f"bundle table entry is missing {required_key!r}")
        manifest_before = (bundle_root / "manifest.json").read_bytes()
        publication_root = pathlib.Path(tmp_dir) / "kinase_publication"
        published = publish_kinase_workflow(
            kinase_result,
            publication_root,
            output_format="csv",
        )
        _assert_published_table_isolation(published, label="kinase publication")
        if (bundle_root / "manifest.json").read_bytes() != manifest_before:
            raise AssertionError("kinase mutation probes changed the bundle manifest")
        return {
            "manifest_version": int(loaded.manifest_version),
            "phospho_entry_logical_type": phospho_entry["logical_type"],
            "phospho_entry_shape": phospho_entry["shape"],
            "processing_state_round_trip": True,
            "mutation_isolation": True,
            "publication_written_key_count": len(published),
        }


def _exercise_signalome_round_trip(kinase_result) -> dict[str, object]:
    config = SignalomeConfig(
        scientific=SignalomeScientificConfig(substrate_support_cutoff=0.0),
        clustering=SignalomeClusteringConfig(module_count=1),
        validation=SignalomeValidationConfig(
            localisation_requirement=LocalisationRequirement(),
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        output=SignalomeOutputConfig(
            network_min_paired_finite_observations=3,
        ),
        mode="exploratory_compatibility",
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=config,
    )
    result = SignalomeWorkflow().run(request)
    if result.module_assignments.table.empty:
        raise AssertionError("installed signalome module assignments are empty")
    if api.SignalomeWorkflow is not SignalomeWorkflow:
        raise AssertionError("installed SignalomeWorkflow aggregate identity changed")
    if api.SignalomeWorkflowRequest is not SignalomeWorkflowRequest:
        raise AssertionError("installed SignalomeWorkflowRequest aggregate identity changed")

    with tempfile.TemporaryDirectory(prefix="phospy-installed-signalome-bundle-") as tmp_dir:
        root = pathlib.Path(tmp_dir)
        bundle_root = root / "signalome_bundle"
        snapshot = SignalomeWorkflowConfigSnapshot.from_request(request)
        written = save_signalome_workflow_bundle(
            result,
            bundle_root,
            config_snapshot=snapshot,
            output_format="csv",
        )
        if written["manifest"] != bundle_root / "manifest.json":
            raise AssertionError("signalome bundle writer returned an unexpected manifest path")
        loaded = load_signalome_workflow_bundle(bundle_root)
        if loaded.config_snapshot != snapshot:
            raise AssertionError("installed signalome bundle config snapshot changed")
        if loaded.result.provenance != result.provenance:
            raise AssertionError("installed signalome bundle provenance changed")
        if loaded.result.dataset.processing_state != result.dataset.processing_state:
            raise AssertionError("installed signalome dataset processing state changed")
        pd.testing.assert_frame_equal(
            loaded.result.module_assignments.table,
            result.module_assignments.table,
            check_dtype=False,
        )
        _assert_frame_snapshot_isolation(
            lambda: loaded.result.dataset.phospho,
            label="reconstructed signalome dataset",
        )
        _assert_frame_snapshot_isolation(
            lambda: loaded.result.module_assignments.table,
            label="reconstructed signalome assignments",
        )
        manifest_before = (bundle_root / "manifest.json").read_bytes()
        published = publish_signalome_workflow(
            result,
            root / "signalome_publication",
            output_format="csv",
        )
        _assert_published_table_isolation(published, label="signalome publication")
        if (bundle_root / "manifest.json").read_bytes() != manifest_before:
            raise AssertionError("signalome mutation probes changed the bundle manifest")
        return {
            "manifest_version": int(loaded.manifest_version),
            "module_assignment_shape": [
                int(result.module_assignments.table.shape[0]),
                int(result.module_assignments.table.shape[1]),
            ],
            "processing_state_round_trip": True,
            "mutation_isolation": True,
            "publication_written_key_count": len(published),
        }


def _assert_frame_snapshot_isolation(getter, *, label: str) -> None:
    baseline = getter()
    mutated = getter()
    if mutated.empty:
        raise AssertionError(f"{label} mutation probe requires a non-empty frame")
    column_dtype = mutated.dtypes.iloc[0]
    current = mutated.iat[0, 0]
    if pd.api.types.is_bool_dtype(column_dtype):
        replacement = not bool(current)
    elif pd.api.types.is_numeric_dtype(column_dtype):
        replacement = 1.0 if pd.isna(current) else float(current) + 1.0
    else:
        replacement = "__consumer_mutation__"
    mutated.iat[0, 0] = replacement
    if not getter().equals(baseline):
        raise AssertionError(f"{label} public frame aliases owner state")


def _assert_published_table_isolation(written, *, label: str) -> None:
    candidates = sorted(
        path for path in written.values() if pathlib.Path(path).suffix == ".csv"
    )
    if not candidates:
        raise AssertionError(f"{label} did not write a CSV table")
    path = pathlib.Path(candidates[0])
    baseline = path.read_bytes()
    frame = pd.read_csv(path)
    if frame.empty:
        raise AssertionError(f"{label} mutation probe requires a non-empty table")
    frame.iloc[0, 0] = "__consumer_mutation__"
    if path.read_bytes() != baseline:
        raise AssertionError(f"{label} in-memory mutation changed persisted bytes")


def _build_dataset(*, include_missing: bool):
    index = ["TSC2;S939;", "GSK3B;S9;", "MAPK14;Y182;", "GSK3A;S21;"]
    phospho = pd.DataFrame(
        {
            "A_1": [10.0, 9.0, 7.5, 8.0],
            "A_2": [10.4, 9.1, 7.7, 8.2],
            "A_3": [9.8, 8.8, 7.8, 7.9],
            "B_1": [12.0, 9.3, 7.4, 9.1],
            "B_2": [12.2, 9.4, 7.6, 9.0],
            "B_3": [11.7, 9.2, 7.3, 9.4],
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
                sample_id="A_3",
                condition="A",
                biological_replicate_id="A_3",
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
            SampleDesignRecord(
                sample_id="B_3",
                condition="B",
                biological_replicate_id="B_3",
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
                sample_id="A_3",
                condition="A",
                biological_replicate_id="A_3",
                block_id="pair_3",
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
            SampleDesignRecord(
                sample_id="B_3",
                condition="B",
                biological_replicate_id="B_3",
                block_id="pair_3",
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
""".replace("__RETIRED_CAPABILITY_ROUTES__", repr(RETIRED_CAPABILITY_ROUTES))
    .replace("__RETIRED_HORIZONTAL_ROUTES__", repr(RETIRED_HORIZONTAL_ROUTES))
    .replace(
        "__RETIRED_ROUTE_IMPORT_PROBE_SOURCE__",
        repr(RETIRED_ROUTE_IMPORT_PROBE_SOURCE),
    )
)


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
    distribution_version: str
    runtime_version: str
    environment_root: Path
    run_directory: Path
    phospy_file: Path
    python_version: str
    requires_python: str
    resource_count: int
    ticket_1_boundary_status: str
    constraint_path: Path | None
    constraint_sha256: str | None
    installation_method: str = "isolated venv; pip install artifact path"
    probe_results: dict[str, Any] = field(default_factory=dict)


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
    distribution_versions = {report.distribution_version for report in reports}
    if len(distribution_versions) != 1:
        raise InstalledDistributionVerificationError(
            "wheel and sdist metadata Version values differ: "
            + ", ".join(sorted(distribution_versions))
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
    distribution_version = _artifact_version(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
    )
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
    probe_path = run_directory / "installed_distribution_probe.py"
    probe_path.write_text(INSTALLED_PROBE_SOURCE, encoding="utf-8", newline="\n")
    result = _run(
        [
            str(environment_python),
            "-I",
            str(probe_path),
            str(repo_root),
            str(environment_root.resolve()),
            artifact_kind,
        ],
        cwd=run_directory,
        repo_root=repo_root,
        context=f"execute installed {artifact_kind} probe",
    )
    payload = _probe_payload(result.stdout, artifact_kind=artifact_kind)
    runtime_version = _required_string(payload, "runtime_version")
    if runtime_version != distribution_version:
        raise InstalledDistributionVerificationError(
            f"installed {artifact_kind} runtime version mismatch for {artifact_path}: "
            f"metadata Version={distribution_version!r}, "
            f"phospy.__version__={runtime_version!r}"
        )
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
        distribution_version=distribution_version,
        runtime_version=runtime_version,
        environment_root=environment_root.resolve(),
        run_directory=run_directory.resolve(),
        phospy_file=phospy_file,
        python_version=_required_string(payload, "python"),
        requires_python=requires_python,
        resource_count=int(resource_count),
        ticket_1_boundary_status=ticket_1_boundary_status,
        constraint_path=None if constraint is None else constraint.resolve(),
        constraint_sha256=None if constraint is None else _sha256_file(constraint),
        installation_method=(
            "isolated venv; pip --no-cache-dir install "
            + ("with dependencies" if install_dependencies else "--no-deps")
            + (" and build isolation" if build_isolation else " --no-build-isolation")
        ),
        probe_results={
            key: value
            for key, value in payload.items()
            if key
            in {
                "public_surface_report",
                "reference_builder_report",
                "import_order_report",
                "resource_report",
                "retired_route_report",
                "prediction_identity_report",
                "scientific_report",
                "enrichment_report",
                "environment_versions",
            }
        },
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
    value = _artifact_metadata_value(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
        field_name="Requires-Python",
    )
    if value is None or not value.strip():
        raise InstalledDistributionVerificationError(
            f"{artifact_kind} metadata is missing Requires-Python: {artifact_path}"
        )
    return value.strip()


def _artifact_version(*, artifact_kind: str, artifact_path: Path) -> str:
    value = _artifact_metadata_value(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
        field_name="Version",
    )
    if value is None or not value.strip():
        raise InstalledDistributionVerificationError(
            f"{artifact_kind} metadata is missing Version: {artifact_path}"
        )
    return value.strip()


def _artifact_metadata_value(
    *,
    artifact_kind: str,
    artifact_path: Path,
    field_name: str,
) -> str | None:
    if artifact_kind == "wheel":
        metadata_text = _wheel_metadata_text(artifact_path)
    elif artifact_kind == "sdist":
        metadata_text = _sdist_metadata_text(artifact_path)
    else:
        raise InstalledDistributionVerificationError(
            f"unsupported artifact kind for metadata verification: {artifact_kind!r}"
        )
    return Parser().parsestr(metadata_text).get(field_name)


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
        "--no-cache-dir",
        arguments[0],
    ]
    # The CLI flag applies only to this pip process. Propagate the equivalent
    # environment setting so pip's build-isolation subprocess also avoids any
    # ambient package cache while installing sdist build requirements.
    environment_overrides = {"PIP_NO_CACHE_DIR": "1"}
    if constraint is not None and arguments[0] == "install":
        command.extend(["-c", str(constraint)])
        environment_overrides["PIP_BUILD_CONSTRAINT"] = str(constraint)
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
                "build_provenance": {
                    "workflow": (
                        "prebuilt wheel and sdist supplied by the repository build "
                        "target; artifact metadata and SHA-256 verified here"
                    ),
                    "dist_dir": str(Path(args.dist_dir).resolve()),
                    "verification_python": str(Path(args.python).resolve()),
                },
                "verified": [
                    {
                        "verification_status": "ok",
                        "artifact_kind": report.artifact_kind,
                        "artifact_path": str(report.artifact_path),
                        "artifact_sha256": report.artifact_sha256,
                        "distribution_version": report.distribution_version,
                        "runtime_version": report.runtime_version,
                        "python": report.python_version,
                        "requires_python": report.requires_python,
                        "environment_root": str(report.environment_root),
                        "run_directory": str(report.run_directory),
                        "installation_method": report.installation_method,
                        "phospy_file": str(report.phospy_file),
                        "resource_count": report.resource_count,
                        "ticket_1_boundary_status": report.ticket_1_boundary_status,
                        "constraint_path": (
                            None
                            if report.constraint_path is None
                            else str(report.constraint_path)
                        ),
                        "constraint_sha256": report.constraint_sha256,
                        "probe_results": report.probe_results,
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
