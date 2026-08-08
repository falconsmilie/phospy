from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._kinase.constants import (
    KINASE_BUNDLE_KIND,
    KINASE_BUNDLE_MANIFEST_VERSION,
)
from phospy.io.bundles._kinase.manifest import (
    KinaseManifestSections,
)
from phospy.io.bundles._kinase.manifest import (
    parse_manifest as parse_kinase_manifest,
)
from phospy.io.bundles._signalome.constants import (
    SIGNALOME_BUNDLE_KIND,
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
)
from phospy.io.bundles._signalome.manifest import (
    SignalomeManifestSections,
)
from phospy.io.bundles._signalome.manifest import (
    parse_manifest as parse_signalome_manifest,
)

ManifestParser = Callable[[dict[str, object]], object]


def _table_entry(path: str) -> dict[str, object]:
    return {
        "path": path,
        "sha256": "0" * 64,
        "byte_size": 0,
        "logical_type": "table",
        "shape": {"rows": 0, "columns": 0},
    }


def _file_entry(path: str) -> dict[str, object]:
    return {
        "path": path,
        "sha256": "0" * 64,
        "byte_size": 0,
        "logical_type": "config_snapshot",
    }


def _kinase_manifest_payload() -> dict[str, object]:
    return {
        "bundle_type": KINASE_BUNDLE_KIND,
        "manifest_version": KINASE_BUNDLE_MANIFEST_VERSION,
        "table_format": "csv",
        "dataset": {
            "metadata": {},
            "tables": {
                "phospho": _table_entry("tables/phospho.csv"),
                "site_metadata": _table_entry("tables/site_metadata.csv"),
                "sample_metadata": None,
                "total": None,
            },
        },
        "resolved_references": {
            "metadata": {},
            "tables": {
                "kinase_substrate_map": _table_entry(
                    "references/kinase_substrate_map.csv"
                ),
                "site_sequences": _table_entry("references/site_sequences.csv"),
            },
        },
        "outputs": {
            "scoring": {
                "tables": {
                    "profile_scores": _table_entry("outputs/profile_scores.csv"),
                    "motif_scores": None,
                    "rank_weighted_fusion_scores": None,
                    "score_fusion_weights": None,
                },
            },
            "prediction": {
                "tables": {
                    "pred_mat": _table_entry("outputs/pred_mat.csv"),
                    "substrate_list": None,
                },
            },
            "activity": {
                "enabled": False,
                "method": None,
                "summary": None,
                "input_semantics": None,
                "profile_metadata": None,
                "membership_selection": None,
                "tables": {
                    "weighted_activity": None,
                    "thresholded_substrate_mean_activity": None,
                    "thresholded_substrate_counts": None,
                    "activity_substrate_counts": None,
                    "target_counts": None,
                    "target_table": None,
                    "statistics_table": None,
                },
            },
        },
        "provenance": {},
        "config_snapshot": _file_entry("config/config.json"),
    }


def _signalome_manifest_payload() -> dict[str, object]:
    return {
        "bundle_type": SIGNALOME_BUNDLE_KIND,
        "manifest_version": SIGNALOME_BUNDLE_MANIFEST_VERSION,
        "table_format": "csv",
        "dataset": {
            "metadata": {},
            "tables": {
                "phospho": _table_entry("tables/phospho.csv"),
                "site_metadata": _table_entry("tables/site_metadata.csv"),
                "sample_metadata": None,
                "total": None,
            },
        },
        "resolved_references": {
            "metadata": {},
            "tables": {
                "kinase_substrate_map": _table_entry(
                    "references/kinase_substrate_map.csv"
                ),
                "site_sequences": _table_entry("references/site_sequences.csv"),
            },
        },
        "upstream_kinase_outputs": {
            "scoring": {
                "tables": {
                    "profile_scores": _table_entry("upstream/profile_scores.csv"),
                    "motif_scores": None,
                    "rank_weighted_fusion_scores": None,
                    "score_fusion_weights": None,
                },
            },
            "prediction": {
                "tables": {
                    "pred_mat": _table_entry("upstream/pred_mat.csv"),
                    "substrate_list": None,
                },
            },
            "activity": {
                "enabled": False,
                "tables": {
                    "weighted_activity": None,
                    "thresholded_substrate_mean_activity": None,
                    "thresholded_substrate_counts": None,
                    "activity_substrate_counts": None,
                    "target_counts": None,
                    "target_table": None,
                },
            },
        },
        "signalome_outputs": {
            "metadata": {
                "kinase_network_nodes_present": False,
                "expanded_signalome_present": False,
                "module_selection_diagnostics": {},
                "score_preconditioning_diagnostics": {},
                "network_correlation_diagnostics": {},
            },
            "tables": {
                "module_assignments": _table_entry("signalome/module_assignments.csv"),
                "signalome_modules": _table_entry("signalome/signalome_modules.csv"),
                "kinase_network_edges": _table_entry(
                    "signalome/kinase_network_edges.csv"
                ),
                "kinase_network_nodes": None,
                "kinase_network_candidate_correlations": None,
                "expanded_signalome": None,
                "site_membership": None,
                "protein_site_context": None,
            },
        },
        "provenance": {},
        "config_snapshot": _file_entry("config/config.json"),
    }


@pytest.mark.parametrize(
    ("payload_factory", "parser", "section_type", "schema_key"),
    [
        (
            _kinase_manifest_payload,
            parse_kinase_manifest,
            KinaseManifestSections,
            "outputs",
        ),
        (
            _signalome_manifest_payload,
            parse_signalome_manifest,
            SignalomeManifestSections,
            "signalome_outputs",
        ),
    ],
)
def test_manifest_parser_stages_accept_minimal_current_payload(
    payload_factory: Callable[[], dict[str, object]],
    parser: ManifestParser,
    section_type: type[object],
    schema_key: str,
) -> None:
    del schema_key

    assert isinstance(parser(payload_factory()), section_type)


@pytest.mark.parametrize(
    ("payload_factory", "parser", "schema_key"),
    [
        (_kinase_manifest_payload, parse_kinase_manifest, "outputs"),
        (_signalome_manifest_payload, parse_signalome_manifest, "signalome_outputs"),
    ],
)
def test_manifest_schema_stage_rejects_missing_required_object(
    payload_factory: Callable[[], dict[str, object]],
    parser: ManifestParser,
    schema_key: str,
) -> None:
    payload = payload_factory()
    payload.pop(schema_key)

    with pytest.raises(PhosPyInputError, match="missing required field"):
        parser(payload)


@pytest.mark.parametrize(
    ("payload_factory", "parser"),
    [
        (_kinase_manifest_payload, parse_kinase_manifest),
        (_signalome_manifest_payload, parse_signalome_manifest),
    ],
)
def test_manifest_path_stage_rejects_parent_directory_traversal(
    payload_factory: Callable[[], dict[str, object]],
    parser: ManifestParser,
) -> None:
    payload = payload_factory()
    mutated = deepcopy(payload)
    dataset = mutated["dataset"]
    assert isinstance(dataset, dict)
    tables = dataset["tables"]
    assert isinstance(tables, dict)
    phospho = tables["phospho"]
    assert isinstance(phospho, dict)
    phospho["path"] = "../outside.csv"

    with pytest.raises(PhosPyInputError, match="parent-directory traversal"):
        parser(mutated)


@pytest.mark.parametrize(
    ("payload_factory", "parser"),
    [
        (_kinase_manifest_payload, parse_kinase_manifest),
        (_signalome_manifest_payload, parse_signalome_manifest),
    ],
)
def test_manifest_file_record_stage_rejects_invalid_digest(
    payload_factory: Callable[[], dict[str, object]],
    parser: ManifestParser,
) -> None:
    payload = payload_factory()
    mutated = deepcopy(payload)
    dataset = mutated["dataset"]
    assert isinstance(dataset, dict)
    tables = dataset["tables"]
    assert isinstance(tables, dict)
    phospho = tables["phospho"]
    assert isinstance(phospho, dict)
    phospho["sha256"] = "not-a-sha"

    with pytest.raises(PhosPyInputError, match="64-character SHA-256"):
        parser(mutated)


@pytest.mark.parametrize(
    ("payload_factory", "parser"),
    [
        (_kinase_manifest_payload, parse_kinase_manifest),
        (_signalome_manifest_payload, parse_signalome_manifest),
    ],
)
def test_manifest_model_assembly_stage_rejects_invalid_metadata_object(
    payload_factory: Callable[[], dict[str, object]],
    parser: ManifestParser,
) -> None:
    payload = payload_factory()
    mutated = deepcopy(payload)
    dataset = mutated["dataset"]
    assert isinstance(dataset, dict)
    dataset["metadata"] = "invalid"

    with pytest.raises(PhosPyInputError, match="bundle manifest.dataset.metadata"):
        parser(mutated)
