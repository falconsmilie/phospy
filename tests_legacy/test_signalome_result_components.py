from __future__ import annotations

from pathlib import Path

import pandas as pd
from phospy.signalomes.serialization import serialize_site_assignments_for_export

from phospy.signalomes import (
    SignalomeCompatibilityView,
    SignalomeCoreResult,
    SignalomeExportAdapter,
    SignalomeFrameBundle,
    SignalomeResult,
    SignalomeVisualizationAdapter,
    build_signalome_result,
)


def _build_component_result() -> SignalomeResult:
    site_ids = [
        "PROTEIN_1;S1;",
        "PROTEIN_1;S2;",
        "PROTEIN_2;S3;",
        "PROTEIN_2;S4;",
    ]
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0, 4.0, 4.0],
            "KINASE_B": [1.1, 1.1, 4.1, 4.1],
            "KINASE_C": [4.0, 4.0, 1.0, 1.0],
        },
        index=site_ids,
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.93, 0.20, 0.25],
            "KINASE_B": [0.90, 0.89, 0.94, 0.92],
            "KINASE_C": [0.10, 0.15, 0.91, 0.90],
        },
        index=site_ids,
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 3.0, 3.1],
            "sample_2": [1.2, 1.0, 2.9, 3.0],
            "sample_3": [0.9, 1.2, 3.2, 2.8],
        },
        index=site_ids,
    )
    return build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        kinase_network_threshold=0.9,
        signalome_cutoff=0.75,
        module_count=2,
    )


def test_signalome_result_exposes_focused_component_adapters() -> None:
    result = _build_component_result()

    assert isinstance(result.core, SignalomeCoreResult)
    assert isinstance(result.frames, SignalomeFrameBundle)
    assert isinstance(result.export, SignalomeExportAdapter)
    assert isinstance(result.visualization, SignalomeVisualizationAdapter)
    assert isinstance(result.compatibility, SignalomeCompatibilityView)


def test_signalome_frame_bundle_matches_legacy_frame_accessors() -> None:
    result = _build_component_result()

    via_bundle = result.frames.to_frames(include_inputs=True)
    via_result = result.to_frames(include_inputs=True)

    assert list(via_bundle) == list(via_result)
    for key in via_bundle:
        pd.testing.assert_frame_equal(via_bundle[key], via_result[key])


def test_signalome_export_adapter_preserves_csv_outputs(tmp_path: Path) -> None:
    result = _build_component_result()

    written = result.export.to_csv(tmp_path)

    assert sorted(written) == [
        "kinase_adjacency_matrix",
        "kinase_correlation_matrix",
        "kinase_module_relationships",
        "kinase_network_edges",
        "kinase_network_nodes",
        "protein_assignments",
        "signalome_modules",
        "site_assignments",
    ]

    reloaded_site_assignments = pd.read_csv(
        written["site_assignments"],
        index_col=0,
    )
    pd.testing.assert_frame_equal(
        reloaded_site_assignments,
        serialize_site_assignments_for_export(result.site_assignments),
    )


def test_signalome_visualization_adapter_matches_result_methods() -> None:
    result = _build_component_result()

    via_adapter_map = result.visualization.to_map_data()
    via_result_map = result.to_map_data()
    pd.testing.assert_frame_equal(via_adapter_map.modules(), via_result_map.modules())
    pd.testing.assert_frame_equal(via_adapter_map.sites(), via_result_map.sites())
    pd.testing.assert_frame_equal(via_adapter_map.kinases(), via_result_map.kinases())
    pd.testing.assert_frame_equal(via_adapter_map.links(), via_result_map.links())

    via_adapter_network = result.visualization.to_network_data()
    via_result_network = result.to_network_data()
    pd.testing.assert_frame_equal(
        via_adapter_network.adjacency(),
        via_result_network.adjacency(),
    )
    pd.testing.assert_frame_equal(
        via_adapter_network.nodes(),
        via_result_network.nodes(),
    )
    pd.testing.assert_frame_equal(
        via_adapter_network.edges(),
        via_result_network.edges(),
    )


def test_signalome_compatibility_view_aliases_match_top_level_properties() -> None:
    result = _build_component_result()
    compatibility = result.compatibility

    pd.testing.assert_frame_equal(
        compatibility.signalome_modules,
        result.signalome_modules,
    )
    pd.testing.assert_frame_equal(
        compatibility.kinase_module_relationships,
        result.kinase_module_relationships,
    )
    pd.testing.assert_frame_equal(
        compatibility.kinase_network_nodes,
        result.kinase_network_nodes,
    )
    assert compatibility.kinase_network == result.kinase_network
