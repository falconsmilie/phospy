from __future__ import annotations

import inspect
from pathlib import Path

import pandas as pd
import pandas.testing as pdt

from phospy.api import (
    AnalysisReadyPhosphoDataset,
    EnrichmentConfig,
    GeneSetCollection,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
)
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.contracts.results import (
    EnrichmentWorkflowResult,
    PhosphositeImportResult,
    SignalomeWorkflowResult,
)
from phospy.science.activities.models import (
    ActivityMethodMetadata,
    ActivityMethodSummary,
    KinaseActivityResult,
)
from phospy.science.enrichment.models import EnrichmentResultRecord
from phospy.science.prediction.models import (
    KinasePredictionResult,
    KinaseScoringResult,
)
from phospy.science.signalomes.constants import (
    CORRELATION_COLUMN,
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
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_SELECTION_POLICY_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
)
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

ACTIVITY_PAYLOAD_FIELDS = {
    "activity_method_id",
    "activity_method_family",
    "activity_method_label",
    "is_ksea",
    "is_phosr_kinase_activity_equivalent",
}
ACTIVITY_SUMMARY_PAYLOAD_FIELDS = {
    "kinases_evaluated",
    "kinase_condition_pairs_evaluated",
    "kinase_condition_pairs_computed",
    "kinase_condition_pairs_insufficient_substrates",
    "kinase_condition_pairs_invalid_background_variance",
    "kinase_condition_pairs_no_finite_background_values",
    "kinase_condition_pairs_no_finite_substrate_values",
}
ENRICHMENT_RESULT_COLUMNS = [
    "term_id",
    "term_name",
    "collection_kind",
    "identifier_kind",
    "input_overlap_count",
    "background_overlap_count",
    "set_size",
    "overlap_identifiers",
    "p_value",
    "adjusted_p_value",
    "correction_method",
    "enrichment_ratio",
]
RESULT_MODEL_SOURCE_FILES = (
    Path("src/phospy/science/activities/models.py"),
    Path("src/phospy/science/enrichment/models.py"),
    Path("src/phospy/contracts/results/base.py"),
    Path("src/phospy/contracts/results/enrichment.py"),
    Path("src/phospy/contracts/results/kinase.py"),
    Path("src/phospy/contracts/results/signalome.py"),
)
FORBIDDEN_RESULT_MODEL_PATTERNS = (
    ".plot(",
    ".to_csv(",
    ".to_excel(",
    ".to_parquet(",
    ".savefig(",
    "matplotlib",
    "plotly",
    "seaborn",
    "write_text(",
    "write_bytes(",
)
FORBIDDEN_PUBLIC_METHOD_PREFIXES = (
    "export",
    "format",
    "plot",
    "publish",
    "render",
    "save",
    "write",
)


def _activity_result() -> KinaseActivityResult:
    kinase_index = pd.Index(["K1"], name="kinase")
    activity_matrix = pd.DataFrame({"c1": [1.0]}, index=kinase_index, dtype=float)
    count_matrix = pd.DataFrame({"c1": [2]}, index=kinase_index, dtype="int64")
    target_site_key = site_key_index_from_display_ids(["MAPK14;Y182;"])[0]
    statistics_table = pd.DataFrame(
        {
            "kinase": ["K1"],
            "condition": ["c1"],
            "z_score": [1.0],
            "p_value": [0.05],
            "q_value": [0.05],
            "n_substrates": [2],
            "n_background_sites": [4],
            "evidence_threshold": [0.5],
            "evidence_threshold_operator": [">="],
            "evidence_threshold_description": ["unit-test threshold"],
            "min_substrates": [1],
            "computability_status": ["computed"],
            "reason": [""],
        }
    )
    return KinaseActivityResult(
        activity_matrix=activity_matrix,
        p_value_matrix=pd.DataFrame({"c1": [0.05]}, index=kinase_index, dtype=float),
        q_value_matrix=pd.DataFrame({"c1": [0.05]}, index=kinase_index, dtype=float),
        confidence_interval_low=pd.DataFrame(
            {"c1": [0.5]},
            index=kinase_index,
            dtype=float,
        ),
        confidence_interval_high=pd.DataFrame(
            {"c1": [1.5]},
            index=kinase_index,
            dtype=float,
        ),
        substrate_count_matrix=count_matrix,
        thresholded_substrate_mean_activity=pd.DataFrame(
            {"c1": [2.0]},
            index=kinase_index,
            dtype=float,
        ),
        thresholded_substrate_counts=pd.Series(
            [2],
            index=kinase_index,
            name="n_substrates",
            dtype="int64",
        ),
        target_counts=pd.Series(
            [2],
            index=kinase_index,
            name="n_targets",
            dtype="int64",
        ),
        target_table=pd.DataFrame(
            {
                "site_id": [str(target_site_key)],
                "site_key": [str(target_site_key)],
                "display_id": ["MAPK14;Y182;"],
                "kinase": ["K1"],
                "score": [0.9],
            }
        ),
        statistics_table=statistics_table,
    )


def _gene_collection() -> GeneSetCollection:
    return GeneSetCollection(
        sets={"mapk_pathway": ("AKT1", "MAPK1")},
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        term_names={"mapk_pathway": "MAPK pathway"},
    )


def _enrichment_result() -> EnrichmentWorkflowResult:
    record = EnrichmentResultRecord(
        term_id="mapk_pathway",
        term_name="MAPK pathway",
        collection_kind="gene_set",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        input_overlap_count=1,
        background_overlap_count=2,
        set_size=2,
        overlap_identifiers=("AKT1",),
        p_value=0.5,
        adjusted_p_value=0.5,
        correction_method="benjamini_hochberg",
        enrichment_ratio=1.5,
    )
    return EnrichmentWorkflowResult(
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        config=EnrichmentConfig(),
        records=(record,),
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["P1;S1;"]
    site_index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame({"sample_a": [1.0]}, index=site_index)
    site_metadata = pd.DataFrame(
        {
            SITE_KEY_COLUMN: site_index.astype(str).tolist(),
            DISPLAY_ID_COLUMN: display_ids,
            **site_key_context_columns(site_index),
            GENE_SYMBOL_COLUMN: ["P1"],
            SITE_COLUMN: ["S1"],
            "site_sequence": [("A" * 15) + "S" + ("A" * 15)],
            PROTEIN_COLUMN: ["P1"],
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False,
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _kinase_result(dataset: AnalysisReadyPhosphoDataset) -> KinaseWorkflowResult:
    site_key = str(dataset.phospho.index[0])
    score_matrix = pd.DataFrame(
        {"K1": [0.8]},
        index=pd.Index([site_key], name=SITE_KEY_COLUMN),
        dtype=float,
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["K1"], "substrate_site": ["P1;S1;"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": [("A" * 15) + "S" + ("A" * 15)]},
                index=pd.Index(["P1;S1;"], name="site_id"),
            ),
        ),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=score_matrix),
    )


def _empty_signalome_assignments() -> SignalomeAssignments:
    columns = [
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
    ]
    return SignalomeAssignments(
        pd.DataFrame(columns=columns, index=pd.Index([], name=SITE_KEY_COLUMN)),
    )


def _signalome_result() -> SignalomeWorkflowResult:
    dataset = _dataset()
    site_key = str(dataset.phospho.index[0])
    return SignalomeWorkflowResult(
        dataset=dataset,
        kinase_result=_kinase_result(dataset),
        module_assignments=_empty_signalome_assignments(),
        signalome_modules=SignalomeModules(
            pd.DataFrame(
                index=pd.Index([], name=MODULE_ID_COLUMN, dtype="int64"),
                columns=pd.Index(["K1"], name="kinase"),
                dtype=float,
            )
        ),
        kinase_network=KinaseNetwork(
            edges=pd.DataFrame(
                columns=[SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN, CORRELATION_COLUMN]
            )
        ),
        expanded_signalome=pd.DataFrame(
            {
                SITE_KEY_COLUMN: [site_key],
                DISPLAY_ID_COLUMN: ["P1;S1;"],
                "kinase": ["K1"],
            }
        ),
    )


def test_activity_payload_helpers_return_expected_fields_and_fresh_dicts() -> None:
    metadata = ActivityMethodMetadata(
        activity_method_id="method_id",
        activity_method_family="family",
        activity_method_label="label",
        is_ksea=False,
        is_phosr_kinase_activity_equivalent=False,
    )
    payload = metadata.to_payload()

    assert set(payload) == ACTIVITY_PAYLOAD_FIELDS
    payload["activity_method_id"] = "mutated"
    assert metadata.to_payload()["activity_method_id"] == "method_id"

    summary = ActivityMethodSummary(
        kinases_evaluated=1,
        kinase_condition_pairs_evaluated=2,
        kinase_condition_pairs_computed=3,
        kinase_condition_pairs_insufficient_substrates=4,
        kinase_condition_pairs_invalid_background_variance=5,
        kinase_condition_pairs_no_finite_background_values=6,
        kinase_condition_pairs_no_finite_substrate_values=7,
    )
    summary_payload = summary.to_payload()

    assert set(summary_payload) == ACTIVITY_SUMMARY_PAYLOAD_FIELDS
    summary_payload["kinases_evaluated"] = 99
    assert summary.to_payload()["kinases_evaluated"] == 1


def test_activity_dataframe_helpers_are_expected_field_snapshots() -> None:
    result = _activity_result()
    activity_snapshot = result.to_dataframe()
    p_value_snapshot = result.p_value_matrix_dataframe()
    q_value_snapshot = result.q_value_matrix_dataframe()
    low_snapshot = result.confidence_interval_low_dataframe()
    high_snapshot = result.confidence_interval_high_dataframe()
    count_snapshot = result.substrate_count_matrix_dataframe()
    threshold_mean_snapshot = result.thresholded_substrate_mean_activity_dataframe()
    target_snapshot = result.target_table_dataframe()
    statistics_snapshot = result.statistics_table_dataframe()

    assert list(activity_snapshot.columns) == ["c1"]
    assert list(count_snapshot.columns) == ["c1"]
    assert list(target_snapshot.columns) == [
        "site_id",
        "site_key",
        "display_id",
        "kinase",
        "score",
    ]
    assert statistics_snapshot is not None
    assert p_value_snapshot is not None
    assert q_value_snapshot is not None
    assert low_snapshot is not None
    assert high_snapshot is not None

    activity_snapshot.iloc[0, 0] = 99.0
    p_value_snapshot.iloc[0, 0] = 0.99
    q_value_snapshot.iloc[0, 0] = 0.99
    low_snapshot.iloc[0, 0] = -99.0
    high_snapshot.iloc[0, 0] = 99.0
    count_snapshot.iloc[0, 0] = 99
    threshold_mean_snapshot.iloc[0, 0] = 99.0
    target_snapshot.loc[0, "score"] = 0.1
    statistics_snapshot.loc[0, "reason"] = "mutated"

    assert result.activity_matrix.iloc[0, 0] == 1.0
    assert result.p_value_matrix is not None
    assert result.p_value_matrix.iloc[0, 0] == 0.05
    assert result.q_value_matrix is not None
    assert result.q_value_matrix.iloc[0, 0] == 0.05
    assert result.confidence_interval_low is not None
    assert result.confidence_interval_low.iloc[0, 0] == 0.5
    assert result.confidence_interval_high is not None
    assert result.confidence_interval_high.iloc[0, 0] == 1.5
    assert result.substrate_count_matrix.iloc[0, 0] == 2
    assert result.thresholded_substrate_mean_activity.iloc[0, 0] == 2.0
    assert result.target_table.loc[0, "score"] == 0.9
    assert result.statistics_table is not None
    assert result.statistics_table.loc[0, "reason"] == ""

    semantics = result.count_field_semantics
    semantics["target_counts"] = "mutated"
    assert result.count_field_semantics["target_counts"] != "mutated"


def test_enrichment_collection_mappings_are_fresh_snapshots() -> None:
    collection = _gene_collection()
    sets_snapshot = collection.sets
    names_snapshot = collection.term_names
    members_snapshot = collection.members_by_set_id

    sets_snapshot["new"] = ("BAD",)
    names_snapshot["mapk_pathway"] = "mutated"
    members_snapshot["mapk_pathway"] = ("BAD",)

    assert "new" not in collection.sets
    assert collection.term_names["mapk_pathway"] == "MAPK pathway"
    assert collection.members_by_set_id["mapk_pathway"] == ("AKT1", "MAPK1")


def test_enrichment_result_dataframe_helper_is_expected_field_snapshot() -> None:
    result = _enrichment_result()
    table = result.to_dataframe()

    assert list(table.columns) == ENRICHMENT_RESULT_COLUMNS
    table.loc[0, "term_id"] = "mutated"
    table.loc[0, "overlap_identifiers"] = ("BAD",)

    assert result.table.loc[0, "term_id"] == "mapk_pathway"
    assert result.result_table.loc[0, "overlap_identifiers"] == ("AKT1",)


def test_phosphosite_import_result_helpers_are_snapshots() -> None:
    result = PhosphositeImportResult(
        phospho_matrix_candidate=pd.DataFrame({"sample_a": [1.0]}, index=["site_1"]),
        site_metadata_candidate=pd.DataFrame({"site_id": ["site_1"]}, index=["site_1"]),
        peptide_evidence=pd.DataFrame({"sample_a": [1.0]}, index=["peptide_1"]),
        sample_column_mapping={"Intensity A": "sample_a"},
    )
    phospho = result.phospho_matrix_candidate
    site_metadata = result.site_metadata_candidate
    evidence = result.peptide_evidence
    mapping = result.sample_column_mapping

    assert evidence is not None
    phospho.iloc[0, 0] = 99.0
    site_metadata.iloc[0, 0] = "mutated"
    evidence.iloc[0, 0] = 99.0
    mapping["new"] = "sample_b"

    assert result.phospho_matrix_candidate.iloc[0, 0] == 1.0
    assert result.site_metadata_candidate.iloc[0, 0] == "site_1"
    assert result.peptide_evidence is not None
    assert result.peptide_evidence.iloc[0, 0] == 1.0
    assert result.sample_column_mapping == {"Intensity A": "sample_a"}


def test_signalome_result_dataframe_helper_is_snapshot() -> None:
    result = _signalome_result()
    expanded = result.to_dataframe()

    assert expanded is not None
    assert list(expanded.columns) == [SITE_KEY_COLUMN, DISPLAY_ID_COLUMN, "kinase"]
    expanded.loc[0, "kinase"] = "mutated"

    assert result.expanded_signalome is not None
    assert result.expanded_signalome.loc[0, "kinase"] == "K1"


def test_primary_result_model_modules_do_not_export_plot_or_report() -> None:
    for source_path in RESULT_MODEL_SOURCE_FILES:
        source = source_path.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_RESULT_MODEL_PATTERNS:
            assert pattern not in source, f"{source_path} contains {pattern}"


def test_result_model_public_helpers_keep_passive_responsibility_names() -> None:
    result_classes = (
        PhosphositeImportResult,
        EnrichmentWorkflowResult,
        SignalomeWorkflowResult,
        ActivityMethodMetadata,
        ActivityMethodSummary,
        KinaseActivityResult,
    )

    for result_class in result_classes:
        helper_names = {
            name
            for name, value in inspect.getmembers(result_class)
            if not name.startswith("_")
            and (
                inspect.isfunction(value)
                or isinstance(value, property)
                or isinstance(value, classmethod)
            )
        }
        forbidden = {
            name
            for name in helper_names
            if name.startswith(FORBIDDEN_PUBLIC_METHOD_PREFIXES)
        }
        assert forbidden == set()


def test_enrichment_result_table_aliases_return_equal_independent_snapshots() -> None:
    result = _enrichment_result()
    table = result.table
    result_table = result.result_table

    pdt.assert_frame_equal(table, result_table)
    table.loc[0, "term_name"] = "mutated"
    assert result.result_table.loc[0, "term_name"] == "MAPK pathway"
