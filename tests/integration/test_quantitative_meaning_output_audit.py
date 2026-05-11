from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.io.bundles._shared.processing_state import processing_state_to_payload
from phospy.io.bundles.kinase import (
    KinaseWorkflowConfigSnapshot,
    save_kinase_workflow_bundle,
)
from phospy.io.publishers.workflows import publish_dataset, publish_kinase_workflow
from phospy.provenance.serialization import to_payload as provenance_to_payload

pytestmark = pytest.mark.integration


def _build_log2_dataset(*, corrected: bool) -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, 7.0, 3.0],
            "sample_b": [31.0, 15.0, 7.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "MDFGLCKEGIKDGATMKLCKRERANWQPWQ",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            "sample_a": [3.0, 1.0, 1.0],
            "sample_b": [7.0, 3.0, 1.0],
        },
        index=pd.Index(["MAPK14", "AKT1", "GSK3B"], name="protein_id"),
    )
    preprocessing_config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2", pseudocount=1.0
        ),
        total_protein_correction=(
            DatasetTotalProteinCorrectionConfig(policy="subtract_log_total")
            if corrected
            else DatasetTotalProteinCorrectionConfig(policy="none")
        ),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total if corrected else None,
            organism=Organism.RAT,
            preprocessing_config=preprocessing_config,
        )
    )


def _build_reference_bundle(dataset: AnalysisReadyPhosphoDataset) -> ReferenceBundle:
    site_ids = list(dataset.phospho.index.astype(str))
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    site_ids[0],
                    site_ids[2],
                    site_ids[1],
                    site_ids[2],
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": dataset.site_metadata.loc[:, "site_sequence"].values},
            index=pd.Index(dataset.site_metadata.index.astype(str), name="site_id"),
        ),
    )


def _build_kinase_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=_build_reference_bundle(dataset),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )


def _load_published_dataset_manifest(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    output_root: Path,
) -> dict[str, object]:
    written = publish_dataset(dataset, output_root, output_format="csv")
    return json.loads(written["dataset.manifest"].read_text(encoding="utf-8"))


def test_bundle_manifest_includes_quantitative_meaning_with_numeric_scale(
    tmp_path: Path,
) -> None:
    dataset = _build_log2_dataset(corrected=True)
    request = _build_kinase_request(dataset)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
        output_format="csv",
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))

    intensity_scale_state = manifest["dataset"]["metadata"]["intensity_scale_state"]
    assert intensity_scale_state["phospho"]["kind"] == "log2"
    assert intensity_scale_state["quantity"] == "phospho_total_log_ratio"
    correction_payload = manifest["dataset"]["metadata"]["processing_state"][
        "total_protein_correction"
    ]
    assert correction_payload["output_scale"] == "log2_ratio"
    assert correction_payload["quantitative_meaning"] == "phospho_total_log_ratio"


def test_preprocessing_report_includes_quantitative_meaning_with_numeric_scale() -> (
    None
):
    datasets = (
        (_build_log2_dataset(corrected=False), "phosphosite_log_abundance"),
        (_build_log2_dataset(corrected=True), "phospho_total_log_ratio"),
    )
    for dataset, expected_quantitative_meaning in datasets:
        report = dataset.preprocessing_report
        assert report is not None
        final_rows = report.operations.loc[
            report.operations.loc[:, "stage"] == "final_dataset_construction"
        ]
        assert final_rows.shape[0] == 1
        parameters = final_rows.iloc[0]["parameters"]
        assert parameters["intensity_scale_label"] == "log2"
        assert parameters["quantitative_meaning"] == expected_quantitative_meaning


def test_corrected_matrix_reports_phospho_total_log_ratio(tmp_path: Path) -> None:
    dataset = _build_log2_dataset(corrected=True)

    assert dataset.intensity_scale_state.label == "log2"
    assert dataset.intensity_scale_state.quantity.value == "phospho_total_log_ratio"
    assert (
        dataset.processing_state.total_protein_correction.quantitative_meaning
        == "phospho_total_log_ratio"
    )

    manifest = _load_published_dataset_manifest(
        dataset,
        output_root=tmp_path / "published_corrected",
    )
    assert manifest["intensity_scale"] == "log2"
    assert manifest["quantitative_meaning"] == "phospho_total_log_ratio"
    assert manifest["processing_state"]["intensity_scale"]["quantity"] == (
        "phospho_total_log_ratio"
    )
    assert manifest["processing_state"]["total_protein_correction"]["output_scale"] == (
        "log2_ratio"
    )
    assert manifest["processing_state"]["total_protein_correction"][
        "quantitative_meaning"
    ] == ("phospho_total_log_ratio")


def test_uncorrected_log2_matrix_reports_phosphosite_log_abundance(
    tmp_path: Path,
) -> None:
    dataset = _build_log2_dataset(corrected=False)

    assert dataset.intensity_scale_state.label == "log2"
    assert dataset.intensity_scale_state.quantity.value == "phosphosite_log_abundance"

    manifest = _load_published_dataset_manifest(
        dataset,
        output_root=tmp_path / "published_uncorrected",
    )
    assert manifest["intensity_scale"] == "log2"
    assert manifest["quantitative_meaning"] == "phosphosite_log_abundance"
    assert manifest["processing_state"]["intensity_scale"]["quantity"] == (
        "phosphosite_log_abundance"
    )
    assert (
        manifest["processing_state"]["total_protein_correction"]["output_scale"] is None
    )


def test_dataset_and_provenance_summaries_include_quantitative_meaning_with_numeric_scale() -> (
    None
):
    datasets = (
        (_build_log2_dataset(corrected=False), "phosphosite_log_abundance"),
        (_build_log2_dataset(corrected=True), "phospho_total_log_ratio"),
    )
    for dataset, expected_quantitative_meaning in datasets:
        assert dataset.provenance is not None
        provenance = provenance_to_payload(dataset.provenance)
        workflow_parameters = provenance["workflow_parameters"]
        assert workflow_parameters["intensity_scale_label"] == "log2"
        assert (
            workflow_parameters["quantitative_meaning"] == expected_quantitative_meaning
        )
        total_correction_stage = next(
            (
                stage
                for stage in provenance["preprocessing_stages"]
                if stage["stage"] == "total_protein_correction"
            ),
            None,
        )
        if total_correction_stage is None:
            continue
        diagnostics = total_correction_stage["diagnostics"] or {}
        if diagnostics.get("output_scale") is not None:
            assert diagnostics["quantitative_meaning"] == expected_quantitative_meaning


def test_exported_matrix_metadata_includes_quantitative_meaning_with_numeric_scale(
    tmp_path: Path,
) -> None:
    datasets = (
        (
            "uncorrected",
            _build_log2_dataset(corrected=False),
            "phosphosite_log_abundance",
        ),
        ("corrected", _build_log2_dataset(corrected=True), "phospho_total_log_ratio"),
    )
    for label, dataset, expected_quantitative_meaning in datasets:
        manifest = _load_published_dataset_manifest(
            dataset,
            output_root=tmp_path / f"published_{label}",
        )
        assert manifest["output_format"] == "csv"
        assert manifest["intensity_scale"] == "log2"
        assert manifest["quantitative_meaning"] == expected_quantitative_meaning
        assert manifest["processing_state"]["intensity_scale"]["quantity"] == (
            expected_quantitative_meaning
        )
        correction_payload = manifest["processing_state"]["total_protein_correction"]
        if correction_payload["output_scale"] is not None:
            assert correction_payload["quantitative_meaning"] == (
                expected_quantitative_meaning
            )


def test_strict_quantitative_meaning_export_audit(
    tmp_path: Path,
) -> None:
    lanes = (
        (
            "uncorrected",
            _build_log2_dataset(corrected=False),
            "phosphosite_log_abundance",
        ),
        ("corrected", _build_log2_dataset(corrected=True), "phospho_total_log_ratio"),
    )
    issues: list[str] = []

    for label, dataset, expected_quantitative_meaning in lanes:
        assert dataset.provenance is not None
        assert dataset.preprocessing_report is not None
        assert dataset.intensity_scale_state.label == "log2"
        assert (
            dataset.intensity_scale_state.quantity.value
            == expected_quantitative_meaning
        )

        request = _build_kinase_request(dataset)
        result = KinaseWorkflow().run(request)
        assert result.provenance is not None

        published = publish_kinase_workflow(
            result,
            tmp_path / f"published_{label}",
            output_format="csv",
        )
        published_dataset_manifest = json.loads(
            published["dataset.manifest"].read_text(encoding="utf-8")
        )
        published_kinase_manifest = json.loads(
            published["kinase.manifest"].read_text(encoding="utf-8")
        )

        bundle_root = tmp_path / f"bundle_{label}"
        save_kinase_workflow_bundle(
            result,
            bundle_root,
            config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
            output_format="csv",
        )
        bundle_manifest = json.loads(
            (bundle_root / "manifest.json").read_text(encoding="utf-8")
        )

        representative_payloads: dict[str, object] = {
            f"{label}.publish.dataset_manifest": published_dataset_manifest,
            f"{label}.publish.kinase_manifest": published_kinase_manifest,
            f"{label}.bundle.kinase_manifest": bundle_manifest,
            f"{label}.serialized.processing_state": processing_state_to_payload(
                dataset.processing_state
            ),
            f"{label}.serialized.dataset_provenance": provenance_to_payload(
                dataset.provenance
            ),
            f"{label}.serialized.kinase_provenance": provenance_to_payload(
                result.provenance
            ),
            f"{label}.report.preprocessing.operations": dataset.preprocessing_report.operations.to_dict(
                orient="records"
            ),
            f"{label}.report.preprocessing.row_counts": dataset.preprocessing_report.row_counts.to_dict(
                orient="records"
            ),
        }
        for payload_name, payload in representative_payloads.items():
            issues.extend(
                _find_scale_without_quantitative_meaning_issues(
                    payload,
                    path=payload_name,
                    inherited_quantitative_meaning=None,
                    expected_quantitative_meaning=expected_quantitative_meaning,
                )
            )

    assert issues == [], "quantitative-meaning export audit failed:\n" + "\n".join(
        issues
    )


_SCALE_INDICATOR_KEYS = frozenset(
    (
        "intensity_scale",
        "intensity_scale_label",
        "scale",
        "input_scale",
        "output_scale",
        "kind",
    )
)
_LOG_SCALE_VALUES = frozenset(("log2", "log2_ratio"))
_MEANING_KEYS = ("quantitative_meaning", "quantity")
_IGNORED_SCALE_AUDIT_SUBTREES = frozenset(("index_structure", "column_index_structure"))


def _find_scale_without_quantitative_meaning_issues(
    value: object,
    *,
    path: str,
    inherited_quantitative_meaning: str | None,
    expected_quantitative_meaning: str,
) -> list[str]:
    issues: list[str] = []
    if isinstance(value, Mapping):
        resolved_quantitative_meaning = inherited_quantitative_meaning
        for meaning_key in _MEANING_KEYS:
            candidate = value.get(meaning_key)
            if _is_non_empty_string(candidate):
                resolved_quantitative_meaning = str(candidate).strip()
                break

        for key, raw in value.items():
            if key not in _SCALE_INDICATOR_KEYS:
                continue
            if not _is_non_empty_string(raw):
                continue
            scale_value = str(raw).strip()
            indicator_path = f"{path}.{key}"
            if not _is_non_empty_string(resolved_quantitative_meaning):
                issues.append(
                    f"{indicator_path}: scale metadata '{scale_value}' is missing "
                    "quantitative_meaning/quantity"
                )
                continue
            if (
                scale_value.lower() in _LOG_SCALE_VALUES
                and resolved_quantitative_meaning != expected_quantitative_meaning
            ):
                issues.append(
                    f"{indicator_path}: scale '{scale_value}' expects "
                    f"quantitative_meaning='{expected_quantitative_meaning}' but found "
                    f"'{resolved_quantitative_meaning}'"
                )

        for key, child in value.items():
            if key in _IGNORED_SCALE_AUDIT_SUBTREES:
                continue
            issues.extend(
                _find_scale_without_quantitative_meaning_issues(
                    child,
                    path=f"{path}.{key}",
                    inherited_quantitative_meaning=resolved_quantitative_meaning,
                    expected_quantitative_meaning=expected_quantitative_meaning,
                )
            )
        return issues

    if isinstance(value, list):
        for index, child in enumerate(value):
            issues.extend(
                _find_scale_without_quantitative_meaning_issues(
                    child,
                    path=f"{path}[{index}]",
                    inherited_quantitative_meaning=inherited_quantitative_meaning,
                    expected_quantitative_meaning=expected_quantitative_meaning,
                )
            )
    return issues


def _is_non_empty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())
