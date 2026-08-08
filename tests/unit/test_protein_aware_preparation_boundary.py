from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt

from phospy import AnalysisReadyDatasetBuilder
from phospy.advanced import DatasetProteinAwarePreparationConfig
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)

ROOT = Path(__file__).resolve().parents[2]
_DIFFERENTIAL_CODE_DIRS = (
    ROOT / "src" / "phospy" / "workflows" / "differential",
    ROOT / "src" / "phospy" / "science" / "differential",
)
_FORBIDDEN_DIFFERENTIAL_PREPARATION_TOKENS = (
    "ProteinAwarePreparationResult",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationStage",
    "protein_aware_preparation",
    "protein_covariate_matrix",
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [1.0, 1.1],
            "A_2": [1.2, 1.3],
            "B_1": [2.0, 2.1],
            "B_2": [2.2, 2.3],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "protein_id": ["P53778", "P31749"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.96],
        },
        index=_phospho().index.copy(),
    )


def _total(*, include_akt1: bool = True) -> pd.DataFrame:
    index = ["P53778", "P31749"] if include_akt1 else ["P53778"]
    rows = {
        "P53778": [10.0, 11.0, 12.0, 13.0],
        "P31749": [20.0, 21.0, 22.0, 23.0],
    }
    return pd.DataFrame(
        {
            sample_id: [rows[row_key][position] for row_key in index]
            for position, sample_id in enumerate(("A_1", "A_2", "B_1", "B_2"))
        },
        index=pd.Index(index, name="protein_id"),
    )


def _build_dataset(
    *,
    preprocessing_config: DatasetPreprocessingConfig,
    total: pd.DataFrame | None = None,
):
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=_total() if total is None else total,
            organism=Organism.HUMAN,
            input_intensity_scale="log2",
            preprocessing_config=preprocessing_config,
        )
    )


def test_protein_aware_preparation_output_is_separate_from_analysis_ready_dataset() -> (
    None
):
    disabled = _build_dataset(preprocessing_config=DatasetPreprocessingConfig())
    prepared = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs",
                protein_mapping_policy="require_unambiguous",
            )
        )
    )

    assert prepared.protein_aware_preparation is not None
    assert prepared.preprocessing_report is not None
    assert (
        prepared.preprocessing_report.protein_aware_preparation
        is prepared.protein_aware_preparation.report
    )
    pdt.assert_frame_equal(prepared.phospho, disabled.phospho)
    pdt.assert_frame_equal(prepared.site_metadata, disabled.site_metadata)
    assert (
        prepared.protein_aware_preparation.matched_pairs_dataframe()
        .loc[:, "site_key"]
        .tolist()
        == prepared.phospho.index.astype(str).tolist()
    )
    covariate_index = prepared.protein_aware_preparation.protein_covariate_matrix_dataframe().index.tolist()
    assert covariate_index == [
        "P53778",
        "P31749",
    ]


def test_protein_aware_diagnostics_are_retained_apart_from_phospho_matrix() -> None:
    disabled = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(),
        total=_total(include_akt1=False),
    )
    prepared = _build_dataset(
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs",
                protein_mapping_policy="allow_missing_with_report",
            )
        ),
        total=_total(include_akt1=False),
    )

    assert prepared.protein_aware_preparation is not None
    diagnostics = (
        prepared.protein_aware_preparation.missing_protein_abundance_diagnostics
    )
    assert diagnostics.to_dict(orient="records") == [
        {
            "site_key": prepared.protein_aware_preparation.report.fallback_site_keys[0],
            "protein_identifier": "P31749",
            "mapping_status": "missing_total_protein_row",
            "reason": "missing_total_protein_row",
        }
    ]
    assert prepared.preprocessing_report is not None
    assert (
        prepared.preprocessing_report.protein_aware_preparation
        is prepared.protein_aware_preparation.report
    )
    covariate_index = prepared.protein_aware_preparation.protein_covariate_matrix_dataframe().index.tolist()
    assert covariate_index == ["P53778"]
    pdt.assert_frame_equal(prepared.phospho, disabled.phospho)
    for diagnostic_column in ("mapping_status", "reason", "protein_covariate_matrix"):
        assert diagnostic_column not in prepared.phospho.columns


def test_differential_domains_do_not_import_or_own_protein_aware_preparation() -> None:
    violations: list[str] = []
    for directory in _DIFFERENTIAL_CODE_DIRS:
        for path in directory.rglob("*.py"):
            source = path.read_text(encoding="utf-8")
            for token in _FORBIDDEN_DIFFERENTIAL_PREPARATION_TOKENS:
                if token not in source:
                    continue
                relative_path = path.relative_to(ROOT).as_posix()
                violations.append(f"{relative_path}: contains {token!r}")

    assert not violations, (
        "protein-aware preparation must stay in dataset preprocessing/building "
        "domains and must not move into differential workflow or result code:\n"
        + "\n".join(violations)
    )
