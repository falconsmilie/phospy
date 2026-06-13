from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy.api.results import (
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.provenance.models import EnvironmentProvenance, RunProvenance
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    ProteinAwarePreparationEligibility,
    ProteinAwareSampleAlignmentDiagnostics,
    ProteinAwareTransformationStateDiagnostics,
)
from phospy.science.datasets.preprocessing.protein_mapping import ProteinMappingStatus


def _sample_alignment() -> ProteinAwareSampleAlignmentDiagnostics:
    return ProteinAwareSampleAlignmentDiagnostics(
        phospho_sample_columns=("sample_a", "sample_b"),
        total_protein_sample_columns=("sample_a", "sample_b"),
        exact_sample_order_match=True,
        sample_order_compatible=True,
        reordered_sample_columns=False,
        allow_reordered_samples=False,
        missing_total_protein_samples=(),
        extra_total_protein_samples=(),
    )


def _transformation_state() -> ProteinAwareTransformationStateDiagnostics:
    return ProteinAwareTransformationStateDiagnostics(
        compatible=True,
        phospho_transformation_state={
            "kind": "log2",
            "transformed": True,
            "established_by": "test.phospho",
        },
        total_protein_transformation_state={
            "kind": "log2",
            "transformed": True,
            "established_by": "test.total",
        },
    )


def _provenance() -> RunProvenance:
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="3.test",
            dependency_versions={},
        ),
        input_tables=(),
        preprocessing_stages=(),
        reference=None,
        workflow_name="protein_aware_preparation",
        workflow_parameters={"policy": "prepare_model_inputs"},
        random_state=None,
        random_seed_policy=None,
        output_tables=(),
    )


def _site_eligibility() -> tuple[ProteinAwareSiteEligibility, ...]:
    return (
        ProteinAwareSiteEligibility(
            site_key="MAPK14;Y182;",
            eligibility=(
                ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
            ),
            mapping_status=ProteinMappingStatus.MATCHED,
            protein_identifier="P53778",
            total_protein_row_key="P53778",
            reasons=("matched_protein_available",),
        ),
        ProteinAwareSiteEligibility(
            site_key="AKT1;T308;",
            eligibility=ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY,
            mapping_status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
            protein_identifier="P31749",
            total_protein_row_key=None,
            reasons=("missing_total_protein_row",),
        ),
        ProteinAwareSiteEligibility(
            site_key="MAPK14;Y182;ambiguous",
            eligibility=ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION,
            mapping_status=ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING,
            protein_identifier=None,
            total_protein_row_key=None,
            reasons=("ambiguous_protein_mapping",),
        ),
    )


def _mapping_diagnostics() -> ProteinAwareMappingDiagnostics:
    return ProteinAwareMappingDiagnostics(
        missing_protein_abundance=pd.DataFrame(
            {
                "site_key": ["AKT1;T308;"],
                "protein_identifier": ["P31749"],
                "mapping_status": ["missing_total_protein_row"],
                "reason": ["missing_total_protein_row"],
            }
        ),
        ambiguous_mapping=pd.DataFrame(
            {
                "site_key": ["MAPK14;Y182;ambiguous"],
                "mapping_status": ["ambiguous_site_protein_mapping"],
                "protein_identifier": [None],
                "candidate_protein_identifiers": [("P53778", "Q16539")],
                "candidate_total_protein_row_keys": [()],
                "reason": ["ambiguous_protein_mapping"],
            }
        ),
    )


def _report() -> ProteinAwarePreparationReport:
    return ProteinAwarePreparationReport(
        site_eligibility=_site_eligibility(),
        mapping_diagnostics=_mapping_diagnostics(),
        sample_alignment=_sample_alignment(),
        transformation_state=_transformation_state(),
        preparation_policy="prepare_model_inputs",
        protein_mapping_policy="allow_missing_with_report",
        policy_parameters={"allow_reordered_samples": False},
        provenance=_provenance(),
    )


def _result() -> ProteinAwarePreparationResult:
    return ProteinAwarePreparationResult(
        matched_pairs=pd.DataFrame(
            {
                "site_key": ["MAPK14;Y182;"],
                "protein_identifier": ["P53778"],
                "total_protein_row_key": ["P53778"],
            }
        ),
        protein_covariate_matrix=pd.DataFrame(
            {"sample_a": [10.0], "sample_b": [11.0]},
            index=pd.Index(["P53778"], name="protein_id"),
        ),
        report=_report(),
    )


def test_protein_aware_preparation_result_construction_with_matched_sites() -> None:
    result = _result()

    assert isinstance(result, ProteinAwarePreparationResult)
    assert result.report.eligible_site_keys == ("MAPK14;Y182;",)
    assert result.report.fallback_site_keys == ("AKT1;T308;",)
    assert result.report.excluded_site_keys == ("MAPK14;Y182;ambiguous",)
    assert result.matched_pairs.loc[0, "total_protein_row_key"] == "P53778"
    assert not hasattr(result, "statistics_table")
    assert not hasattr(result, "contrast_results")


def test_protein_aware_preparation_missing_protein_diagnostics_are_machine_readable() -> (
    None
):
    result = _result()

    diagnostics = result.missing_protein_abundance_diagnostics

    assert diagnostics.to_dict(orient="records") == [
        {
            "site_key": "AKT1;T308;",
            "protein_identifier": "P31749",
            "mapping_status": "missing_total_protein_row",
            "reason": "missing_total_protein_row",
        }
    ]


def test_protein_aware_preparation_ambiguous_mapping_diagnostics_are_machine_readable() -> (
    None
):
    result = _result()

    diagnostics = result.ambiguous_mapping_diagnostics

    assert diagnostics.loc[0, "site_key"] == "MAPK14;Y182;ambiguous"
    assert diagnostics.loc[0, "mapping_status"] == "ambiguous_site_protein_mapping"
    assert diagnostics.loc[0, "candidate_protein_identifiers"] == (
        "P53778",
        "Q16539",
    )


def test_protein_aware_preparation_covariate_matrix_shape_is_aligned() -> None:
    result = _result()

    covariates = result.protein_covariate_matrix

    assert covariates.shape == (1, 2)
    assert covariates.columns.tolist() == ["sample_a", "sample_b"]
    assert covariates.index.tolist() == ["P53778"]


def test_protein_aware_preparation_site_eligibility_table_is_not_warning_text() -> None:
    result = _result()

    table = result.site_eligibility_table

    assert table["site_key"].tolist() == [
        "MAPK14;Y182;",
        "AKT1;T308;",
        "MAPK14;Y182;ambiguous",
    ]
    assert table["eligibility"].tolist() == [
        "eligible_for_protein_aware_preparation",
        "fallback_to_phospho_only",
        "excluded_from_preparation",
    ]
    assert table.loc[1, "reasons"] == ("missing_total_protein_row",)


def test_protein_aware_preparation_serialisation_representation() -> None:
    result = _result()

    payload = result.to_payload()

    assert payload["protein_covariate_matrix_shape"] == [1, 2]
    assert payload["matched_pairs"] == [
        {
            "site_key": "MAPK14;Y182;",
            "protein_identifier": "P53778",
            "total_protein_row_key": "P53778",
        }
    ]
    assert payload["report"]["mapping_diagnostics"]["ambiguous_mapping"][0][
        "candidate_protein_identifiers"
    ] == ["P53778", "Q16539"]
    assert "model" not in payload


def test_protein_aware_preparation_policy_and_provenance_fields() -> None:
    result = _result()

    assert result.preparation_policy == "prepare_model_inputs"
    assert result.protein_mapping_policy == "allow_missing_with_report"
    assert result.provenance is not None
    assert result.provenance.workflow_name == "protein_aware_preparation"
    assert result.to_payload()["provenance"]["workflow_name"] == (
        "protein_aware_preparation"
    )


def test_protein_aware_preparation_public_exports_are_defensive_snapshots() -> None:
    result = _result()

    pairs = result.matched_pairs_dataframe()
    covariates = result.protein_covariate_matrix_dataframe()
    eligibility = result.report.site_eligibility_dataframe()

    pairs.loc[0, "protein_identifier"] = "CHANGED"
    covariates.iloc[0, 0] = 999.0
    eligibility.loc[0, "eligibility"] = "CHANGED"

    assert result.matched_pairs_dataframe().loc[0, "protein_identifier"] == "P53778"
    assert float(result.protein_covariate_matrix_dataframe().iloc[0, 0]) == 10.0
    assert (
        result.report.site_eligibility_dataframe().loc[0, "eligibility"]
        == "eligible_for_protein_aware_preparation"
    )


def test_protein_aware_preparation_mapping_diagnostic_tables_are_exportable() -> None:
    diagnostics = _mapping_diagnostics()

    pdt.assert_frame_equal(
        diagnostics.missing_protein_abundance_dataframe(),
        diagnostics.missing_protein_abundance,
    )
    pdt.assert_frame_equal(
        diagnostics.ambiguous_mapping_dataframe(),
        diagnostics.ambiguous_mapping,
    )
