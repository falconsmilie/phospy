from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats

import phospy.science.differential as differential_public
import phospy.science.differential.aggregation as aggregation_public
from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)
from phospy.errors import PhosPyInputError
from phospy.science.differential.aggregation import (
    PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY,
    PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
)
from phospy.science.differential.aggregation.experimental import (
    EXPERIMENTAL_INTERNAL_API,
    EXPERIMENTAL_INTERNAL_REASON,
    WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE,
    PeptideToSiteAggregator,
)
from phospy.science.differential.aggregation.models import (
    PEPTIDE_DIFFERENTIAL_STATISTIC_DISTRIBUTION_MODERATED_T,
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES,
    PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL,
    PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID,
    PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT,
    PeptideDifferentialEstimateTable,
)

ROOT = Path(__file__).resolve().parents[2]

WITHDRAWN_PRODUCTION_EXPORTS = {
    "PeptideDifferentialEstimateTable",
    "PeptideToSiteAggregationConfig",
    "PeptideToSiteAggregationExecutor",
    "PeptideToSiteAggregationResult",
    "PeptideToSiteAggregator",
    "PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY",
    "PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT",
    "PEPTIDE_DIFFERENTIAL_STATISTIC_DISTRIBUTION_MODERATED_T",
    "PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC",
    "PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT",
    "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE",
    "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH",
    "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P",
    "SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES",
    "SUPPORTED_PEPTIDE_TO_SITE_UNCERTAINTY_METHODS",
    "signed_z_from_t_statistic",
    "signed_z_from_two_sided_p_value",
}


def _estimate_frame(*, mapping_policy: str) -> pd.DataFrame:
    statistic = 2.5
    moderated_degrees_of_freedom = 4.0
    return pd.DataFrame(
        {
            "site_id": ["MAPK1;S10;"],
            "peptide_id": ["pep_1"],
            "contrast_id": ["B_vs_A"],
            "contrast_orientation": ["B_minus_A"],
            "effect_scale": ["log2_fold_change"],
            "effect_unit": ["log2_ratio"],
            "model_estimator_id": ["limma_moderated_ols"],
            "statistic_distribution": [
                PEPTIDE_DIFFERENTIAL_STATISTIC_DISTRIBUTION_MODERATED_T
            ],
            "uncertainty_method_version": ["limma_ebayes_moderated_t_v1"],
            "effect": [1.25],
            "standard_error": [0.5],
            "statistic": [statistic],
            "p_value": [
                float(
                    2.0
                    * stats.t.sf(
                        abs(statistic),
                        df=moderated_degrees_of_freedom,
                    )
                )
            ],
            "residual_degrees_of_freedom": [3.0],
            "moderated_degrees_of_freedom": [moderated_degrees_of_freedom],
            "source_experiment_id": ["run_1"],
            "dependence_policy": [
                PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES
            ],
            "peptide_to_site_mapping_policy": [mapping_policy],
        }
    )


def test_posthoc_route_is_absent_from_supported_science_facades() -> None:
    assert WITHDRAWN_PRODUCTION_EXPORTS.isdisjoint(aggregation_public.__all__)
    assert WITHDRAWN_PRODUCTION_EXPORTS.isdisjoint(differential_public.__all__)
    assert not hasattr(aggregation_public, "PeptideToSiteAggregator")
    assert not hasattr(differential_public, "PeptideToSiteAggregator")
    assert PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS == (
        "unsupported_withdrawn_posthoc_estimate_combination_v1"
    )
    assert "supported_typed_estimate_combination" not in (
        PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS
    )


def test_compatibility_shell_is_internal_and_fails_closed() -> None:
    assert EXPERIMENTAL_INTERNAL_API is True
    assert "withdrawn from public support" in EXPERIMENTAL_INTERNAL_REASON
    assert "coherent combined effect/inference" in (
        WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE
    )
    assert "executable mapping semantics" in (
        WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE
    )

    aggregator = PeptideToSiteAggregator()
    assert aggregator.experimental_internal_api is True
    assert aggregator.scientific_support_status == (
        "unsupported_withdrawn_posthoc_estimate_combination_v1"
    )

    with pytest.raises(
        PhosPyInputError,
        match="withdrawn from public support",
    ):
        aggregator.run(pd.DataFrame({"logFC": [-9.75]}), contrast_name="B_vs_A")

    with pytest.raises(
        PhosPyInputError,
        match="coherent combined effect/inference",
    ):
        aggregator.run_estimates(
            estimates=pd.DataFrame({"uncertainty_statistic": [2.411115]}),
        )

    with pytest.raises(
        PhosPyInputError,
        match="executable peptide-to-site mapping semantics",
    ):
        aggregator.run_table(
            peptide_differential_table=pd.DataFrame({"P.Value": [0.015904]}),
            evidence=object(),
        )


@pytest.mark.parametrize(
    "mapping_policy",
    (
        PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT,
        PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL,
    ),
)
def test_non_executable_mapping_policies_are_not_silently_accepted(
    mapping_policy: str,
) -> None:
    with pytest.raises(PhosPyInputError, match="peptide_to_site_mapping_policy"):
        PeptideDifferentialEstimateTable(_estimate_frame(mapping_policy=mapping_policy))


def test_explicit_site_id_remains_parseable_for_internal_future_work_only() -> None:
    table = PeptideDifferentialEstimateTable(
        _estimate_frame(mapping_policy=PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID)
    )

    assert table.site_ids == ("MAPK1;S10;",)


def test_documentation_boundary_no_longer_claims_supported_posthoc_route() -> None:
    docs_text = "\n".join(
        (ROOT / relative_path).read_text(encoding="utf-8")
        for relative_path in (
            "README.md",
            "docs/api/differential-analysis.md",
            "docs/scientific-coverage.md",
            "docs/adr/adr_0020_peptide_evidence_and_site_level_resolution_policy.md",
            ("docs/adr/adr_0041_peptide_to_site_differential_uncertainty_policy.md"),
            "docs/maintenance.md",
            "CHANGELOG.md",
            "docs/release-notes.md",
        )
    )
    normalized = " ".join(docs_text.lower().split())

    assert "supported_typed_estimate_combination_v2" not in normalized
    assert "advanced post-hoc `peptidetositeaggregator` support" not in normalized
    assert "site_result = peptidetositeaggregator().run" not in normalized
    assert "withdrawn from public support" in normalized
    assert "future public support requires" in normalized
    assert "coherent combined estimand" in normalized
    assert "executable mapping semantics" in normalized
    assert "sample-intensity level before differential model fitting" in normalized


def test_installed_distribution_probe_covers_withdrawn_public_boundary() -> None:
    verifier_source = (
        ROOT / "scripts" / "verify_installed_distributions.py"
    ).read_text(encoding="utf-8")

    assert "_verify_withdrawn_peptide_to_site_boundary" in verifier_source
    assert "PeptideToSiteAggregator" in verifier_source
    assert "unsupported_withdrawn_posthoc_estimate_combination_v1" in verifier_source


def test_uncertainty_handling_stays_out_of_validators_interpreters_and_assemblers() -> (
    None
):
    forbidden_files = (
        ROOT / "src" / "phospy" / "validation" / "workflows" / "differential.py",
        ROOT / "src" / "phospy" / "workflows" / "differential" / "interpreter.py",
        ROOT / "src" / "phospy" / "workflows" / "differential" / "result_assembly.py",
    )
    forbidden_fragments = (
        "PeptideDifferentialEstimateTable",
        "PeptideToSiteAggregator",
        "signed_z_from_t_statistic",
        "signed_z_from_two_sided_p_value",
        "stouffer",
    )

    for path in forbidden_files:
        source = path.read_text(encoding="utf-8")
        for fragment in forbidden_fragments:
            assert fragment not in source

    executor_source = (
        ROOT
        / "src"
        / "phospy"
        / "science"
        / "differential"
        / "aggregation"
        / "executor.py"
    ).read_text(encoding="utf-8")
    assert "signed_z_from_t_statistic" in executor_source
    assert "adjust_p_values" in executor_source


def test_preferred_lane_resolves_peptide_intensities_before_core_differential_model() -> (
    None
):
    peptide_evidence = pd.DataFrame(
        {
            "peptide_row_id": ["pep_1", "pep_2"],
            "site_id": ["MAPK1;S10;", "GSK3B;S9;"],
            "unique_feature_id": ["feat_1", "feat_2"],
            "gene_symbol": ["MAPK1", "GSK3B"],
            "protein_accession": ["P28482", "P49841"],
            "site_string": ["S10", "S9"],
            "A_1": [100.0, 80.0],
            "A_2": [105.0, 84.0],
            "B_1": [200.0, 82.0],
            "B_2": [210.0, 83.0],
            "peptide_sequence": ["AAASAAA", "AAASAAA"],
            "modified_peptide_sequence": ["AAA[pS]AAA", "AAA[pS]AAA"],
            "multi_site": [False, False],
            "provenance_source": ["fixture", "fixture"],
            "localisation_confidence": [0.99, 0.98],
            "site_sequence": [
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
        }
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=peptide_evidence,
            peptide_evidence_sample_intensity_columns=("A_1", "A_2", "B_1", "B_2"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )
    assert dataset.provenance is not None
    peptide_resolution = dataset.provenance.workflow_parameters[
        "peptide_evidence_resolution"
    ]
    assert isinstance(peptide_resolution, Mapping)
    assert peptide_resolution["aggregation_policy"] is not None
    assert (
        PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY
        == "sample_intensity_resolution_before_differential_model"
    )

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
            ),
        )
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
        )
    )
    table = result.table_for("B_vs_A")

    assert "B_vs_A" in result.contrast_tables
    assert float(table.loc[table["display_id"] == "MAPK1;S10;", "logFC"].iloc[0]) > 0
