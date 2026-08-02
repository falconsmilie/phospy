from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
    TechnicalReplicatePolicy,
)
from phospy.api.configs import (
    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
)
from phospy.contracts.result_caveats import result_caveats_from_payloads
from phospy.errors import WorkflowValidationError
from phospy.workflows.differential.caveats import (
    DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ("MAPK14", "AKT1", "GSK3B", "RPS6")
_SITES = ("Y182", "T308", "S9", "S235")


def _dataset(sample_ids: tuple[str, ...]) -> AnalysisReadyPhosphoDataset:
    site_index = protein_site_key_index(
        protein_identifiers=_GENES,
        sites=_SITES,
    )
    values = {
        sample_id: [
            10.0 + feature_position + (0.2 * sample_position)
            for feature_position in range(len(site_index))
        ]
        for sample_position, sample_id in enumerate(sample_ids)
    }
    phospho = pd.DataFrame(values, index=site_index.copy())
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(_GENES, _SITES, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": list(_GENES),
            "site": list(_SITES),
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "protein_id": list(_GENES),
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _contrast() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    design: ExperimentalDesign,
    config: DifferentialAnalysisConfig | None = None,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=_contrast(),
        config=DifferentialAnalysisConfig() if config is None else config,
    )


def _one_vs_many_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
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


def test_production_rejects_one_vs_many_single_biological_replicate() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="insufficient replicate counts.*condition='A'.*required=2",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(("A_1", "B_1", "B_2")),
                design=_one_vs_many_design(),
            )
        )


def test_minimum_one_integer_without_exploratory_profile_is_rejected() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="minimum_condition_replicates=1.*exploratory_single_replicate",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(("A_1", "B_1", "B_2")),
                design=_one_vs_many_design(),
                config=DifferentialAnalysisConfig(minimum_condition_replicates=1),
            )
        )


def test_exploratory_profile_allows_one_vs_many_and_serializes_caveat() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(("A_1", "B_1", "B_2")),
            design=_one_vs_many_design(),
            config=DifferentialAnalysisConfig(
                reliability_profile=(
                    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
                ),
            ),
        )
    )

    caveats = {
        caveat.code: caveat
        for caveat in result.caveats
        if caveat.code == DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE
    }
    assert set(caveats) == {DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE}
    caveat = caveats[DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE]
    assert caveat.severity == "warning"
    assert caveat.details["reliability_profile"] == (
        DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
    )
    assert caveat.details["production_supported_inference"] is False
    assert caveat.details["computable_model_output"] is True
    assert caveat.details["contrasted_conditions_below_production_minimum"] == ("A",)

    assert result.policy_provenance is not None
    replicate_policy = result.policy_provenance.replicates
    assert replicate_policy.reliability_profile == (
        DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
    )
    assert replicate_policy.minimum_condition_replicates == 1

    payload = result.to_payload()
    json.dumps(payload, allow_nan=False)
    serialized_caveats = payload["caveats"]
    assert isinstance(serialized_caveats, list)
    serialized = {
        item["code"]: item for item in serialized_caveats if isinstance(item, dict)
    }
    assert (
        serialized[DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE]["details"][
            "inferential_support"
        ]
        == "exploratory_only"
    )
    reconstructed = result_caveats_from_payloads(serialized_caveats)
    reconstructed_by_code = {item.code: item for item in reconstructed}
    assert DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE in (
        reconstructed_by_code
    )
    assert (
        reconstructed_by_code[
            DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE
        ].details["production_supported_inference"]
        is False
    )


def test_technical_replicate_ids_without_biological_ids_do_not_inflate_counts() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_T1",
                condition="A",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A_T2",
                condition="A",
                technical_replicate_id="T2",
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

    with pytest.raises(
        WorkflowValidationError,
        match="technical replicates are not independent biological replicates",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(("A_T1", "A_T2", "B_1", "B_2")),
                design=design,
            )
        )


def test_technical_aggregation_policy_cannot_make_one_biological_replicate_production() -> (
    None
):
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_T1",
                condition="A",
                biological_replicate_id="A_r1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A_T2",
                condition="A",
                biological_replicate_id="A_r1",
                technical_replicate_id="T2",
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

    with pytest.raises(
        WorkflowValidationError,
        match="insufficient replicate counts.*condition='A'.*replicates=1",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(("A_T1", "A_T2", "B_1", "B_2")),
                design=design,
                config=DifferentialAnalysisConfig(
                    technical_replicate_policy=TechnicalReplicatePolicy.MEAN
                ),
            )
        )


def test_production_fixed_block_requires_two_biological_replicates_per_condition() -> (
    None
):
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                block_id="pair_1",
            ),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="insufficient replicate counts.*condition='B'.*required=2",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(("A_1", "B_1")),
                design=design,
                config=DifferentialAnalysisConfig(
                    paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK
                ),
            )
        )


def test_exploratory_fixed_block_still_requires_estimable_model() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                block_id="pair_1",
            ),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="residual degrees of freedom must be positive",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(("A_1", "B_1")),
                design=design,
                config=DifferentialAnalysisConfig(
                    reliability_profile=(
                        DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
                    ),
                    paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
                ),
            )
        )


def test_production_fixed_block_two_biological_replicates_passes() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                block_id="pair_2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                block_id="pair_2",
            ),
        )
    )

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(("A_1", "A_2", "B_1", "B_2")),
            design=design,
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK
            ),
        )
    )

    assert result.policy_provenance is not None
    assert result.policy_provenance.replicates.reliability_profile == "production"
    assert result.policy_provenance.replicates.condition_replicate_counts == (
        ("A", 2),
        ("B", 2),
    )


def test_replicate_reliability_policy_is_not_owned_by_numerical_execution() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    executor_source = (
        repo_root / "src" / "phospy" / "science" / "differential" / "executor.py"
    ).read_text(encoding="utf-8")
    linear_model_source = (
        repo_root / "src" / "phospy" / "science" / "differential" / "linear_model.py"
    ).read_text(encoding="utf-8")

    assert "reliability_profile" not in executor_source
    assert "minimum_condition_replicates" not in executor_source
    assert "reliability_profile" not in linear_model_source
    assert "minimum_condition_replicates" not in linear_model_source
