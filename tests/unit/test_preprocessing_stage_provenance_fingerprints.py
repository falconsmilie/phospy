from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.advanced import (
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.provenance.hashing import fingerprint_optional_table, hash_table_tolerance
from phospy.provenance.models import PreprocessingStageProvenance, TableFingerprint
from phospy.science.datasets.preprocessing.diagnostics_normalization import (
    _StageDiagnosticsDefaultsResolver,
    _StageDiagnosticsNormalizer,
)
from phospy.science.datasets.preprocessing.fingerprints import (
    _hash_stage_table_fingerprints,
    _StageFingerprintService,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
)
from phospy.science.datasets.preprocessing.trace_builder import _StageTraceBuilder
from phospy.science.transformations.quantitative_contracts import (
    preserve_quantitative_contract,
)


def _base_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [15.0, 7.0],
            "sample_b": [31.0, 15.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _base_site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "protein_id": ["MAPK14", "AKT1"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=index.copy(),
    )


def _base_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"comparison_group": ["treated", "control"]},
        index=columns.copy(),
    )


def _base_total(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            str(columns[0]): [3.0, 1.0],
            str(columns[1]): [7.0, 3.0],
        },
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )


def _build_dataset(
    *,
    phospho: pd.DataFrame | None = None,
    site_metadata: pd.DataFrame | None = None,
    sample_metadata: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
    preprocessing_config: DatasetPreprocessingConfig,
):
    resolved_phospho = _base_phospho() if phospho is None else phospho
    resolved_site_metadata = (
        _base_site_metadata(resolved_phospho.index)
        if site_metadata is None
        else site_metadata
    )
    request_kwargs: dict[str, object] = {
        "phospho": resolved_phospho,
        "site_metadata": resolved_site_metadata,
        "sample_metadata": sample_metadata,
        "total": total,
        "organism": Organism.RAT,
        "preprocessing_config": preprocessing_config,
    }
    policy = str(preprocessing_config.intensity_transform.policy).strip().lower()
    if policy != "log2":
        request_kwargs["input_intensity_scale"] = "linear"
    return AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(**request_kwargs))


def _stage(
    dataset,
    stage_name: str,
) -> PreprocessingStageProvenance:
    assert dataset.provenance is not None
    return next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == stage_name
    )


def _hash_by_name(
    fingerprints: tuple[TableFingerprint, ...],
    table_name: str,
) -> str:
    return next(
        item.tolerance_hash_value for item in fingerprints if item.name == table_name
    )


def test_stage_fingerprint_service_reproduces_existing_exact_and_tolerance_hashes() -> (
    None
):
    phospho = _base_phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_base_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.default(),
    )

    bundle = _StageFingerprintService().run(
        stage_key="fake_stage",
        previous=state,
        current=state,
        consumed_input_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
        produced_output_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
    )

    expected_fingerprint = fingerprint_optional_table(
        phospho,
        name="dataset.phospho",
    )
    assert expected_fingerprint is not None
    assert bundle.consumed_input_tables == (expected_fingerprint,)
    assert bundle.produced_output_tables == (expected_fingerprint,)
    assert bundle.phospho_input_hash == hash_table_tolerance(
        phospho,
        name="fake_stage.input.phospho",
    )
    assert bundle.phospho_output_hash == hash_table_tolerance(
        phospho,
        name="fake_stage.output.phospho",
    )
    assert bundle.input_hash == _hash_stage_table_fingerprints(
        stage_key="fake_stage",
        direction="input",
        table_fingerprints=(expected_fingerprint,),
    )


def test_stage_trace_builder_reproduces_pipeline_trace_payload() -> None:
    class FakeStage:
        stage_key = "fake_stage"

        def validate_before_quantitative_contract(
            self,
            state: PreprocessingState,
        ) -> None:
            del state
            return None

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={"diagnostics": {"policy": "fake"}},
            )

    contract = PreprocessingStageMetadata(
        stage_key="fake_stage",
        display_label="fake_stage",
        operation_name=lambda _plan: "fake_operation",
        serialize_parameters=lambda _plan: {"mode": "test"},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        quantitative_contract=preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )
    phospho = _base_phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_base_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("fake_stage",)),
    )

    _, pipeline_trace = PreprocessingPipeline(
        stage_registry=(FakeStage(),),
        stage_contract_registry=(contract,),
    ).run_with_trace(state)
    stage_result = FakeStage().run(state)
    diagnostics = _StageDiagnosticsNormalizer().run(
        stage_key="fake_stage",
        raw=stage_result.diagnostics,
        defaults=_StageDiagnosticsDefaultsResolver().run(
            previous=state,
            current=stage_result.state,
        ),
    )
    fingerprints = _StageFingerprintService().run(
        stage_key="fake_stage",
        previous=state,
        current=stage_result.state,
        consumed_input_tables=contract.consumed_input_tables,
        produced_output_tables=contract.produced_output_tables,
    )
    trace = _StageTraceBuilder().run(
        stage_key="fake_stage",
        contract=contract,
        interpreted_contract=contract.interpret(state.plan),
        previous=state,
        current=stage_result.state,
        stage_result=stage_result,
        diagnostics=diagnostics,
        fingerprints=fingerprints,
        intensity_transformation_event=None,
    )

    assert trace == pipeline_trace[0]


def test_site_metadata_change_updates_site_matrix_stage_provenance() -> None:
    phospho = _base_phospho()
    site_metadata_a = _base_site_metadata(phospho.index)
    site_metadata_b = site_metadata_a.copy(deep=True)
    site_metadata_b.loc["AKT1;T308;", "site"] = "S473"

    config = DatasetPreprocessingConfig(
        site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
    )
    built_a = _build_dataset(
        phospho=phospho,
        site_metadata=site_metadata_a,
        preprocessing_config=config,
    )
    built_b = _build_dataset(
        phospho=phospho,
        site_metadata=site_metadata_b,
        preprocessing_config=config,
    )

    stage_a = _stage(built_a, "site_matrix")
    stage_b = _stage(built_b, "site_matrix")
    assert _hash_by_name(
        stage_a.consumed_input_tables, "dataset.site_metadata"
    ) != _hash_by_name(stage_b.consumed_input_tables, "dataset.site_metadata")


def test_sample_metadata_change_updates_comparisons_stage_provenance() -> None:
    phospho = _base_phospho()
    sample_metadata_a = _base_sample_metadata(phospho.columns)
    sample_metadata_b = sample_metadata_a.copy(deep=True)
    sample_metadata_b.loc["sample_b", "comparison_group"] = "vehicle"

    config = DatasetPreprocessingConfig(
        comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs")
    )
    built_a = _build_dataset(
        phospho=phospho,
        sample_metadata=sample_metadata_a,
        preprocessing_config=config,
    )
    built_b = _build_dataset(
        phospho=phospho,
        sample_metadata=sample_metadata_b,
        preprocessing_config=config,
    )

    stage_a = _stage(built_a, "comparisons")
    stage_b = _stage(built_b, "comparisons")
    assert _hash_by_name(
        stage_a.consumed_input_tables, "dataset.sample_metadata"
    ) != _hash_by_name(stage_b.consumed_input_tables, "dataset.sample_metadata")
    assert _hash_by_name(
        stage_a.produced_output_tables, "dataset.comparisons"
    ) != _hash_by_name(stage_b.produced_output_tables, "dataset.comparisons")


def test_total_input_change_updates_total_correction_stage_provenance() -> None:
    phospho = _base_phospho()
    total_a = _base_total(phospho.columns)
    total_b = total_a.copy(deep=True)
    total_b.loc["AKT1", "sample_b"] = 4.0
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        ),
        total_protein_correction=DatasetTotalProteinCorrectionConfig(
            policy="subtract_log_total"
        ),
    )
    built_a = _build_dataset(
        phospho=phospho,
        total=total_a,
        preprocessing_config=config,
    )
    built_b = _build_dataset(
        phospho=phospho,
        total=total_b,
        preprocessing_config=config,
    )

    stage_a = _stage(built_a, "total_protein_correction")
    stage_b = _stage(built_b, "total_protein_correction")
    assert _hash_by_name(
        stage_a.consumed_input_tables, "dataset.total"
    ) != _hash_by_name(stage_b.consumed_input_tables, "dataset.total")


def test_stage_configuration_change_updates_stage_provenance() -> None:
    phospho = _base_phospho()
    config_a = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        )
    )
    config_b = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=2.0,
        )
    )
    built_a = _build_dataset(phospho=phospho, preprocessing_config=config_a)
    built_b = _build_dataset(phospho=phospho, preprocessing_config=config_b)

    stage_a = _stage(built_a, "intensity_transform")
    stage_b = _stage(built_b, "intensity_transform")
    assert stage_a.parameters != stage_b.parameters
    assert _hash_by_name(
        stage_a.produced_output_tables, "dataset.phospho"
    ) != _hash_by_name(stage_b.produced_output_tables, "dataset.phospho")


def test_missing_data_missingness_mask_hash_is_stable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, float("nan")],
            "sample_b": [float("nan"), 15.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    config = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=1,
            input_scale="linear",
        )
    )

    built_a = _build_dataset(phospho=phospho, preprocessing_config=config)
    built_b = _build_dataset(phospho=phospho, preprocessing_config=config)

    stage_a = _stage(built_a, "missing_data")
    stage_b = _stage(built_b, "missing_data")
    diagnostics_a = stage_a.diagnostics or {}
    diagnostics_b = stage_b.diagnostics or {}

    assert diagnostics_a.get("missingness_mask_hash") == diagnostics_b.get(
        "missingness_mask_hash"
    )
