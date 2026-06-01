from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import (
    InvalidTransformationStateError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
)
from phospy.science.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.references.models import Organism
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.science.transformations.contracts import TransformationResult
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentSource,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    establish_intensity_scale_state,
)
from phospy.science.transformations.transformers import IdentityTransformer
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
)


def _site_key() -> str:
    return encode_site_key(
        build_protein_scoped_site_key(
            organism="rat",
            protein_namespace="protein_id",
            protein_identifier="P28482",
            residue="Y",
            position=182,
            field_name=(
                "tests.unit.test_dataset_transformation_state_establishment.site_key"
            ),
            error_type=ValueError,
        )
    )


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([_site_key()], name="site_key"),
    )


def _total() -> pd.DataFrame:
    return pd.DataFrame({"sample_a": [2.0]}, index=["MAPK14"])


def _site_metadata() -> pd.DataFrame:
    site_key = _site_key()
    return pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": ["MAPK14;Y182;"],
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=pd.Index([site_key], name="site_key"),
    )


def _processing_state_for(intensity_scale_state: IntensityScaleState):
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )


def test_resolver_fails_when_state_is_unknown_and_no_transformer_is_configured() -> (
    None
):
    resolver = DatasetIntensityScaleResolver()

    with pytest.raises(
        TransformationStateEstablishmentError, match="unable to establish"
    ):
        resolver.run(
            phospho=_phospho(),
            total=None,
        )


def test_resolver_uses_configured_transformer_to_establish_state() -> None:
    class DeclaredLog2Transformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState(
                    phospho=MatrixIntensityScaleState.log2(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixIntensityScaleState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetIntensityScaleResolver(transformer=DeclaredLog2Transformer())

    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
    )

    assert resolved.intensity_scale_state.phospho.kind.value == "log2"
    assert resolved.intensity_scale_state.total is not None
    assert resolved.intensity_scale_state.total.kind.value == "log2"
    assert resolved.intensity_scale_state.is_established
    assert resolved.intensity_scale_state.established_via is not None
    assert "DeclaredLog2Transformer" in resolved.intensity_scale_state.established_via


def test_resolver_rejects_invalid_transformer_state_contract() -> None:
    class MismatchedKindTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState(
                    phospho=MatrixIntensityScaleState.linear(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixIntensityScaleState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetIntensityScaleResolver(transformer=MismatchedKindTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="configured transformer produced an invalid intensity scale state",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
        )


def test_resolver_rejects_transformer_that_changes_total_matrix_presence() -> None:
    class DropsTotalTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=None,
                state=IntensityScaleState.raw(has_total_matrix=False),
            )

    resolver = DatasetIntensityScaleResolver(transformer=DropsTotalTransformer())

    with pytest.raises(
        TransformationStateEstablishmentError,
        match="changed total-matrix presence",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
        )


def test_resolver_translates_expected_transformer_contract_errors_with_cause() -> None:
    class FailingTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            raise ValueError("unsupported matrix shape")

    resolver = DatasetIntensityScaleResolver(transformer=FailingTransformer())

    with pytest.raises(
        TransformerExecutionError, match="configured transformer failed"
    ) as exc_info:
        resolver.run(
            phospho=_phospho(),
            total=None,
        )
    assert isinstance(exc_info.value.__cause__, ValueError)


def test_resolver_does_not_translate_unexpected_transformer_errors() -> None:
    class FailingTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            raise RuntimeError("boom")

    resolver = DatasetIntensityScaleResolver(transformer=FailingTransformer())

    with pytest.raises(RuntimeError, match="boom"):
        resolver.run(
            phospho=_phospho(),
            total=None,
        )


def test_dataset_boundary_rejects_declared_intensity_scale_state_bypass() -> None:
    intensity_scale_state = IntensityScaleState.raw(has_total_matrix=False)
    with pytest.raises(
        TransformationValidationError,
        match="must be established through a supported PhosPy path",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            intensity_scale_state=intensity_scale_state,
            processing_state=_processing_state_for(intensity_scale_state),
        )


def test_dataset_boundary_accepts_supported_established_state() -> None:
    declared_linear = IntensityScaleState.raw(has_total_matrix=False)
    supported_state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
            expected_scale_kind=IntensityScaleKind.LINEAR,
            declared_input_scale_state=declared_linear,
            input_declaration_source="tests.unit",
        )
        .intensity_scale_state
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_state,
        processing_state=_processing_state_for(supported_state),
    )
    assert dataset.intensity_scale_state.is_established


def test_direct_mint_established_raw_is_rejected() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        IntensityScaleState.established_raw(has_total_matrix=False)


def test_direct_establishment_function_call_is_rejected() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        establish_intensity_scale_state(
            IntensityScaleState.raw(has_total_matrix=False),
            established_via="phospy.science.datasets.builders.transformation_resolver",
        )


def test_fake_authority_object_is_rejected_even_with_supported_source() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        establish_intensity_scale_state(
            IntensityScaleState.raw(has_total_matrix=False),
            established_via="phospy.science.datasets.builders.transformation_resolver",
            _authority=object(),
        )


def test_dataset_boundary_distinguishes_declared_from_supported_state() -> None:
    declared = IntensityScaleState.raw(has_total_matrix=False)
    supported = supported_linear_intensity_scale_state(has_total_matrix=False)
    assert not declared.is_established
    assert supported.is_established


def test_identity_transformer_is_strict_passthrough_establisher() -> None:
    phospho = _phospho().copy(deep=True)
    total = _total().copy(deep=True)

    result = IdentityTransformer().run(phospho=phospho, total=total)

    pdt.assert_frame_equal(result.phospho, phospho)
    pdt.assert_frame_equal(result.total, total)
    assert not result.state.is_established
    assert result.state.kind.value == "linear"


def test_identity_transformer_cannot_establish_scale_without_explicit_declaration() -> (
    None
):
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="pass-through/identity transformer cannot establish scientific input scale",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
        )


def test_identity_transformer_exposes_scale_capabilities() -> None:
    transformer = IdentityTransformer()
    assert transformer.preserves_input_scale_state is True
    assert transformer.changes_numeric_values is False
    assert transformer.requires_established_input_state is False


def test_identity_transformer_preserves_declared_linear_state() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    declared_linear = IntensityScaleState.raw(has_total_matrix=True)
    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
        expected_scale_kind=IntensityScaleKind.LINEAR,
        declared_input_scale_state=declared_linear,
    )
    assert resolved.intensity_scale_state.kind.value == "linear"
    assert resolved.intensity_scale_state.total is not None
    assert resolved.intensity_scale_state.total.kind.value == "linear"


def test_identity_transformer_preserves_already_declared_log2_state() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    declared_log2 = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="trusted.input"),
        total=MatrixIntensityScaleState.log2(established_by="trusted.input"),
    )
    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
        expected_scale_kind=IntensityScaleKind.LOG2,
        declared_input_scale_state=declared_log2,
    )
    assert resolved.intensity_scale_state.phospho.kind.value == "log2"
    assert resolved.intensity_scale_state.total is not None
    assert resolved.intensity_scale_state.total.kind.value == "log2"
    assert resolved.intensity_scale_state.phospho.established_by == "trusted.input"


def test_capability_preserving_transformer_preserves_declared_state() -> None:
    class CapabilityPreservingTransformer:
        preserves_input_scale_state = True
        changes_numeric_values = False
        requires_established_input_state = False

        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState.raw(has_total_matrix=total is not None),
            )

    resolver = DatasetIntensityScaleResolver(
        transformer=CapabilityPreservingTransformer()
    )
    declared_log2 = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="trusted.input"),
        total=MatrixIntensityScaleState.log2(established_by="trusted.input"),
    )
    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
        expected_scale_kind=IntensityScaleKind.LOG2,
        declared_input_scale_state=declared_log2,
    )

    assert resolved.intensity_scale_state.kind is IntensityScaleKind.LOG2
    assert resolved.intensity_scale_state.total is not None
    assert resolved.intensity_scale_state.total.kind is IntensityScaleKind.LOG2
    assert (
        resolved.intensity_scale_state.establishment_mode
        is IntensityScaleEstablishmentMode.DECLARED
    )


def test_identity_transformer_cannot_establish_log2_from_unknown_state() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="missing intensity state evidence for expected 'log2'",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
            expected_scale_kind=IntensityScaleKind.LOG2,
        )


def test_capability_non_preserving_transformer_rejects_declared_state() -> None:
    class CapabilityNonPreservingTransformer:
        preserves_input_scale_state = False
        changes_numeric_values = True
        requires_established_input_state = False

        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState.raw(has_total_matrix=total is not None),
            )

    resolver = DatasetIntensityScaleResolver(
        transformer=CapabilityNonPreservingTransformer()
    )
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="unsupported identity state establishment",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
            declared_input_scale_state=IntensityScaleState.raw(has_total_matrix=True),
        )


def test_capability_non_preserving_transformer_establishes_new_state() -> None:
    class CapabilityNonPreservingTransformer:
        preserves_input_scale_state = False
        changes_numeric_values = True
        requires_established_input_state = False

        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState(
                    phospho=MatrixIntensityScaleState.log2(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixIntensityScaleState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetIntensityScaleResolver(
        transformer=CapabilityNonPreservingTransformer()
    )
    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
        expected_scale_kind=IntensityScaleKind.LOG2,
    )

    assert resolved.intensity_scale_state.kind is IntensityScaleKind.LOG2
    assert (
        resolved.intensity_scale_state.establishment_mode
        is IntensityScaleEstablishmentMode.DERIVED
    )


def test_resolver_rejects_declared_state_with_non_identity_transformer() -> None:
    class DeclaredLog2Transformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState(
                    phospho=MatrixIntensityScaleState.log2(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixIntensityScaleState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetIntensityScaleResolver(transformer=DeclaredLog2Transformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="unsupported identity state establishment",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
            declared_input_scale_state=IntensityScaleState.raw(has_total_matrix=True),
        )


def test_resolver_distinguishes_mismatched_expected_state() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="mismatched expected intensity state",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
            expected_scale_kind=IntensityScaleKind.LOG2,
            declared_input_scale_state=IntensityScaleState.raw(has_total_matrix=True),
        )


def test_builder_succeeds_when_input_explicitly_declares_already_log2_values() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    assert built.intensity_scale_state.kind is IntensityScaleKind.LOG2


def test_builder_fails_on_expected_log2_without_transform_or_declaration() -> None:
    class NoopPreprocessor:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            site_metadata: pd.DataFrame,
            sample_metadata: pd.DataFrame | None,
            total: pd.DataFrame | None,
            plan: PreprocessingPlan,
        ) -> PreprocessedDatasetBuildTables:
            return PreprocessedDatasetBuildTables(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                preprocessing_trace=None,
            )

    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )
    builder = AnalysisReadyDatasetBuilder(
        executor=DatasetBuildExecutor(preprocessor=NoopPreprocessor())
    )

    with pytest.raises(
        TransformationStateEstablishmentError,
        match="missing intensity state evidence for expected 'log2'",
    ):
        builder.run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(policy="log2")
                ),
            )
        )


def test_bundle_reconstruction_lane_establishes_state() -> None:
    state = intensity_scale_state_from_payload(
        {
            "phospho": {
                "kind": "linear",
                "transformed": False,
                "established_by": "bundle.fixture",
            },
            "total": None,
            "quantity": "phosphosite_abundance",
        }
    )

    assert state.is_established
    assert state.established_via == "phospy.io.bundles._shared.intensity_scale_state"
    provenance = state.establishment_provenance
    assert provenance is not None
    assert (
        provenance.source
        is IntensityScaleEstablishmentSource.RESTORED_FROM_TRUSTED_PROVENANCE
    )


def test_processing_state_payload_rejects_minimal_total_correction_without_explicit_quantitative_meaning() -> (
    None
):
    payload = {
        "intensity_scale": {
            "phospho": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "total": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "quantity": "phospho_total_log_ratio",
        },
        "missing_data": {
            "policy": "forbid",
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": False,
        },
        "normalisation": {"policy": "none"},
        "total_protein_correction": {
            "policy": "subtract_log_total",
            "applied": True,
            "requires_log_scale": True,
            "diagnostics": {
                "diagnostics_schema_version": 1,
                "policy": "subtract_log_total",
                "requested_policy": "subtract_log_total",
                "resolved_policy": "subtract_log_total",
                "quantitative_meaning": "phospho_total_log_ratio",
            },
        },
        "site_matrix": {
            "policy": "as_input",
            "constructed": False,
            "missing_data_policy": "drop_any_missing",
            "minimum_observed_values": None,
            "duplicate_site_policy": "max_mean_signal",
        },
        "comparisons": {
            "policy": "none",
            "sample_group_column": "comparison_group",
            "pairs": None,
        },
    }

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning is required"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_payload_rejects_applied_total_correction_without_versioned_diagnostics() -> (
    None
):
    payload = {
        "intensity_scale": {
            "phospho": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "total": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "quantity": "phospho_total_log_ratio",
        },
        "missing_data": {
            "policy": "forbid",
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": False,
        },
        "normalisation": {"policy": "none"},
        "total_protein_correction": {
            "policy": "subtract_log_total",
            "applied": True,
            "requires_log_scale": True,
            "quantitative_meaning": "phospho_total_log_ratio",
            "diagnostics": None,
        },
        "site_matrix": {
            "policy": "as_input",
            "constructed": False,
            "missing_data_policy": "drop_any_missing",
            "minimum_observed_values": None,
            "duplicate_site_policy": "max_mean_signal",
        },
        "comparisons": {
            "policy": "none",
            "sample_group_column": "comparison_group",
            "pairs": None,
        },
    }

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics "
            "must be an object with"
        ),
    ):
        processing_state_from_payload(payload)
