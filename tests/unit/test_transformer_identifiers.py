from __future__ import annotations

import json
from dataclasses import replace

import pandas as pd

from phospy.datasets.construction.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.datasets.construction.executor import DatasetBuildExecutor
from phospy.datasets.construction.intensity_scale import (
    DatasetIntensityScaleResolver,
)
from phospy.preprocessing.contracts import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.preprocessing.contracts.policies import IntensityTransformPolicy
from phospy.quantitative.contracts import TransformationResult
from phospy.quantitative.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.quantitative.transformers import IdentityTransformer, Log2Transformer
from phospy.run_lineage.serialization import to_payload
from phospy.taxonomy import Organism
from tests.support.dataset_preprocessor_fakes import (
    ConformingDatasetPreprocessorFake,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

_BUILT_IN_IDENTITY_IDENTIFIER = (
    "phospy.science.transformations.transformers.identity.IdentityTransformer"
)
_BUILT_IN_LOG2_IDENTIFIER = (
    "phospy.science.transformations.transformers.log2.Log2Transformer"
)
_CUSTOM_IDENTITY_IDENTIFIER = "test_transformer_identifiers.CustomIdentityTransformer"
_CUSTOM_LOG2_IDENTIFIER = "test_transformer_identifiers.CustomLog2Transformer"


class CustomIdentityTransformer(IdentityTransformer):
    pass


class CustomLog2Transformer(Log2Transformer):
    pass


class ExplicitCustomLog2Transformer(Log2Transformer):
    stable_identifier = "user.transformers.explicit-log2"


class StandaloneTransformer:
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
                    established_by="tests.standalone"
                ),
                total=(
                    MatrixIntensityScaleState.log2(established_by="tests.standalone")
                    if total is not None
                    else None
                ),
            ),
        )


def _payload_bytes(payload: object) -> bytes:
    return json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")


def _phospho() -> pd.DataFrame:
    index = site_key_index_from_display_ids(["MAPK14;Y182;"])
    return pd.DataFrame({"sample_a": [3.0]}, index=index)


def _site_metadata() -> pd.DataFrame:
    index = _phospho().index
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;"],
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=index.copy(),
    )


def _declared_linear_state() -> IntensityScaleState:
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(established_by="tests.declared"),
        total=None,
    )


def _resolver_payload(transformer: object) -> tuple[str | None, bytes]:
    resolver = DatasetIntensityScaleResolver(transformer=transformer)  # type: ignore[arg-type]
    if bool(getattr(transformer, "preserves_input_scale_state", False)):
        resolved = resolver.run(
            phospho=_phospho(),
            total=None,
            declared_input_scale_state=_declared_linear_state(),
            declared_input_establishment_mode=(
                IntensityScaleEstablishmentMode.IDENTITY
            ),
        )
    else:
        resolved = resolver.run(phospho=_phospho(), total=None)
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert provenance is not None
    return (
        resolved.intensity_scale_state.established_via,
        _payload_bytes(provenance.to_payload()),
    )


def test_built_in_and_subclass_resolver_identifiers_are_stable() -> None:
    built_in_identity_source, built_in_identity_payload = _resolver_payload(
        IdentityTransformer()
    )
    custom_identity_source, custom_identity_payload = _resolver_payload(
        CustomIdentityTransformer()
    )
    built_in_log2_source, built_in_log2_payload = _resolver_payload(
        Log2Transformer(pseudocount=1.0)
    )
    custom_log2_source, custom_log2_payload = _resolver_payload(
        CustomLog2Transformer(pseudocount=1.0)
    )

    assert built_in_identity_source == _BUILT_IN_IDENTITY_IDENTIFIER
    assert custom_identity_source == _CUSTOM_IDENTITY_IDENTIFIER
    assert built_in_log2_source == _BUILT_IN_LOG2_IDENTIFIER
    assert custom_log2_source == _CUSTOM_LOG2_IDENTIFIER
    assert custom_identity_payload == built_in_identity_payload.replace(
        _BUILT_IN_IDENTITY_IDENTIFIER.encode(),
        _CUSTOM_IDENTITY_IDENTIFIER.encode(),
    )
    assert custom_log2_payload == built_in_log2_payload.replace(
        _BUILT_IN_LOG2_IDENTIFIER.encode(),
        _CUSTOM_LOG2_IDENTIFIER.encode(),
    )


def test_log2_direct_event_and_provenance_bytes_are_stable() -> None:
    built_in = Log2Transformer(pseudocount=1.0).run(phospho=_phospho())
    custom = CustomLog2Transformer(pseudocount=1.0).run(phospho=_phospho())
    built_in_event = built_in.intensity_transformation_event
    custom_event = custom.intensity_transformation_event

    assert built_in_event is not None
    assert custom_event is not None
    assert built_in_event.transformer_name == _BUILT_IN_LOG2_IDENTIFIER
    assert custom_event.transformer_name == _CUSTOM_LOG2_IDENTIFIER
    assert _payload_bytes(custom_event.to_payload()) == _payload_bytes(
        built_in_event.to_payload()
    ).replace(
        _BUILT_IN_LOG2_IDENTIFIER.encode(),
        _CUSTOM_LOG2_IDENTIFIER.encode(),
    )
    assert _payload_bytes(custom.provenance) == _payload_bytes(
        built_in.provenance
    ).replace(
        _BUILT_IN_LOG2_IDENTIFIER.encode(),
        _CUSTOM_LOG2_IDENTIFIER.encode(),
    )


def test_custom_fallback_and_explicit_identifier_precedence_are_preserved() -> None:
    standalone_source, _ = _resolver_payload(StandaloneTransformer())
    assert standalone_source == ("test_transformer_identifiers.StandaloneTransformer")

    transformer = ExplicitCustomLog2Transformer(pseudocount=1.0)
    direct = transformer.run(phospho=_phospho())
    assert direct.intensity_transformation_event is not None
    assert (
        direct.intensity_transformation_event.transformer_name
        == "user.transformers.explicit-log2"
    )
    resolved = DatasetIntensityScaleResolver(transformer=transformer).run(
        phospho=_phospho(),
        total=None,
        establishment_transformer_name="request.explicit-log2",
    )
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert resolved.intensity_scale_state.established_via == (
        "user.transformers.explicit-log2"
    )
    assert provenance is not None
    assert provenance.transformer_name == "request.explicit-log2"


def test_custom_log2_full_dataset_provenance_bytes_are_stable() -> None:
    built_in_bytes = _dataset_provenance_bytes(Log2Transformer(pseudocount=1.0))
    custom_bytes = _dataset_provenance_bytes(CustomLog2Transformer(pseudocount=1.0))

    assert built_in_bytes.count(_BUILT_IN_LOG2_IDENTIFIER.encode()) == 2
    assert custom_bytes == built_in_bytes.replace(
        _BUILT_IN_LOG2_IDENTIFIER.encode(),
        _CUSTOM_LOG2_IDENTIFIER.encode(),
    )


def _dataset_provenance_bytes(transformer: Log2Transformer) -> bytes:
    transformed = transformer.run(phospho=_phospho())
    event = transformed.intensity_transformation_event
    assert event is not None
    event = replace(
        event,
        input_fingerprint="input-fingerprint",
        output_fingerprint="output-fingerprint",
    )
    execution = PreprocessingStageExecution(
        stage=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        operation="log2",
        parameters={"pseudocount": 1.0},
        input_shape=(1, 1),
        output_shape=(1, 1),
        input_hash="input-hash",
        output_hash="output-hash",
        phospho_input_hash="input-fingerprint",
        phospho_output_hash="output-fingerprint",
        intensity_transformation_event=event,
    )
    preprocessed = PreprocessedDatasetBuildTables(
        phospho=transformed.phospho,
        site_metadata=_site_metadata(),
        sample_metadata=None,
        total=None,
        preprocessing_trace=(execution,),
        intensity_transformation_event=event,
    )
    built = DatasetBuildExecutor(
        preprocessor=ConformingDatasetPreprocessorFake(result=preprocessed)
    ).run(
        InterpretedDatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=None,
            total=None,
            organism=Organism.RAT,
            preprocessing_plan=PreprocessingPlan(
                intensity_transform_policy=IntensityTransformPolicy.LOG2,
                intensity_transform_pseudocount=1.0,
                stage_order=(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,),
            ),
        )
    )
    assert built.provenance is not None
    return _payload_bytes(to_payload(built.provenance))
