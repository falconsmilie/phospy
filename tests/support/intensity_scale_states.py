from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.preprocessing import build_dataset_processing_state
from phospy.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.transformations.models import IntensityScaleState
from phospy.transformations.transformers import IdentityTransformer


def supported_linear_intensity_scale_state(
    *,
    has_total_matrix: bool,
) -> IntensityScaleState:
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["GENEA;S1;"], name="site_id"),
    )
    total = None
    if has_total_matrix:
        total = pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["GENEA"], name="protein_id"),
        )
    return (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=total,
        )
        .intensity_scale_state
    )


def supported_linear_processing_state(*, has_total_matrix: bool):
    intensity_scale_state = supported_linear_intensity_scale_state(
        has_total_matrix=has_total_matrix
    )
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )
