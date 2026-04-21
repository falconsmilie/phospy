from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.transformation_resolver import (
    DatasetTransformationResolver,
)
from phospy.transformations.models import TransformationState
from phospy.transformations.transformers import IdentityTransformer


def supported_linear_state(*, has_total_matrix: bool) -> TransformationState:
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
        DatasetTransformationResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=total,
        )
        .transformation_state
    )
