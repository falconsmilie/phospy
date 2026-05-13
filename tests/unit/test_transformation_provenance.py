from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.transformations.models import (
    IntensityScaleEstablishmentMode,
    QuantitativeMeaning,
)
from phospy.transformations.transformers import IdentityTransformer


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 2.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )


def test_establishment_provenance_payload_exposes_mode_and_scale() -> None:
    state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
        )
        .intensity_scale_state
    )
    provenance = state.establishment_provenance

    assert provenance is not None
    payload = provenance.to_payload()
    assert payload["scale"] == "linear"
    assert (
        payload["establishment_mode"] == IntensityScaleEstablishmentMode.IDENTITY.value
    )
    assert payload["diagnostic_warnings"] == []


def test_quantitative_meaning_update_preserves_establishment_provenance() -> None:
    state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
        )
        .intensity_scale_state
    )
    updated = state.with_quantitative_meaning(QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE)

    assert updated.establishment_mode is IntensityScaleEstablishmentMode.IDENTITY
    assert updated.establishment_provenance == state.establishment_provenance
