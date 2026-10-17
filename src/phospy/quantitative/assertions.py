"""Published passive assertion payloads for quantitative state."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.provenance.models import (
    JsonValue,
)
from phospy.quantitative.state_coherence import observe_numeric_domain


def numeric_domain_assertion_payload(
    frame: pd.DataFrame,
    *,
    table_name: str,
) -> Mapping[str, JsonValue]:
    """Return the ordered numeric-domain assertion payload for ``frame``."""

    observation = observe_numeric_domain(frame, table_name=table_name)
    return {
        "table": observation.table_name,
        "observed_numeric_domain": observation.observed_domain.value,
        "value_count": observation.value_count,
        "negative_count": observation.negative_count,
        "zero_count": observation.zero_count,
        "positive_count": observation.positive_count,
        "min": observation.minimum,
        "max": observation.maximum,
    }


__all__ = ["numeric_domain_assertion_payload"]
