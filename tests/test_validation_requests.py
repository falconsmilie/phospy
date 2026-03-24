from __future__ import annotations

import pytest

from phospy.validation.errors import RequestValidationError
from phospy.validation.requests import (
    CorePipelineRequest,
)


def test_core_pipeline_request_requires_existing_paths(tmp_path) -> None:
    total_path = tmp_path / "missing_total.tsv"
    phospho_path = tmp_path / "missing_phospho.tsv"

    with pytest.raises(RequestValidationError, match="Path does not exist"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
        )


def test_core_pipeline_request_rejects_duplicate_comparisons(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    request = CorePipelineRequest.validate_request(
        total_path=total_path,
        phospho_path=phospho_path,
        kinase_activity_threshold=0.7,
        kinase_activity_min_substrates=5,
        kinase_activity_top_n_substrates=12,
    )

    assert request.kinase_activity_threshold == 0.7
    assert request.kinase_activity_min_substrates == 5
    assert request.kinase_activity_top_n_substrates == 12
