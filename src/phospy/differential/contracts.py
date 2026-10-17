"""Public request contract for differential analysis."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.differential.config import DifferentialAnalysisConfig
from phospy.differential.design.models import Contrast, ExperimentalDesign


@dataclass(frozen=True, slots=True, eq=False)
class DifferentialAnalysisRequest:
    """Request for the public differential workflow.

    Construction stores the payload only. Dataset eligibility, design/sample
    alignment, contrast validity, replicate requirements, and config coherence
    are validated by ``DifferentialAnalysisWorkflow.run(...)`` before
    interpretation or statistical execution.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    config: DifferentialAnalysisConfig = field(
        default_factory=DifferentialAnalysisConfig
    )


__all__ = ["DifferentialAnalysisRequest"]
