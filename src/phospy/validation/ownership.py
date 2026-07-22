"""Validation ownership registry.

Each major validation rule has one primary owner. Higher-level layers may
compose owners, but should not duplicate the same rule checks.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ValidationRuleOwner:
    """One documented ownership entry for a major validation rule."""

    rule: str
    owner: str


VALIDATION_RULE_OWNERS: tuple[ValidationRuleOwner, ...] = (
    ValidationRuleOwner(
        rule="dataset build request input source types",
        owner="DatasetInputSourceValidator.run",
    ),
    ValidationRuleOwner(
        rule="dataset build preprocessing config policy",
        owner="DatasetPreprocessingConfigValidator.run",
    ),
    ValidationRuleOwner(
        rule="kinase workflow request config policy",
        owner="KinaseWorkflowConfigValidator.run",
    ),
    ValidationRuleOwner(
        rule="signalome workflow request config policy",
        owner="SignalomeConfigValidator.run",
    ),
    ValidationRuleOwner(
        rule="reference input compatibility (preset/bundle vs dataset organism)",
        owner="ReferenceCompatibilityValidator.run",
    ),
    ValidationRuleOwner(
        rule="reference bundle structural contract",
        owner="ReferenceBundleValidator.run",
    ),
    ValidationRuleOwner(
        rule="kinase library resource structural contract",
        owner="KinaseLibraryResourceValidator.run",
    ),
    ValidationRuleOwner(
        rule="analysis-ready dataset structural contract",
        owner="AnalysisReadyPhosphoDataset._init_analysis_ready_tables",
    ),
    ValidationRuleOwner(
        rule="dataset/intensity-scale-state coherence",
        owner="IntensityScaleStateValidator.run",
    ),
    ValidationRuleOwner(
        rule="sequence-aware workflow centred site-sequence context",
        owner="enforce_centred_site_sequence_context",
    ),
    ValidationRuleOwner(
        rule="signalome result expanded_signalome field type/ownership",
        owner="SignalomeWorkflowResult.__post_init__",
    ),
)


__all__ = ["VALIDATION_RULE_OWNERS", "ValidationRuleOwner"]
