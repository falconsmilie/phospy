"""Signalome bundle payload parsing and normalization helpers."""

from phospy.io.bundles._signalome.compatibility.config import (
    signalome_config_from_payload,
)
from phospy.io.bundles._signalome.compatibility.diagnostics import (
    signalome_alignment_diagnostics_from_payload,
    signalome_alignment_diagnostics_to_payload,
    signalome_module_selection_diagnostics_from_payload,
    signalome_module_selection_diagnostics_to_payload,
    signalome_network_correlation_diagnostics_from_payload,
    signalome_network_correlation_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_from_payload,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.io.bundles._signalome.compatibility.tables import (
    normalize_module_assignments_table,
)

__all__ = [
    "normalize_module_assignments_table",
    "signalome_alignment_diagnostics_from_payload",
    "signalome_alignment_diagnostics_to_payload",
    "signalome_config_from_payload",
    "signalome_module_selection_diagnostics_from_payload",
    "signalome_module_selection_diagnostics_to_payload",
    "signalome_network_correlation_diagnostics_from_payload",
    "signalome_network_correlation_diagnostics_to_payload",
    "signalome_score_preconditioning_diagnostics_from_payload",
    "signalome_score_preconditioning_diagnostics_to_payload",
]
