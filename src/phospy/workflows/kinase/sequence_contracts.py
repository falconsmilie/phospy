"""Kinase workflow sequence-context contract selection."""

from __future__ import annotations

from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.kinase_library_models import KinaseLibraryResource
from phospy.validation.identity_contracts import SequenceContextContract
from phospy.workflows.kinase.scoring_mode_contracts import (
    kinase_scoring_mode_input_contract,
)


def kinase_sequence_context_contract(
    *,
    scoring_mode: str,
    kinase_library_resource: object,
) -> SequenceContextContract | None:
    """Return the workflow sequence contract selected by kinase scoring config."""

    mode_contract = kinase_scoring_mode_input_contract(scoring_mode)
    if not mode_contract.requires_kinase_library_resource:
        return None
    if not isinstance(kinase_library_resource, KinaseLibraryResource):
        return None
    sequence_window = kinase_library_resource.sequence_window
    upstream = int(sequence_window.upstream_residues)
    downstream = int(sequence_window.downstream_residues)
    window_length = upstream + 1 + downstream
    return SequenceContextContract(
        requires_site_sequence=True,
        requires_centered_site=True,
        required_window_length=window_length,
        center_index=upstream,
        allowed_residues=frozenset("ACDEFGHIKLMNPQRSTVWY"),
        allow_terminal_padding=False,
        allow_lowercase=False,
        allow_modified_residue_symbols=False,
        required_center_residues=frozenset({"S", "T", "Y"}),
        requires_known_sequence_source=True,
        contract_id=f"kinase_library_fixed_{window_length}aa_centered_window",
    )


def dataset_sequence_source_label(dataset: AnalysisReadyPhosphoDataset) -> str | None:
    """Return a known dataset sequence source label when durable provenance exists."""

    site_sequence_resolution = dataset.processing_state.site_sequence_resolution
    if bool(site_sequence_resolution.configured):
        return "dataset_site_sequence_resolution"
    preprocessing_report = dataset.preprocessing_report
    if preprocessing_report is not None:
        report = preprocessing_report.site_sequence_resolution_summary()
        if report is not None and int(report.final_sequence_complete_sites) > 0:
            return "dataset_preprocessing_report"
    provenance = dataset.provenance
    if provenance is not None:
        derivation = provenance.workflow_parameters.get("site_sequence_derivation")
        if isinstance(derivation, dict):
            return "dataset_builder_site_sequence_derivation"
    return None


__all__ = [
    "dataset_sequence_source_label",
    "kinase_sequence_context_contract",
]
