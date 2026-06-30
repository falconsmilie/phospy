"""FragPipe/PTMProphet importer orchestration."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from phospy.contracts.requests import PhosphositeImportRequest
from phospy.contracts.results import PhosphositeImportResult
from phospy.errors.input import PhosPyInputError
from phospy.io.readers.fragpipe.columns import _resolve_fragpipe_columns
from phospy.io.readers.fragpipe.constants import (
    _ADAPTED_GENE_SYMBOL_COLUMN,
    _ADAPTED_LOCALISATION_COLUMN,
    _ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN,
    _ADAPTED_PEPTIDE_SEQUENCE_COLUMN,
    _ADAPTED_PEPTIDE_SITE_STRING_COLUMN,
    _ADAPTED_PROTEIN_ACCESSION_COLUMN,
    _ADAPTED_PROTEIN_ID_COLUMN,
    _ADAPTED_ROW_ID_COLUMN,
    _ADAPTED_SITE_COLUMN,
    _ADAPTED_SITE_SEQUENCE_COLUMN,
    _ADAPTED_UNIQUE_FEATURE_ID_COLUMN,
)
from phospy.io.readers.fragpipe.filtering import _apply_flag_policies
from phospy.io.readers.fragpipe.models import FragPipePTMProphetImportRequest
from phospy.io.readers.fragpipe.normalization import _adapt_fragpipe_source
from phospy.io.readers.fragpipe.raw import (
    read_fragpipe_source,
    require_fragpipe_source_columns,
)
from phospy.io.readers.fragpipe.reporting import _augment_mapped_result
from phospy.io.readers.importers import MappedPhosphositeTableImporter
from phospy.validation.datasets.fragpipe import (
    validate_fragpipe_flag_policy,
    validate_ptmprophet_position_reference,
)


class FragPipePTMProphetImporter:
    """Import FragPipe/PTMProphet phosphosite output into PhosPy candidates."""

    def __init__(
        self,
        *,
        mapped_importer: MappedPhosphositeTableImporter | None = None,
    ) -> None:
        self._mapped_importer = mapped_importer or MappedPhosphositeTableImporter()

    def run(
        self,
        request: FragPipePTMProphetImportRequest,
    ) -> PhosphositeImportResult:
        if not isinstance(request, FragPipePTMProphetImportRequest):
            raise PhosPyInputError(
                "FragPipe importer input must be a FragPipePTMProphetImportRequest"
            )
        contaminant_policy = validate_fragpipe_flag_policy(
            request.contaminant_policy,
            field_name="fragpipe import request contaminant_policy",
        )
        decoy_policy = validate_fragpipe_flag_policy(
            request.decoy_policy,
            field_name="fragpipe import request decoy_policy",
        )
        position_reference = validate_ptmprophet_position_reference(
            request.ptmprophet_position_reference,
            field_name="fragpipe import request ptmprophet_position_reference",
        )
        source = read_fragpipe_source(request.source)
        require_fragpipe_source_columns(source)
        resolved = _resolve_fragpipe_columns(
            source,
            request.column_mapping,
            intensity_column_prefixes=request.intensity_column_prefixes,
        )
        filtered, flags, filter_diagnostics, filter_warnings = _apply_flag_policies(
            source,
            resolved=resolved,
            contaminant_policy=contaminant_policy,
            decoy_policy=decoy_policy,
        )
        if filtered.empty:
            raise PhosPyInputError(
                "FragPipe importer removed all rows after contaminant/decoy filtering"
            )

        adapted, adapter_diagnostics, adapter_warnings = _adapt_fragpipe_source(
            filtered,
            resolved=resolved,
            ptmprophet_position_reference=position_reference,
        )
        mapped_result = self._mapped_importer.run(
            PhosphositeImportRequest(
                source=adapted,
                sample_intensity_columns=resolved.intensity_columns,
                gene_symbol_column=_ADAPTED_GENE_SYMBOL_COLUMN,
                site_column=_ADAPTED_SITE_COLUMN,
                row_id_column=_ADAPTED_ROW_ID_COLUMN,
                protein_id_column=_ADAPTED_PROTEIN_ID_COLUMN,
                protein_accession_column=_ADAPTED_PROTEIN_ACCESSION_COLUMN,
                site_sequence_column=(
                    _ADAPTED_SITE_SEQUENCE_COLUMN
                    if _ADAPTED_SITE_SEQUENCE_COLUMN in adapted.columns
                    else None
                ),
                localisation_confidence_column=_ADAPTED_LOCALISATION_COLUMN,
                localisation_confidence_scale="probability",
                unique_feature_id_column=_ADAPTED_UNIQUE_FEATURE_ID_COLUMN,
                peptide_sequence_column=_ADAPTED_PEPTIDE_SEQUENCE_COLUMN,
                modified_peptide_sequence_column=(
                    _ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN
                ),
                peptide_site_string_column=_ADAPTED_PEPTIDE_SITE_STRING_COLUMN,
                source_name=request.source_name,
            )
        )
        return _augment_mapped_result(
            mapped_result,
            adapted=adapted,
            flags=flags,
            contaminant_policy=contaminant_policy,
            decoy_policy=decoy_policy,
            resolved=resolved,
            filter_diagnostics=filter_diagnostics,
            adapter_diagnostics=adapter_diagnostics,
            warnings=filter_warnings + adapter_warnings,
        )


__all__ = ["FragPipePTMProphetImporter"]
