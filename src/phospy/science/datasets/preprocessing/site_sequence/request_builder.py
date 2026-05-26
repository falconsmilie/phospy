"""Row-level resolution request construction for site-sequence stage."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.science.datasets.preprocessing.policy_models import (
    SiteSequenceResolutionMode,
)
from phospy.science.datasets.preprocessing.site_sequence.models import (
    SiteSequenceResolutionContext,
    SiteSequenceResolutionRowRequest,
)
from phospy.science.sequences import PhosphositeSequenceResolutionRequest


class SiteSequenceResolutionRequestBuilder:
    """Build row-level resolver requests and mode-driven skip decisions."""

    def build(
        self,
        *,
        site_metadata: pd.DataFrame,
        normalized_existing_site_sequence: pd.Series,
        context: SiteSequenceResolutionContext,
    ) -> tuple[SiteSequenceResolutionRowRequest, ...]:
        self.require_columns(
            site_metadata=site_metadata,
            accession_column=context.accession_column,
            site_column=context.site_column,
        )
        requests: list[SiteSequenceResolutionRowRequest] = []
        for row_index, row_id in enumerate(site_metadata.index.tolist()):
            row_key = str(row_id)
            existing_value = normalized_existing_site_sequence.loc[row_id]
            has_existing = not bool(pd.isna(existing_value))
            existing_site_sequence = (
                None if pd.isna(existing_value) else str(existing_value)
            )
            skip_status: str | None = None
            skip_reason: str | None = None
            skip_action: str | None = None
            resolver_request: PhosphositeSequenceResolutionRequest | None = None

            if (
                has_existing
                and context.mode is SiteSequenceResolutionMode.FILL_MISSING_ONLY
            ):
                skip_status = "preserved_existing"
                skip_reason = "existing site_sequence preserved (fill_missing_only)"
                skip_action = "preserve_existing"
            elif (
                not has_existing
            ) and context.mode is SiteSequenceResolutionMode.VALIDATE_EXISTING_ONLY:
                skip_status = "missing_existing_sequence"
                skip_reason = "missing site_sequence preserved (validate_existing_only)"
                skip_action = "skip_missing_existing"
            else:
                resolver_request = PhosphositeSequenceResolutionRequest(
                    accession=site_metadata.at[row_id, context.accession_column],
                    site_token=site_metadata.at[row_id, context.site_column],
                    flank_size=context.flank_size,
                )
            requests.append(
                SiteSequenceResolutionRowRequest(
                    row_index=int(row_index),
                    row_id=row_id,
                    row_key=row_key,
                    site_id=row_key,
                    has_existing=has_existing,
                    existing_site_sequence=existing_site_sequence,
                    resolver_request=resolver_request,
                    skip_status=skip_status,
                    skip_reason=skip_reason,
                    skip_action=skip_action,
                )
            )
        return tuple(requests)

    @staticmethod
    def require_columns(
        *,
        site_metadata: pd.DataFrame,
        accession_column: str,
        site_column: str,
    ) -> None:
        missing_columns = [
            column
            for column in (accession_column, site_column)
            if column not in site_metadata.columns
        ]
        if not missing_columns:
            return
        joined = ", ".join(missing_columns)
        raise UnsupportedInputFormatError(
            "dataset preprocessing stage 'site_sequence_resolution' requires "
            f"site_metadata column(s): {joined}"
        )


__all__ = ["SiteSequenceResolutionRequestBuilder"]
