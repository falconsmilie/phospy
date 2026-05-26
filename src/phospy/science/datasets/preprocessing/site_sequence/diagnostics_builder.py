"""Diagnostics builder collaborator for site-sequence resolution stage."""

from __future__ import annotations

from collections import defaultdict

from phospy.science.datasets.preprocessing.site_sequence.models import (
    SiteSequenceConflictResolution,
    SiteSequenceResolutionContext,
    SiteSequenceResolutionDiagnostics,
    SiteSequenceResolutionRowRequest,
)


class SiteSequenceDiagnosticsBuilder:
    """Collect and shape stage diagnostics while preserving prior schema."""

    def __init__(
        self,
        *,
        context: SiteSequenceResolutionContext,
        resolver_version: str,
    ) -> None:
        self._context = context
        self._resolver_version = resolver_version
        self._resolved_site_count = 0
        self._unresolved_site_count = 0
        self._unresolved_counts_by_reason: dict[str, int] = defaultdict(int)
        self._existing_sequence_conflict_count = 0
        self._filled_missing_count = 0
        self._replaced_existing_count = 0
        self._preserved_existing_count = 0
        self._row_status: list[dict[str, object]] = []
        self._row_diagnostics: list[dict[str, object]] = []
        self._has_conflict_errors = False

    def record_pre_resolution_skip(self, row: SiteSequenceResolutionRowRequest) -> None:
        if (
            row.skip_status is None
            or row.skip_reason is None
            or row.skip_action is None
        ):
            return
        if row.skip_status == "preserved_existing":
            self._preserved_existing_count += 1
            resolved_site_sequence = row.existing_site_sequence
        else:
            resolved_site_sequence = None
        self._append_row_status(
            row_key=row.row_key,
            status=row.skip_status,
            reason=row.skip_reason,
        )
        self._append_row_diagnostic(
            row=row,
            status=row.skip_status,
            existing_site_sequence=row.existing_site_sequence,
            fasta_site_sequence=None,
            resolved_site_sequence=resolved_site_sequence,
            action=row.skip_action,
            reason=row.skip_reason,
        )

    def record_resolved_missing_sequence(
        self,
        *,
        row: SiteSequenceResolutionRowRequest,
    ) -> None:
        self._resolved_site_count += 1
        self._unresolved_site_count += 1
        self._unresolved_counts_by_reason["resolved_missing_sequence"] += 1
        if row.has_existing:
            self._preserved_existing_count += 1
        status = "resolved_missing_sequence"
        reason = "resolver reported success without sequence output"
        self._append_row_status(
            row_key=row.row_key,
            status=status,
            reason=reason,
        )
        self._append_row_diagnostic(
            row=row,
            status=status,
            existing_site_sequence=row.existing_site_sequence,
            fasta_site_sequence=None,
            resolved_site_sequence=row.existing_site_sequence,
            action="preserve_existing" if row.has_existing else "unresolved",
            reason=reason,
        )

    def record_existing_resolution(
        self,
        *,
        row: SiteSequenceResolutionRowRequest,
        fasta_site_sequence: str,
        conflict: SiteSequenceConflictResolution,
    ) -> None:
        self._resolved_site_count += 1
        if conflict.is_conflict:
            self._existing_sequence_conflict_count += 1
            if conflict.action in {"preserve_existing", "error"}:
                self._preserved_existing_count += 1
            if conflict.action == "replace_existing":
                self._replaced_existing_count += 1
            if conflict.is_error:
                self._has_conflict_errors = True
        else:
            self._preserved_existing_count += 1
        self._append_row_status(
            row_key=row.row_key,
            status=conflict.status,
            reason=conflict.reason,
        )
        self._append_row_diagnostic(
            row=row,
            status=conflict.status,
            existing_site_sequence=row.existing_site_sequence,
            fasta_site_sequence=fasta_site_sequence,
            resolved_site_sequence=conflict.resolved_site_sequence,
            action=conflict.action,
            reason=conflict.reason,
        )

    def record_fill_missing(
        self,
        *,
        row: SiteSequenceResolutionRowRequest,
        resolved_site_sequence: str,
    ) -> None:
        self._resolved_site_count += 1
        self._filled_missing_count += 1
        status = "resolved"
        reason = "missing site_sequence resolved from FASTA"
        self._append_row_status(
            row_key=row.row_key,
            status=status,
            reason=reason,
        )
        self._append_row_diagnostic(
            row=row,
            status=status,
            existing_site_sequence=None,
            fasta_site_sequence=resolved_site_sequence,
            resolved_site_sequence=resolved_site_sequence,
            action="fill_missing",
            reason=reason,
        )

    def record_unresolved(
        self,
        *,
        row: SiteSequenceResolutionRowRequest,
        status: str,
        reason: str | None,
    ) -> None:
        self._unresolved_site_count += 1
        self._unresolved_counts_by_reason[status] += 1
        if row.has_existing:
            self._preserved_existing_count += 1
        self._append_row_status(
            row_key=row.row_key,
            status=status,
            reason=reason,
        )
        self._append_row_diagnostic(
            row=row,
            status=status,
            existing_site_sequence=row.existing_site_sequence,
            fasta_site_sequence=None,
            resolved_site_sequence=row.existing_site_sequence,
            action="preserve_existing" if row.has_existing else "unresolved",
            reason=reason,
        )

    def build(
        self,
        *,
        accession_column: str,
        site_column: str,
    ) -> SiteSequenceResolutionDiagnostics:
        payload = {
            "fasta_source_path": self._context.fasta_source_path,
            "fasta_source_label": self._context.fasta_source_label,
            "fasta_sha256": self._context.fasta_sha256,
            "resolver_version": self._resolver_version,
            "flank_size": int(self._context.flank_size),
            "mode": self._context.mode.value,
            "conflict_policy": self._context.conflict_policy.value,
            "accession_column": accession_column,
            "site_column": site_column,
            "resolved_site_count": int(self._resolved_site_count),
            "unresolved_site_count": int(self._unresolved_site_count),
            "unresolved_counts_by_reason": {
                key: int(self._unresolved_counts_by_reason[key])
                for key in sorted(self._unresolved_counts_by_reason)
            },
            "existing_sequence_conflict_count": int(
                self._existing_sequence_conflict_count
            ),
            "filled_missing_count": int(self._filled_missing_count),
            "replaced_existing_count": int(self._replaced_existing_count),
            "preserved_existing_count": int(self._preserved_existing_count),
            "row_status": self._row_status,
            "row_diagnostics": self._row_diagnostics,
        }
        return SiteSequenceResolutionDiagnostics(
            payload=payload,
            has_conflict_errors=self._has_conflict_errors,
        )

    def _append_row_status(
        self, *, row_key: str, status: str, reason: str | None
    ) -> None:
        self._row_status.append({"row_id": row_key, "status": status, "reason": reason})

    def _append_row_diagnostic(
        self,
        *,
        row: SiteSequenceResolutionRowRequest,
        status: str,
        existing_site_sequence: str | None,
        fasta_site_sequence: str | None,
        resolved_site_sequence: str | None,
        action: str,
        reason: str | None,
    ) -> None:
        self._row_diagnostics.append(
            {
                "row_index": int(row.row_index),
                "row_id": row.row_key,
                "site_id": row.site_id,
                "status": status,
                "existing_site_sequence": existing_site_sequence,
                "fasta_site_sequence": fasta_site_sequence,
                "resolved_site_sequence": resolved_site_sequence,
                "action": action,
                "reason": reason,
                "conflict_policy": self._context.conflict_policy.value,
                "resolver_version": self._resolver_version,
                "fasta_source_path": self._context.fasta_source_path,
                "fasta_sha256": self._context.fasta_sha256,
            }
        )


__all__ = ["SiteSequenceDiagnosticsBuilder"]
