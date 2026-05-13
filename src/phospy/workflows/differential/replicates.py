"""Technical-replicate policy resolution for differential workflows."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import ExperimentalDesign, SampleDesignRecord
from phospy.science.differential.policy_models import TechnicalReplicatePolicy


@dataclass(frozen=True, slots=True)
class _ResolvedReplicateGroup:
    condition: str
    biological_replicate_id: str
    output_sample_id: str
    input_sample_ids: tuple[str, ...]
    technical_replicate_ids: tuple[str, ...]
    batch: str | None
    block: str | None


@dataclass(frozen=True, slots=True)
class TechnicalReplicateResolution:
    """Resolved dataset/design pair after technical-replicate policy handling."""

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    workflow_provenance: Mapping[str, object] | None = None


class TechnicalReplicateResolver:
    """Apply explicit technical-replicate policy before design-matrix assembly."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> TechnicalReplicateResolution:
        if not isinstance(technical_replicate_policy, TechnicalReplicatePolicy):
            raise WorkflowValidationError(
                "differential workflow request technical_replicate_policy must be "
                "TechnicalReplicatePolicy"
            )

        repeated_groups = self._find_repeated_biological_groups(design)
        if not repeated_groups:
            return TechnicalReplicateResolution(
                dataset=dataset,
                design=design,
                workflow_provenance=None,
            )

        if technical_replicate_policy == TechnicalReplicatePolicy.REJECT:
            rendered_groups = ", ".join(
                f"{condition}:{biological_replicate_id}"
                for condition, biological_replicate_id in repeated_groups
            )
            raise WorkflowValidationError(
                "experimental design contains repeated biological_replicate_id values "
                "within condition groups (condition:biological_replicate_id="
                f"{rendered_groups}). Technical replicates require explicit aggregation; "
                "set technical_replicate_policy='mean' or "
                "technical_replicate_policy='median'."
            )

        groups = self._build_resolved_groups(design=design)
        aggregated_dataset = self._aggregate_dataset(
            dataset=dataset,
            groups=groups,
            technical_replicate_policy=technical_replicate_policy,
        )
        aggregated_design = ExperimentalDesign(
            samples=tuple(
                SampleDesignRecord(
                    sample_id=group.output_sample_id,
                    condition=group.condition,
                    biological_replicate_id=group.biological_replicate_id,
                    technical_replicate_id=None,
                    batch=group.batch,
                    block=group.block,
                )
                for group in groups
            )
        )
        provenance_groups: list[dict[str, object]] = []
        for group in groups:
            provenance_groups.append(
                {
                    "condition": group.condition,
                    "biological_replicate_id": group.biological_replicate_id,
                    "output_sample_id": group.output_sample_id,
                    "input_sample_ids": list(group.input_sample_ids),
                    "technical_replicate_ids": list(group.technical_replicate_ids),
                    "n_technical_replicates": len(group.input_sample_ids),
                }
            )
        workflow_provenance: dict[str, object] = {
            "technical_replicate_policy": technical_replicate_policy.value,
            "groups": provenance_groups,
        }
        return TechnicalReplicateResolution(
            dataset=aggregated_dataset,
            design=aggregated_design,
            workflow_provenance=workflow_provenance,
        )

    @staticmethod
    def _find_repeated_biological_groups(
        design: ExperimentalDesign,
    ) -> tuple[tuple[str, str], ...]:
        group_counts: Counter[tuple[str, str]] = Counter()
        for record in design.samples:
            if record.biological_replicate_id is None:
                continue
            group_counts[(record.condition, record.biological_replicate_id)] += 1
        repeated = sorted(key for key, count in group_counts.items() if count > 1)
        return tuple(repeated)

    @staticmethod
    def _build_resolved_groups(
        *, design: ExperimentalDesign
    ) -> tuple[_ResolvedReplicateGroup, ...]:
        missing_biological_id_samples = [
            record.sample_id
            for record in design.samples
            if record.biological_replicate_id is None
        ]
        if missing_biological_id_samples:
            raise WorkflowValidationError(
                "technical_replicate_policy aggregation requires "
                "biological_replicate_id on every design sample; missing for: "
                + ", ".join(missing_biological_id_samples)
            )

        grouped_records: dict[tuple[str, str], list[SampleDesignRecord]] = {}
        for record in design.samples:
            biological_replicate_id = record.biological_replicate_id
            if biological_replicate_id is None:
                raise WorkflowValidationError(
                    "technical_replicate_policy aggregation requires "
                    "biological_replicate_id on every design sample"
                )
            group_key = (record.condition, biological_replicate_id)
            if group_key not in grouped_records:
                grouped_records[group_key] = []
            grouped_records[group_key].append(record)

        biological_id_counts: Counter[str] = Counter(
            biological_replicate_id
            for _, biological_replicate_id in grouped_records.keys()
        )
        resolved_groups: list[_ResolvedReplicateGroup] = []
        for (condition, biological_replicate_id), records in grouped_records.items():
            output_sample_id = biological_replicate_id
            if biological_id_counts[biological_replicate_id] > 1:
                output_sample_id = f"{condition}__{biological_replicate_id}"
            technical_ids = tuple(
                technical_id
                for technical_id in dict.fromkeys(
                    record.technical_replicate_id for record in records
                )
                if technical_id is not None
            )
            batch = TechnicalReplicateResolver._require_consistent_optional_field(
                records=tuple(records),
                field_name="batch",
            )
            block = TechnicalReplicateResolver._require_consistent_optional_field(
                records=tuple(records),
                field_name="block",
            )
            resolved_groups.append(
                _ResolvedReplicateGroup(
                    condition=condition,
                    biological_replicate_id=biological_replicate_id,
                    output_sample_id=output_sample_id,
                    input_sample_ids=tuple(record.sample_id for record in records),
                    technical_replicate_ids=technical_ids,
                    batch=batch,
                    block=block,
                )
            )
        return tuple(resolved_groups)

    @staticmethod
    def _require_consistent_optional_field(
        *,
        records: tuple[SampleDesignRecord, ...],
        field_name: str,
    ) -> str | None:
        values = [getattr(record, field_name) for record in records]
        first = values[0]
        if any(value != first for value in values):
            sample_ids = ", ".join(record.sample_id for record in records)
            biological_replicate_id = records[0].biological_replicate_id
            raise WorkflowValidationError(
                "technical replicate aggregation requires consistent "
                f"{field_name!r} within each condition+biological_replicate_id group; "
                f"condition={records[0].condition!r}, "
                f"biological_replicate_id={biological_replicate_id!r}, "
                f"samples={sample_ids}"
            )
        return first

    @staticmethod
    def _aggregate_dataset(
        *,
        dataset: AnalysisReadyPhosphoDataset,
        groups: tuple[_ResolvedReplicateGroup, ...],
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> AnalysisReadyPhosphoDataset:
        phospho = dataset._borrow_phospho_frame()
        required_sample_ids = [
            sample_id for group in groups for sample_id in group.input_sample_ids
        ]
        missing = sorted(set(required_sample_ids) - set(phospho.columns.astype(str)))
        if missing:
            raise WorkflowValidationError(
                "technical replicate aggregation design samples are missing from "
                "dataset.phospho columns: " + ", ".join(missing)
            )

        aggregated_phospho = TechnicalReplicateResolver._aggregate_numeric_matrix(
            matrix=phospho,
            groups=groups,
            technical_replicate_policy=technical_replicate_policy,
        )
        total = dataset._borrow_total_frame()
        aggregated_total = None
        if total is not None:
            aggregated_total = TechnicalReplicateResolver._aggregate_numeric_matrix(
                matrix=total,
                groups=groups,
                technical_replicate_policy=technical_replicate_policy,
            )
        sample_metadata = dataset._borrow_sample_metadata_frame()
        aggregated_sample_metadata = None
        if sample_metadata is not None:
            aggregated_sample_metadata = TechnicalReplicateResolver._aggregate_metadata(
                metadata=sample_metadata,
                groups=groups,
            )
        comparisons = dataset._borrow_comparisons_frame()
        return AnalysisReadyPhosphoDataset._from_owned(
            phospho=aggregated_phospho,
            site_metadata=dataset._borrow_site_metadata_frame().copy(deep=True),
            intensity_scale_state=dataset.intensity_scale_state,
            processing_state=dataset.processing_state,
            sample_metadata=aggregated_sample_metadata,
            total=None if aggregated_total is None else aggregated_total,
            comparisons=None if comparisons is None else comparisons.copy(deep=True),
            organism=dataset.organism,
            preprocessing_report=dataset.preprocessing_report,
            provenance=dataset.provenance,
        )

    @staticmethod
    def _aggregate_numeric_matrix(
        *,
        matrix: pd.DataFrame,
        groups: tuple[_ResolvedReplicateGroup, ...],
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> pd.DataFrame:
        aggregated_columns: list[pd.Series] = []
        output_labels: list[str] = []
        for group in groups:
            source = matrix.loc[:, list(group.input_sample_ids)]
            if technical_replicate_policy == TechnicalReplicatePolicy.MEAN:
                collapsed = source.mean(axis=1)
            elif technical_replicate_policy == TechnicalReplicatePolicy.MEDIAN:
                collapsed = source.median(axis=1)
            else:
                raise WorkflowValidationError(
                    "technical_replicate_policy must be 'mean' or 'median' for "
                    "aggregation"
                )
            collapsed.name = group.output_sample_id
            aggregated_columns.append(collapsed.astype(float))
            output_labels.append(group.output_sample_id)
        aggregated = pd.concat(aggregated_columns, axis=1)
        aggregated.columns = pd.Index(output_labels)
        aggregated.index = matrix.index.copy()
        return aggregated

    @staticmethod
    def _aggregate_metadata(
        *,
        metadata: pd.DataFrame,
        groups: tuple[_ResolvedReplicateGroup, ...],
    ) -> pd.DataFrame:
        rows: list[pd.Series] = []
        for group in groups:
            rows.append(metadata.loc[group.input_sample_ids[0], :].copy(deep=True))
        aggregated = pd.DataFrame(
            rows,
            index=pd.Index(
                [group.output_sample_id for group in groups],
                name=metadata.index.name,
            ),
            columns=metadata.columns.copy(),
        )
        return aggregated


__all__ = [
    "TechnicalReplicateResolution",
    "TechnicalReplicateResolver",
]
