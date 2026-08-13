from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.api import (
    EnrichmentDerivedQuantitativeSetProvenance,
    EnrichmentDerivedSetMissingValueRule,
    EnrichmentDerivedSetSourceResultKind,
    EnrichmentDerivedSetThresholdDirection,
    EnrichmentDerivedSetValueMeaning,
    EnrichmentDerivedSetValueScale,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    InputIntensityScaleEvidence,
)
from phospy.provenance import fingerprint_table


def _source_result_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["AKT1", "MAPK1"],
            "log2_fc": [1.25, 1.10],
            "adjusted_p_value": [0.01, 0.03],
        }
    )


def _derived_quantitative_provenance(
    source_result_table: pd.DataFrame,
) -> EnrichmentDerivedQuantitativeSetProvenance:
    return EnrichmentDerivedQuantitativeSetProvenance(
        source_result_fingerprint=fingerprint_table(
            source_result_table,
            name="differential.contrast_table[stim_vs_ctrl]",
        ),
        source_result_kind=EnrichmentDerivedSetSourceResultKind.CONTRAST,
        source_profile_or_contrast="stim_vs_ctrl",
        identifier_namespace="gene_symbol",
        threshold=1.0,
        direction=EnrichmentDerivedSetThresholdDirection.GREATER_THAN_OR_EQUAL,
        missing_value_rule=(
            EnrichmentDerivedSetMissingValueRule.TREAT_MISSING_AS_NOT_SELECTED
        ),
        quantitative_scale=EnrichmentDerivedSetValueScale.LOG2,
        quantitative_meaning=(
            EnrichmentDerivedSetValueMeaning.CONTRAST_LOG2_FOLD_CHANGE
        ),
        software_version="phospy-example",
    )


def _collection() -> GeneSetCollection:
    return GeneSetCollection(
        sets={
            "kinase_response": ("AKT1", "MAPK1", "MTOR"),
            "cell_cycle": ("CDK1", "CDK2", "MAPK1"),
        },
        identifier_kind="gene_symbol",
        term_names={
            "kinase_response": "Kinase response",
            "cell_cycle": "Cell cycle",
        },
        source_name="example in-memory gene sets",
    )


def _request(
    *,
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance,
    input_table: pd.DataFrame | None = None,
) -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind="gene_symbol",
        set_collection=_collection(),
        input_table=input_table,
        selected_identifiers=None if input_table is not None else ("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        selected_identifier_provenance=selected_identifier_provenance,
    )


def main() -> None:
    source_result_table = _source_result_table()
    manual = EnrichmentIdentifierSetProvenance(
        source_type=EnrichmentIdentifierSetSourceType.MANUAL,
        source_label="curated manual hit list",
        identifier_count=2,
    )
    phospy_derived = EnrichmentIdentifierSetProvenance(
        source_type=EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        source_label="differential significant genes",
        identifier_count=2,
        upstream_workflow_id="differential-workflow",
        upstream_result_id="contrast-a-vs-b",
        input_intensity_scale_evidence=InputIntensityScaleEvidence(
            input_intensity_scale="log2",
            input_intensity_scale_evidence_level="observed_transformation",
            input_intensity_scale_source=(
                "phospy.science.transformations.transformers.log2"
            ),
            input_intensity_scale_source_detail="log2_transform",
        ),
        derived_quantitative_provenance=_derived_quantitative_provenance(
            source_result_table
        ),
    )

    workflow = EnrichmentWorkflow()
    manual_result = workflow.run(_request(selected_identifier_provenance=manual))
    derived_result = workflow.run(
        _request(
            selected_identifier_provenance=phospy_derived,
            input_table=source_result_table,
        )
    )

    print("Enrichment provenance demo")
    print(f"Manual provenance: {manual_result.selected_identifier_provenance}")
    print(f"PhosPy-derived provenance: {derived_result.selected_identifier_provenance}")
    print(f"Statistics unchanged: {manual_result.records == derived_result.records}")
    assert derived_result.provenance is not None
    selected_provenance_payload = derived_result.provenance.workflow_parameters[
        "selected_identifier_provenance"
    ]
    if not isinstance(selected_provenance_payload, Mapping):
        raise RuntimeError("selected identifier provenance payload was not serialized")
    print(f"Run-provenance source type: {selected_provenance_payload['source_type']}")


if __name__ == "__main__":
    main()
