from __future__ import annotations

from phospy.api import (
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)
from phospy.provenance.models import InputIntensityScaleEvidence


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
) -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind="gene_symbol",
        set_collection=_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        selected_identifier_provenance=selected_identifier_provenance,
    )


def main() -> None:
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
    )

    workflow = EnrichmentWorkflow()
    manual_result = workflow.run(_request(selected_identifier_provenance=manual))
    derived_result = workflow.run(
        _request(selected_identifier_provenance=phospy_derived)
    )

    print("Enrichment provenance demo")
    print(f"Manual provenance: {manual_result.selected_identifier_provenance}")
    print(f"PhosPy-derived provenance: {derived_result.selected_identifier_provenance}")
    print(f"Statistics unchanged: {manual_result.records == derived_result.records}")
    assert derived_result.provenance is not None
    print(
        "Run-provenance source type: "
        f"{derived_result.provenance.workflow_parameters['selected_identifier_provenance']['source_type']}"
    )


if __name__ == "__main__":
    main()
