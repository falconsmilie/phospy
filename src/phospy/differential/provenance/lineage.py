"""Run-lineage composition for Differential-derived quantitative datasets."""

from phospy.datasets.derived_quantitative import DerivedQuantitativeDataProvenance
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.models import EnvironmentProvenance
from phospy.references.provenance.context import ReferenceContext
from phospy.run_lineage.models import RunProvenance


def build_derived_quantitative_run_provenance(
    *,
    lineage: DerivedQuantitativeDataProvenance,
    environment: EnvironmentProvenance | None = None,
    reference_context: ReferenceContext | None = None,
) -> RunProvenance:
    """Build truthful run provenance for a derived quantitative dataset."""

    resolved_environment = (
        collect_environment_provenance() if environment is None else environment
    )
    return RunProvenance(
        environment=resolved_environment,
        input_tables=lineage.parent_dataset_fingerprints,
        preprocessing_stages=(),
        reference=None,
        workflow_name=lineage.derivation_type,
        workflow_parameters={
            "derived_quantitative_data": lineage.to_payload(),
            "construction": {
                "method": lineage.implementation,
                "dataset_type": lineage.derived_dataset_type,
                "parent_dataset_type": lineage.parent_dataset_type,
                "model_constructor": (
                    "DerivedAnalysisReadyPhosphoDataset._from_owned_derived_tables"
                ),
                "source_dataset_provenance_reused": False,
                "source_dataset_preprocessing_report_reused": False,
            },
        },
        random_state=None,
        random_seed_policy=None,
        output_tables=lineage.derived_dataset_fingerprints,
        scientific_policies=(),
        reference_context=reference_context,
    )


__all__ = ["build_derived_quantitative_run_provenance"]
