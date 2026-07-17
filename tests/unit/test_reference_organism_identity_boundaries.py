from __future__ import annotations

from typing import cast

import pandas as pd
import pytest

from phospy.api import AnalysisReadyDatasetBuilder
from phospy.contracts.requests import DatasetBuildRequest
from phospy.errors import DatasetValidationError, UnsupportedOrganismError
from phospy.errors.references import ReferenceCompatibilityError
from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
    DerivedSampleMapping,
)
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import ReferenceProvenance, RunProvenance
from phospy.provenance.reference_context import ReferenceContext
from phospy.provenance.serialization import from_payload, to_payload
from phospy.science.datasets.builders.provenance_assembler import (
    DatasetRunProvenanceAssembler,
)
from phospy.science.datasets.derived_quantitative import (
    DerivedAnalysisReadyPhosphoDataset,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.references.models import (
    Organism,
    ReferenceBundle,
    ReferencePreset,
)
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.sites.site_keys import ProteinScopedPhosphositeKey
from phospy.science.transformations.models import IntensityScaleState
from phospy.validation.references.compatibility import ReferenceCompatibilityValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def test_reference_context_aliases_store_canonical_organism_enum_once() -> None:
    context = _context(organism="Rattus norvegicus")

    assert context.organism is Organism.RAT
    assert context.to_payload()["organism"] == "rat"

    restored = ReferenceContext.from_payload(
        {
            "organism": "10116",
            "protein_namespace": "gene_symbol",
            "source_name": "unit-reference",
            "source_version": "v1",
            "proteome_version": None,
            "reference_table_sha256": "a" * 64,
        }
    )
    assert restored.organism is Organism.RAT


def test_reference_context_rejects_unsupported_reference_organism() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="reference_context\\.organism has unsupported organism 'yeast'",
    ):
        _context(organism=cast(Organism, "yeast"))


def test_malformed_direct_site_key_construction_rejects_free_organism_string() -> None:
    with pytest.raises(ValueError, match="site_key\\.organism has unsupported"):
        ProteinScopedPhosphositeKey(
            organism=cast(Organism, "yeast"),
            protein_namespace="protein_id",
            protein_identifier="P1",
            residue="S",
            position=1,
        )


def test_direct_dataset_constructor_rejects_provenance_reference_context_mismatch() -> (
    None
):
    payload = _analysis_ready_payload(organism=Organism.RAT)
    payload["provenance"] = _run_provenance(reference_context=_context("human"))

    with pytest.raises(DatasetValidationError) as exc_info:
        _dataset_from_payload(payload)

    message = str(exc_info.value)
    assert "dataset organism identity conflict" in message
    assert "dataset.organism='rat'" in message
    assert "dataset.provenance.reference_context.organism='human'" in message


def test_trusted_dataset_constructor_rejects_provenance_reference_context_mismatch() -> (
    None
):
    payload = _analysis_ready_payload(organism=Organism.RAT)

    with pytest.raises(DatasetValidationError, match="dataset organism identity"):
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=_payload_phospho(payload),
            site_metadata=_payload_site_metadata(payload),
            organism=_payload_organism(payload),
            intensity_scale_state=_payload_intensity_scale_state(payload),
            processing_state=_payload_processing_state(payload),
            provenance=_run_provenance(reference_context=_context("human")),
        )


def test_derived_dataset_constructor_rejects_provenance_reference_context_mismatch() -> (
    None
):
    payload = _analysis_ready_payload(organism=Organism.RAT)
    phospho = payload["phospho"]
    assert isinstance(phospho, pd.DataFrame)

    with pytest.raises(DatasetValidationError, match="dataset organism identity"):
        DerivedAnalysisReadyPhosphoDataset.from_owned_derived_tables(
            phospho=phospho,
            site_metadata=_payload_site_metadata(payload),
            organism=_payload_organism(payload),
            intensity_scale_state=_payload_intensity_scale_state(payload),
            processing_state=_payload_processing_state(payload),
            derived_lineage=_derived_lineage(phospho),
            provenance=_run_provenance(reference_context=_context("human")),
        )


def test_builder_dataset_construction_rejects_bad_assembled_reference_context(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _bad_run(
        self: DatasetRunProvenanceAssembler, **_kwargs: object
    ) -> RunProvenance:
        return _run_provenance(reference_context=_context("human"))

    monkeypatch.setattr(DatasetRunProvenanceAssembler, "run", _bad_run)

    with pytest.raises(DatasetValidationError, match="dataset organism identity"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=pd.DataFrame(
                    {"sample_a": [1.0], "sample_b": [2.0]},
                    index=pd.Index(["MAPK14;Y182;"], name="site_id"),
                ),
                site_metadata=pd.DataFrame(
                    {
                        "gene_symbol": ["MAPK14"],
                        "site": ["Y182"],
                        "site_sequence": ["AAAAAAAYAAAAAAA"],
                        "protein_id": ["MAPK14"],
                        "localisation_confidence": [0.95],
                    },
                    index=pd.Index(["MAPK14;Y182;"], name="site_id"),
                ),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_deserialized_run_provenance_rejects_reference_context_organism_mismatch() -> (
    None
):
    provenance = _run_provenance(
        reference_context=_context("rat"),
        reference=ReferenceProvenance(
            source_type="explicit",
            organism=Organism.RAT,
            bundle_id=None,
            table_fingerprints=(),
        ),
    )
    payload = to_payload(provenance)
    reference_context = payload["reference_context"]
    assert isinstance(reference_context, dict)
    reference_context["organism"] = "human"

    with pytest.raises(
        ReferenceValidationError,
        match="Run provenance reference-context organism mismatch",
    ):
        from_payload(payload)


def test_reference_bundle_rejects_provenance_organism_mismatch() -> None:
    with pytest.raises(ReferenceValidationError) as exc_info:
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=_kinase_substrate_map(),
            site_sequences=_site_sequences(),
            provenance=ReferenceProvenance(
                source_type="explicit",
                organism=Organism.HUMAN,
                bundle_id=None,
                table_fingerprints=(),
            ),
        )

    message = str(exc_info.value)
    assert "reference bundle organism identity conflict" in message
    assert "references.organism='rat'" in message
    assert "references.provenance.organism='human'" in message


def test_reference_provenance_rejects_embedded_context_organism_mismatch() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="Reference provenance organism mismatch",
    ):
        ReferenceProvenance(
            source_type="explicit",
            organism=Organism.RAT,
            bundle_id=None,
            table_fingerprints=(),
            reference_context=_context("human"),
        )


def test_reference_auto_and_explicit_preset_mismatch_remain_defense_in_depth() -> None:
    assert (
        ReferenceCompatibilityValidator().resolve_preset_organism(
            preset=ReferencePreset.AUTO,
            dataset_organism="Rattus norvegicus",
        )
        is Organism.RAT
    )

    with pytest.raises(ReferenceCompatibilityError) as exc_info:
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.RAT,
        )
    message = str(exc_info.value)
    assert "dataset.organism='rat'" in message
    assert "references='human' resolved_to='human'" in message

    with pytest.raises(UnsupportedOrganismError, match="supported bundled organisms"):
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=Organism.HUMAN,
        )


def _context(organism: object = "rat") -> ReferenceContext:
    return ReferenceContext(
        organism=organism,
        protein_namespace="gene_symbol",
        source_name="unit-reference",
        source_version="v1",
        proteome_version=None,
        reference_table_sha256="a" * 64,
    )


def _run_provenance(
    *,
    reference_context: ReferenceContext | None,
    reference: ReferenceProvenance | None = None,
) -> RunProvenance:
    return RunProvenance(
        environment=collect_environment_provenance(),
        input_tables=(),
        preprocessing_stages=(),
        reference=reference,
        workflow_name="unit-test",
        workflow_parameters={},
        random_state=None,
        random_seed_policy=None,
        output_tables=(),
        reference_context=reference_context,
    )


def _analysis_ready_payload(*, organism: Organism) -> dict[str, object]:
    site_key = (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        "protein_identifier=MAPK14|residue=Y|position=182"
    )
    index = pd.Index([site_key], name="site_key")
    return {
        "phospho": pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [2.0]},
            index=index.copy(),
        ),
        "site_metadata": pd.DataFrame(
            {
                "site_key": [site_key],
                "display_id": ["MAPK14;Y182;"],
                "organism": ["rat"],
                "protein_namespace": ["protein_id"],
                "protein_identifier": ["MAPK14"],
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["AAAAAAAYAAAAAAA"],
                "protein_id": ["MAPK14"],
            },
            index=index.copy(),
        ),
        "organism": organism,
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        "processing_state": supported_linear_processing_state(has_total_matrix=False),
    }


def _dataset_from_payload(
    payload: dict[str, object],
) -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyPhosphoDataset(
        phospho=_payload_phospho(payload),
        site_metadata=_payload_site_metadata(payload),
        organism=_payload_organism(payload),
        intensity_scale_state=_payload_intensity_scale_state(payload),
        processing_state=_payload_processing_state(payload),
        provenance=cast(RunProvenance | None, payload.get("provenance")),
    )


def _payload_phospho(payload: dict[str, object]) -> pd.DataFrame:
    return cast(pd.DataFrame, payload["phospho"])


def _payload_site_metadata(payload: dict[str, object]) -> pd.DataFrame:
    return cast(pd.DataFrame, payload["site_metadata"])


def _payload_organism(payload: dict[str, object]) -> Organism:
    return cast(Organism, payload["organism"])


def _payload_intensity_scale_state(
    payload: dict[str, object],
) -> IntensityScaleState:
    return cast(IntensityScaleState, payload["intensity_scale_state"])


def _payload_processing_state(payload: dict[str, object]) -> DatasetProcessingState:
    return cast(DatasetProcessingState, payload["processing_state"])


def _derived_lineage(phospho: pd.DataFrame) -> DerivedQuantitativeDataProvenance:
    fingerprints = (fingerprint_table(phospho, name="dataset.phospho"),)
    return DerivedQuantitativeDataProvenance(
        derivation_type="unit-derived",
        parent_dataset_type="AnalysisReadyPhosphoDataset",
        derived_dataset_type="DerivedAnalysisReadyPhosphoDataset",
        parent_dataset_fingerprints=fingerprints,
        derived_dataset_fingerprints=fingerprints,
        sample_mapping=(
            DerivedSampleMapping(
                output_sample_id="sample_a",
                input_sample_ids=("sample_a",),
                condition="unit",
                biological_replicate_id="bio_a",
            ),
            DerivedSampleMapping(
                output_sample_id="sample_b",
                input_sample_ids=("sample_b",),
                condition="unit",
                biological_replicate_id="bio_b",
            ),
        ),
        aggregation_method="identity",
        input_intensity_scale="linear",
        output_intensity_scale="linear",
        quantitative_meaning="phosphosite intensity",
        missingness_policy={},
        matrices_transformed={"phospho": False},
        implementation="unit-test",
        implementation_version="test",
    )


def _kinase_substrate_map() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinase": ["AKT1"],
            "substrate_site": ["MAPK14;Y182;"],
        }
    )


def _site_sequences() -> pd.DataFrame:
    return pd.DataFrame(
        {"site_sequence": ["AAAAAAAYAAAAAAA"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
