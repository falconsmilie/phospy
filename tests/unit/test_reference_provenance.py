from __future__ import annotations

import pandas as pd

from phospy.provenance.models import (
    EnvironmentProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.serialization import from_payload, to_payload
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.references.resolution import ReferenceResolver
from phospy.references.resources import bundled_reference_name_for_organism


def test_explicit_reference_bundle_defaults_to_explicit_provenance() -> None:
    bundle = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_A"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )

    assert bundle.provenance is not None
    assert bundle.provenance.source_type == "explicit"
    assert bundle.provenance.bundle_id is None
    fingerprint_names = {item.name for item in bundle.provenance.table_fingerprints}
    assert fingerprint_names == {
        "references.kinase_substrate_map",
        "references.site_sequences",
    }


def test_bundled_reference_resolution_sets_bundled_provenance() -> None:
    resolved = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert resolved.provenance is not None
    assert resolved.provenance.source_type == "bundled"
    assert resolved.provenance.bundle_id == bundled_reference_name_for_organism(
        Organism.RAT
    )
    assert resolved.provenance.organism == Organism.RAT.value


def test_reference_resolver_keeps_explicit_bundle_identity_and_provenance() -> None:
    bundle = ReferenceBundle(
        organism=Organism.MOUSE,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["AKT1"], "substrate_site": ["AKT1;T308;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_B"]},
            index=pd.Index(["AKT1;T308;"], name="site_id"),
        ),
    )

    resolved = ReferenceResolver().run(
        bundle,
        dataset_organism=Organism.MOUSE,
    )
    assert resolved is bundle
    assert resolved.provenance is not None
    assert resolved.provenance.source_type == "explicit"


def test_explicit_reference_bundle_provenance_includes_identifier_normalisation() -> (
    None
):
    bundle = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": [" akt1 "], "substrate_site": [" mapk1 ; s123 "]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_A"]},
            index=pd.Index([" mapk1 ; s123 "], name="site_id"),
        ),
    )

    assert bundle.provenance is not None
    report = bundle.provenance.identifier_normalisation
    assert report is not None
    assert report.changed_identifier_count >= 2
    assert report.original_row_count == 2
    assert report.normalised_row_count == 2


def test_reference_provenance_serialization_round_trip_preserves_identifier_normalisation() -> (
    None
):
    bundle = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": [" akt1 "], "substrate_site": [" mapk1 ; s123 "]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_A"]},
            index=pd.Index([" mapk1 ; s123 "], name="site_id"),
        ),
    )
    assert bundle.provenance is not None

    run_provenance = RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="3.13",
            dependency_versions={},
            platform={},
        ),
        input_tables=(),
        preprocessing_stages=(),
        reference=bundle.provenance,
        workflow_name="unit_test",
        workflow_parameters={},
        random_state=None,
        random_seed_policy=None,
        output_tables=(),
    )

    payload = to_payload(run_provenance)
    restored = from_payload(payload)

    assert restored.reference is not None
    assert restored.reference.identifier_normalisation is not None
    assert (
        restored.reference.identifier_normalisation.changed_identifier_count
        == bundle.provenance.identifier_normalisation.changed_identifier_count
    )


def test_reference_provenance_from_payload_supports_legacy_missing_identifier_normalisation() -> (
    None
):
    table = TableFingerprint(
        name="references.kinase_substrate_map",
        rows=1,
        columns=2,
        index_name=None,
        column_names=("kinase", "substrate_site"),
        dtypes=("object", "object"),
        hash_algorithm="sha256",
        hash_value="a" * 64,
    )
    run_provenance = RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="3.13",
            dependency_versions={},
            platform={},
        ),
        input_tables=(),
        preprocessing_stages=(),
        reference=ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["AKT1"], "substrate_site": ["MAPK1;S123;"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["SEQ_A"]},
                index=pd.Index(["MAPK1;S123;"], name="site_id"),
            ),
        ).provenance,
        workflow_name="unit_test",
        workflow_parameters={},
        random_state=None,
        random_seed_policy=None,
        output_tables=(table,),
    )
    payload = to_payload(run_provenance)
    assert isinstance(payload["reference"], dict)
    payload_reference = dict(payload["reference"])
    payload_reference.pop("identifier_normalisation", None)
    payload["reference"] = payload_reference

    restored = from_payload(payload)
    assert restored.reference is not None
    assert restored.reference.identifier_normalisation is None
