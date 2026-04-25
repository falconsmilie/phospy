from __future__ import annotations

import pandas as pd

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
