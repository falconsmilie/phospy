from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy.api import (
    Organism,
    ReferenceBundleBuilder,
    ReferenceBundleBuildRequest,
)


def main() -> None:
    with TemporaryDirectory() as tmp_dir:
        root = Path(tmp_dir)
        kinase_path = root / "mouse_kinase_substrates.csv"
        sequence_path = root / "mouse_site_sequences.csv"
        pd.DataFrame(
            {
                "kinase": ["Map2k1", "Akt1"],
                "site_id": ["Mapk1;Y185;", "Gsk3b;S9;"],
                "organism": ["mouse", "mouse"],
            }
        ).to_csv(kinase_path, index=False)
        pd.DataFrame(
            {
                "site_id": ["Mapk1;Y185;", "Gsk3b;S9;"],
                "site_sequence": [_window("Y"), _window("S")],
                "gene_symbol": ["Mapk1", "Gsk3b"],
                "protein_id": ["MAPK1_MOUSE", "GSK3B_MOUSE"],
                "organism": ["mouse", "mouse"],
            }
        ).to_csv(sequence_path, index=False)

        references = ReferenceBundleBuilder().run(
            ReferenceBundleBuildRequest(
                organism=Organism.MOUSE,
                kinase_substrate_path=kinase_path,
                site_sequence_path=sequence_path,
                source_name="synthetic local mouse reference",
                source_version="demo-v1",
                retrieved_at="2026-06-11",
                license="synthetic demo data",
                redistribution_status="redistributable synthetic fixture",
                identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
            )
        )

    print("Reference bundle builder demo")
    print(f"Organism: {references.organism.value}")
    print(f"Kinase-substrate rows: {len(references.kinase_substrate_map)}")
    print(f"Sequence rows: {len(references.site_sequences)}")
    print(
        f"Source type: {references.provenance.source_type if references.provenance else None}"
    )


def _window(center: str) -> str:
    return f"{'A' * 15}{center}{'A' * 15}"


if __name__ == "__main__":
    main()
