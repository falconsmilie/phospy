from __future__ import annotations

import pandas as pd

from phospy.io.readers import (
    MaxQuantColumnMapping,
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)


def main() -> None:
    source = pd.DataFrame(
        {
            "Leading proteins": ["P28482", "P31749", "CON__P99999"],
            "Gene names": ["MAPK1", "AKT1", "CONGENE"],
            "Amino acid": ["S", "S", "S"],
            "Positions within proteins": ["10", "473", "20"],
            "Localization prob": [0.95, 0.91, 0.99],
            "Sequence": ["AAAAASAAAA", "BBBBBSBBBB", "CCCCCSCCCC"],
            "Modified sequence": [
                "AAAAA(ph)SAAAA",
                "BBBBB(ph)SBBBB",
                "CCCCC(ph)SCCCC",
            ],
            "Intensity Control": [10.0, 20.0, 30.0],
            "Intensity Stim": [12.0, 21.0, 31.0],
            "Potential contaminant": ["", "", "+"],
            "Reverse": ["", "", ""],
        }
    )

    imported = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=source,
            column_mapping=MaxQuantColumnMapping(
                protein_accession="Leading proteins",
                gene_symbol="Gene names",
                amino_acid="Amino acid",
                site_position="Positions within proteins",
                localisation_confidence="Localization prob",
                peptide_sequence="Sequence",
                modified_peptide_sequence="Modified sequence",
                intensity_columns={
                    "Intensity Control": "control",
                    "Intensity Stim": "stim",
                },
                potential_contaminant="Potential contaminant",
                reverse="Reverse",
            ),
            contaminant_policy="remove",
            reverse_policy="remove",
        )
    )

    print("MaxQuant importer candidate matrix:")
    print(imported.phospho_matrix_candidate)
    print("MaxQuant filtering diagnostics:")
    print(imported.diagnostics["maxquant"]["filtering"])


if __name__ == "__main__":
    main()
