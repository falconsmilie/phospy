# Phosphosite Importers

Importers translate upstream search-engine phosphosite tables into PhosPy
dataset-builder input candidates. They do not construct
`AnalysisReadyPhosphoDataset` objects and they do not infer differential design.

Use importers when a table comes from an upstream search or quantification tool
and needs consistent column mapping before the strict dataset builder runs.

## Current Support Summary

Current public phosphosite importer support is:

- `MappedPhosphositeTableImporter` for caller-supplied explicit column mapping
- `MaxQuantPhosphositeImporter` for MaxQuant-style `Phospho (STY)Sites.txt`
  exports
- `FragPipePTMProphetImporter` for FragPipe/Philosopher peptide or site tables
  with PTMProphet localisation probabilities

Spectronaut and DIA-NN phosphosite importers are not currently implemented as
dedicated semantic importers. A caller may still use
`MappedPhosphositeTableImporter` for a manually mapped compatible table, but
that is generic mapped-table support, not Spectronaut/DIA-NN support and not
upstream statistical result import.

## Responsibility Boundary

Importers own:

- explicit source-column mapping
- sample intensity column extraction
- localisation-confidence normalisation to `localisation_confidence`
- site metadata candidate assembly
- optional peptide-evidence candidate assembly
- duplicate and multi-site diagnostics

`AnalysisReadyDatasetBuilder` still owns:

- input-source and analysis-ready validation
- site-key derivation
- site-sequence handling
- missing-value and localisation policy enforcement
- duplicate-site handling
- peptide-evidence resolution into site-level matrices

Importers never infer sample groups, contrasts, batches, or differential
designs from column names.

## MaxQuant Phosphosite Importer

Use `MaxQuantPhosphositeImporter` for MaxQuant-style
`Phospho (STY)Sites.txt` exports. The importer detects common MaxQuant column
names, but every scientific field can be overridden with
`MaxQuantColumnMapping`; do not rely on identical export headers across
MaxQuant versions or processing templates.

```python
from phospy.api import Organism
from phospy.io.readers import (
    MaxQuantColumnMapping,
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)

import_result = MaxQuantPhosphositeImporter().run(
    MaxQuantPhosphositeImportRequest(
        source="Phospho (STY)Sites.txt",
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
        localisation_confidence_scale="probability",
    )
)
```

The localisation output column is `localisation_confidence`. It is numeric and
normalised to a probability in `[0.0, 1.0]` by the shared importer helpers.
MaxQuant probability strings such as `S(0.95)` are parsed by the MaxQuant
adapter; raw score-difference columns are not converted into probabilities.
If you map a score-like column, it must already be threshold-ready on the
configured scale.

Contaminant and reverse handling is explicit:

- `contaminant_policy="remove"` and `reverse_policy="remove"` drop flagged rows.
- `contaminant_policy="flag"` and `reverse_policy="flag"` retain rows and add
  `maxquant_potential_contaminant` / `maxquant_reverse` metadata columns.
- `contaminant_policy="error"` and `reverse_policy="error"` fail if flagged rows
  are present.

The importer generates protein/site-scoped source row IDs by default and stores
protein context in `protein_accession` and `protein_id` where possible. It does
not use display IDs such as `GENE;S123;` as row identity.

Multi-site MaxQuant rows are retained as peptide evidence. Hand them to the
builder peptide-evidence lane when site-level resolution is required:

```python
request = import_result.to_dataset_build_request(
    site_resolution_mode="peptide_evidence",
    multi_site_policy="split",
    organism=Organism.HUMAN,
    input_intensity_scale="linear",
)
```

## FragPipe/PTMProphet Importer

Use `FragPipePTMProphetImporter` for FragPipe/Philosopher peptide or site
tables that include PTMProphet localisation probabilities. The importer parses
protein accessions, modified peptide strings, PTMProphet site probabilities,
protein positions, intensity columns, and contaminant/decoy flags before
delegating to the common PhosPy importer result contract.

```python
from phospy.io.readers import (
    FragPipeColumnMapping,
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
)

import_result = FragPipePTMProphetImporter().run(
    FragPipePTMProphetImportRequest(
        source="combined_peptide.tsv",
        column_mapping=FragPipeColumnMapping(
            protein_accession="Protein",
            gene_symbol="Gene",
            peptide_sequence="Peptide",
            modified_peptide_sequence="Modified Peptide",
            protein_start="Protein Start",
            ptmprophet_probabilities="PTMProphet Probability",
            intensity_columns={
                "Intensity Control": "control",
                "Intensity Stim": "stim",
            },
            contaminant="Contaminant",
            decoy="Decoy",
        ),
        contaminant_policy="remove",
        decoy_policy="remove",
    )
)
```

By default, PTMProphet positions such as `S4(0.95)` are interpreted as
peptide-relative positions and translated with `protein_start`. If your table
already reports protein-relative tokens such as `S473(0.98)`, set
`ptmprophet_position_reference="protein"`.

The importer does not accept localisation strings as opaque labels. Supported
probability tokens include forms such as `S4(0.95)` and `S4:0.95`; malformed
strings fail fast. Modified peptides are parsed for common phospho annotations
such as `[pS]`, `S[+79.9663]`, and `(ph)S`.

Ambiguous PTMProphet evidence is represented explicitly. If a single-phospho
peptide has tied top candidates, the importer emits a joint multi-site row
such as `S10,T11` rather than silently selecting the first site. It also adds
`fragpipe_ptmprophet_candidate_sites`,
`fragpipe_ptmprophet_site_probabilities`, and
`fragpipe_ptmprophet_ambiguous` metadata columns for auditability.

The localisation output column remains the shared
`localisation_confidence` probability column, so it can be enforced by
`DatasetLocalisationConfig` during dataset building.

## Generic Column-Mapped Importer

The foundation importer is intentionally generic. The current dedicated
tool-specific phosphosite importers are MaxQuant and FragPipe/PTMProphet. Any
future dedicated Spectronaut or DIA-NN importer should stay a small adapter
that maps known columns into this contract, emits `PhosphositeImportResult`
candidates, and still relies on dataset-builder validation.

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import Organism, PhosphositeImportRequest
from phospy.io.readers import MappedPhosphositeTableImporter

source = pd.DataFrame(
    {
        "gene": ["MAPK14", "GSK3B"],
        "site": ["Y182", "S9"],
        "protein": ["MAPK14", "GSK3B"],
        "sequence_window": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        ],
        "localisation_percent": [95.0, 92.0],
        "intensity_a": [10.0, 20.0],
        "intensity_b": [11.0, 21.0],
    }
)

import_result = MappedPhosphositeTableImporter().run(
    PhosphositeImportRequest(
        source=source,
        sample_intensity_columns={
            "intensity_a": "sample_a",
            "intensity_b": "sample_b",
        },
        gene_symbol_column="gene",
        site_column="site",
        protein_id_column="protein",
        site_sequence_column="sequence_window",
        localisation_confidence_column="localisation_percent",
        localisation_confidence_scale="percent",
    )
)

dataset = AnalysisReadyDatasetBuilder().run(
    import_result.to_dataset_build_request(
        organism=Organism.RAT,
        input_intensity_scale="linear",
    )
)
```

`import_result.phospho_matrix_candidate` and
`import_result.site_metadata_candidate` are defensive snapshots. Inspect
`import_result.warnings` and `import_result.diagnostics` before running the
builder when upstream output may contain duplicate site rows, multi-site
peptides, missing localisation scores, or invalid localisation scores.

## Peptide Evidence

When peptide-level evidence columns are mapped, importer output can be handed to
the builder peptide-evidence lane:

```python
request = import_result.to_dataset_build_request(
    site_resolution_mode="peptide_evidence",
    multi_site_policy="split",
    organism=Organism.RAT,
    input_intensity_scale="linear",
)

dataset = AnalysisReadyDatasetBuilder().run(request)
```

The explicit `multi_site_policy` is required. Ambiguous localisation or
multi-site rows are retained in importer output and reported; they are not
silently dropped by importers.
