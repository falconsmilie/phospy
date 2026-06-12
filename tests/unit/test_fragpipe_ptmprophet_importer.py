from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api import PhosphositeImportResult
from phospy.errors import PhosPyInputError
from phospy.io.readers import (
    FragPipeColumnMapping,
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "fragpipe"


def _import_fixture(**kwargs: object) -> PhosphositeImportResult:
    return FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=FIXTURES / "ptmprophet_sites.tsv",
            **kwargs,
        )
    )


def test_fragpipe_ptmprophet_importer_reads_single_site_peptide() -> None:
    result = _import_fixture()

    assert isinstance(result, PhosphositeImportResult)
    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["Control", "Stim"]
    assert phospho.shape == (3, 2)
    assert float(phospho.iloc[0]["Control"]) == pytest.approx(10.0)
    metadata = result.site_metadata_candidate
    first = metadata.iloc[0]
    assert first["gene_symbol"] == "MAPK1"
    assert first["site"] == "S10"
    assert first["protein_accession"] == "P28482"
    assert first["localisation_confidence"] == pytest.approx(0.95)
    assert first["fragpipe_ptmprophet_candidate_sites"] == "S10"
    assert first["fragpipe_ptmprophet_site_probabilities"] == "S10:0.95"
    assert bool(first["fragpipe_ptmprophet_ambiguous"]) is False
    assert result.localisation_confidence_column == "localisation_confidence"
    assert result.diagnostics["fragpipe"]["filtering"]["removed_rows"] == 1


def test_fragpipe_ptmprophet_importer_retains_multi_site_peptide_evidence() -> None:
    result = _import_fixture()

    metadata = result.site_metadata_candidate
    assert metadata.iloc[1]["site"] == "S473,T475"
    assert metadata.iloc[1]["localisation_confidence"] == pytest.approx(0.97)
    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.iloc[1]["site_string"] == "S473;T475"
    assert evidence.iloc[1]["site_id"] == "AKT1;S473,T475;"
    assert bool(evidence.iloc[1]["multi_site"]) is True
    assert evidence.iloc[1]["modified_peptide_sequence"] == "AA[pS]A[pT]AA"
    assert evidence.iloc[1]["fragpipe_modified_peptide_phospho_count"] == 2
    assert any("multi-site candidates" in warning for warning in result.warnings)


def test_fragpipe_ptmprophet_ambiguous_localisation_is_joint_not_first_site() -> None:
    result = _import_fixture()

    metadata = result.site_metadata_candidate
    ambiguous = metadata.iloc[2]
    assert ambiguous["gene_symbol"] == "GSK3B"
    assert ambiguous["site"] == "S10,T11"
    assert ambiguous["localisation_confidence"] == pytest.approx(0.50)
    assert ambiguous["fragpipe_ptmprophet_candidate_sites"] == "S10;T11"
    assert ambiguous["fragpipe_ptmprophet_site_probabilities"] == ("S10:0.5;T11:0.5")
    assert bool(ambiguous["fragpipe_ptmprophet_ambiguous"]) is True
    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.iloc[2]["site_string"] == "S10;T11"
    assert evidence.iloc[2]["site_id"] == "GSK3B;S10,T11;"
    assert bool(evidence.iloc[2]["multi_site"]) is True
    assert any("ambiguous localisation" in warning for warning in result.warnings)
    assert (
        result.diagnostics["fragpipe"]["adaptation"]["ambiguous_localisation_rows"] == 1
    )


def test_fragpipe_ptmprophet_malformed_localisation_string_fails() -> None:
    source = pd.read_csv(FIXTURES / "ptmprophet_sites.tsv", sep="\t")
    source.loc[0, "PTMProphet Probability"] = "S4=not_a_probability"

    with pytest.raises(PhosPyInputError, match="malformed PTMProphet"):
        FragPipePTMProphetImporter().run(FragPipePTMProphetImportRequest(source=source))


def test_fragpipe_ptmprophet_contaminants_can_be_flagged() -> None:
    result = _import_fixture(contaminant_policy="flag")

    metadata = result.site_metadata_candidate
    assert metadata.shape[0] == 4
    assert metadata.iloc[3]["protein_accession"] == "P99999"
    assert metadata.loc[:, "fragpipe_contaminant"].tolist() == [
        False,
        False,
        False,
        True,
    ]
    evidence = result.peptide_evidence
    assert evidence is not None
    assert "fragpipe_contaminant" in evidence.columns
    assert bool(evidence.iloc[3]["fragpipe_contaminant"]) is True


def test_fragpipe_ptmprophet_supports_protein_position_localisation_strings() -> None:
    source = pd.DataFrame(
        {
            "protein": ["sp|Q16539|MK14_HUMAN"],
            "gene": ["MAPK14"],
            "peptide": ["AAAAAYAAAA"],
            "mod": ["AAAAA[pY]AAAA"],
            "ptm": ["Y182(0.93)"],
            "raw control": ["100.0"],
            "raw stim": ["120.0"],
        }
    )

    result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=source,
            column_mapping=FragPipeColumnMapping(
                protein_accession="protein",
                gene_symbol="gene",
                peptide_sequence="peptide",
                modified_peptide_sequence="mod",
                ptmprophet_probabilities="ptm",
                intensity_columns={
                    "raw control": "control",
                    "raw stim": "stim",
                },
            ),
            ptmprophet_position_reference="protein",
            source_name="custom_fragpipe",
        )
    )

    assert result.source_name == "custom_fragpipe"
    assert result.sample_column_mapping == {
        "raw control": "control",
        "raw stim": "stim",
    }
    assert result.site_metadata_candidate.iloc[0]["site"] == "Y182"
    assert result.site_metadata_candidate.iloc[0][
        "localisation_confidence"
    ] == pytest.approx(0.93)
