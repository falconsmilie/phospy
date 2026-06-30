"""FragPipe/PTMProphet importer constants."""

from __future__ import annotations

import re

_ADAPTED_ROW_ID_COLUMN = "__phospy_fragpipe_row_id"
_ADAPTED_PROTEIN_ACCESSION_COLUMN = "__phospy_fragpipe_protein_accession"
_ADAPTED_PROTEIN_ID_COLUMN = "__phospy_fragpipe_protein_id"
_ADAPTED_GENE_SYMBOL_COLUMN = "__phospy_fragpipe_gene_symbol"
_ADAPTED_SITE_COLUMN = "__phospy_fragpipe_site"
_ADAPTED_SITE_SEQUENCE_COLUMN = "__phospy_fragpipe_site_sequence"
_ADAPTED_LOCALISATION_COLUMN = "__phospy_fragpipe_localisation_confidence"
_ADAPTED_UNIQUE_FEATURE_ID_COLUMN = "__phospy_fragpipe_feature_id"
_ADAPTED_PEPTIDE_SEQUENCE_COLUMN = "__phospy_fragpipe_peptide_sequence"
_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN = (
    "__phospy_fragpipe_modified_peptide_sequence"
)
_ADAPTED_PEPTIDE_SITE_STRING_COLUMN = "__phospy_fragpipe_peptide_site_string"
_ADAPTED_CANDIDATE_SITES_COLUMN = "fragpipe_ptmprophet_candidate_sites"
_ADAPTED_SITE_PROBABILITIES_COLUMN = "fragpipe_ptmprophet_site_probabilities"
_ADAPTED_AMBIGUOUS_COLUMN = "fragpipe_ptmprophet_ambiguous"
_ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN = "fragpipe_modified_peptide_phospho_count"
_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN = "fragpipe_contaminant"
_FRAGPIPE_DECOY_OUTPUT_COLUMN = "fragpipe_decoy"
_DEFAULT_INTENSITY_PREFIXES = (
    "Intensity ",
    "LFQ Intensity ",
    "MaxLFQ Intensity ",
    "Abundance ",
    "Area ",
)
_PROTEIN_ACCESSION_CANDIDATES = (
    "Protein",
    "Protein ID",
    "Protein IDs",
    "Protein ID(s)",
    "Proteins",
    "Mapped Proteins",
    "Protein Accession",
    "Protein Accession(s)",
)
_GENE_SYMBOL_CANDIDATES = (
    "Gene",
    "Gene Name",
    "Gene names",
    "Genes",
    "Mapped Genes",
)
_PEPTIDE_SEQUENCE_CANDIDATES = (
    "Peptide",
    "Peptide Sequence",
    "Sequence",
)
_MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES = (
    "Modified Peptide",
    "Modified Peptide Sequence",
    "Modified Sequence",
    "Modified Peptide Sequence With Flanking AAs",
)
_PTMPROPHET_PROBABILITY_CANDIDATES = (
    "PTMProphet Probability",
    "PTMProphet Probabilities",
    "PTMProphet Site Probabilities",
    "PTMProphet Localization",
    "PTMProphet Localisation",
    "Localization Probability",
    "Localisation Probability",
    "Best Localization",
    "Best Localisation",
)
_PROTEIN_START_CANDIDATES = (
    "Protein Start",
    "Protein Start Position",
    "Start",
    "Start Position",
    "Peptide Start",
    "Mapped Start",
)
_SITE_CANDIDATES = (
    "Site",
    "Phosphosite",
    "Phospho Site",
    "Modified Site",
)
_SITE_SEQUENCE_CANDIDATES = (
    "Sequence Window",
    "Sequence window",
    "Window Sequence",
    "Site Sequence",
)
_UNIQUE_FEATURE_CANDIDATES = (
    "Spectrum",
    "Spectrum ID",
    "PSM ID",
    "Peptide ID",
    "Index",
    "id",
    "ID",
)
_CONTAMINANT_CANDIDATES = (
    "Contaminant",
    "Potential contaminant",
    "Potential Contaminant",
    "Is Contaminant",
)
_DECOY_CANDIDATES = (
    "Decoy",
    "Reverse",
    "Is Decoy",
    "Protein Decoy",
)
_ROW_ID_CANDIDATES: tuple[str, ...] = ()
_NUMERIC_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_POSITIONED_LOCALISATION_PATTERNS = (
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*([1-9][0-9]*)\s*\(\s*({_NUMERIC_PATTERN})\s*\)"
    ),
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*([1-9][0-9]*)\s*[:=]\s*({_NUMERIC_PATTERN})"
    ),
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*@\s*([1-9][0-9]*)\s*[:=]\s*({_NUMERIC_PATTERN})"
    ),
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*\(\s*([1-9][0-9]*)\s*[,;:]\s*({_NUMERIC_PATTERN})\s*\)"
    ),
)
_RESIDUE_ONLY_LOCALISATION_PATTERN = re.compile(
    rf"(?<![A-Za-z0-9])([STYsty])\s*\(\s*({_NUMERIC_PATTERN})\s*\)"
)
_CONTAMINANT_PREFIXES = ("CON__", "CON_", "CONTAM_", "CONTAMINANT_")
_DECOY_PREFIXES = ("REV__", "REV_", "DECOY__", "DECOY_")


__all__ = [
    "_ADAPTED_AMBIGUOUS_COLUMN",
    "_ADAPTED_CANDIDATE_SITES_COLUMN",
    "_ADAPTED_GENE_SYMBOL_COLUMN",
    "_ADAPTED_LOCALISATION_COLUMN",
    "_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN",
    "_ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN",
    "_ADAPTED_PEPTIDE_SEQUENCE_COLUMN",
    "_ADAPTED_PEPTIDE_SITE_STRING_COLUMN",
    "_ADAPTED_PROTEIN_ACCESSION_COLUMN",
    "_ADAPTED_PROTEIN_ID_COLUMN",
    "_ADAPTED_ROW_ID_COLUMN",
    "_ADAPTED_SITE_COLUMN",
    "_ADAPTED_SITE_PROBABILITIES_COLUMN",
    "_ADAPTED_SITE_SEQUENCE_COLUMN",
    "_ADAPTED_UNIQUE_FEATURE_ID_COLUMN",
    "_CONTAMINANT_CANDIDATES",
    "_CONTAMINANT_PREFIXES",
    "_DECOY_CANDIDATES",
    "_DECOY_PREFIXES",
    "_DEFAULT_INTENSITY_PREFIXES",
    "_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN",
    "_FRAGPIPE_DECOY_OUTPUT_COLUMN",
    "_GENE_SYMBOL_CANDIDATES",
    "_MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES",
    "_PEPTIDE_SEQUENCE_CANDIDATES",
    "_POSITIONED_LOCALISATION_PATTERNS",
    "_PROTEIN_ACCESSION_CANDIDATES",
    "_PROTEIN_START_CANDIDATES",
    "_PTMPROPHET_PROBABILITY_CANDIDATES",
    "_RESIDUE_ONLY_LOCALISATION_PATTERN",
    "_ROW_ID_CANDIDATES",
    "_SITE_CANDIDATES",
    "_SITE_SEQUENCE_CANDIDATES",
    "_UNIQUE_FEATURE_CANDIDATES",
]
