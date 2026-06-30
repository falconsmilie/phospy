"""FragPipe/PTMProphet phosphosite domain conversion helpers."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.io.readers._table_parsing import (
    first_list_token,
    is_missing,
    required_text,
)
from phospy.io.readers.fragpipe.constants import (
    _CONTAMINANT_PREFIXES,
    _DECOY_PREFIXES,
    _POSITIONED_LOCALISATION_PATTERNS,
    _RESIDUE_ONLY_LOCALISATION_PATTERN,
)
from phospy.io.readers.fragpipe.models import (
    _LocalisationCandidate,
    _ProteinSiteCandidate,
    _ResolvedFragPipeColumns,
    _SiteCall,
)
from phospy.science.evidence.modified_peptides import ModifiedResidue
from phospy.science.evidence.multi_site import parse_phospho_site_tokens
from phospy.validation.datasets.fragpipe import (
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE,
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN,
)


def _resolve_site_call(
    row: Mapping[str, object],
    *,
    resolved: _ResolvedFragPipeColumns,
    modified_phospho_sites: tuple[ModifiedResidue, ...],
    ptmprophet_position_reference: str,
    row_position: int,
) -> _SiteCall:
    localisation_candidates = _parse_ptmprophet_localisation_candidates(
        row[resolved.ptmprophet_probabilities],
        modified_phospho_sites=modified_phospho_sites,
        field_name=f"FragPipe {resolved.ptmprophet_probabilities}",
        row_position=row_position,
    )
    protein_start: int | None = None
    protein_start_was_resolved = False

    def resolve_protein_start() -> int | None:
        nonlocal protein_start, protein_start_was_resolved
        if not protein_start_was_resolved:
            protein_start = _resolve_protein_start(
                row,
                resolved=resolved,
                row_position=row_position,
            )
            protein_start_was_resolved = True
        return protein_start

    if resolved.site is not None:
        protein_start_for_localisation = (
            resolve_protein_start()
            if localisation_candidates
            and ptmprophet_position_reference
            == FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE
            else None
        )
        selected = tuple(
            _ProteinSiteCandidate(
                residue=token.residue,
                protein_position=token.position,
                probability=_probability_for_site_token(
                    token.token,
                    candidates=localisation_candidates,
                    protein_start=protein_start_for_localisation,
                    ptmprophet_position_reference=ptmprophet_position_reference,
                    row_position=row_position,
                ),
            )
            for token in parse_phospho_site_tokens(
                row[resolved.site],
                field_name=f"FragPipe {resolved.site} row_position={row_position}",
            )
        )
        all_candidates = _convert_candidates_to_protein_sites(
            localisation_candidates,
            protein_start=protein_start_for_localisation,
            ptmprophet_position_reference=ptmprophet_position_reference,
            row_position=row_position,
        )
        if not all_candidates:
            all_candidates = selected
        _, ambiguous = _select_localised_sites(
            all_candidates,
            phospho_site_count=len(modified_phospho_sites),
        )
        return _site_call_from_candidates(
            selected,
            all_candidates=all_candidates,
            ambiguous=ambiguous,
            phospho_site_count=len(modified_phospho_sites),
        )

    if localisation_candidates:
        protein_start_for_localisation = (
            resolve_protein_start()
            if ptmprophet_position_reference
            == FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE
            else None
        )
        protein_candidates = _convert_candidates_to_protein_sites(
            localisation_candidates,
            protein_start=protein_start_for_localisation,
            ptmprophet_position_reference=ptmprophet_position_reference,
            row_position=row_position,
        )
        selected, ambiguous = _select_localised_sites(
            protein_candidates,
            phospho_site_count=len(modified_phospho_sites),
        )
        return _site_call_from_candidates(
            selected,
            all_candidates=protein_candidates,
            ambiguous=ambiguous,
            phospho_site_count=len(modified_phospho_sites),
        )

    if not modified_phospho_sites:
        raise PhosPyInputError(
            "FragPipe importer could not extract phosphosite evidence from "
            f"modified peptide or PTMProphet localisation; row_position={row_position}"
        )
    fallback_candidates = tuple(
        _protein_site_from_peptide_position(
            residue=site.residue,
            peptide_position=site.position,
            probability=None,
            protein_start=resolve_protein_start(),
            row_position=row_position,
        )
        for site in modified_phospho_sites
    )
    return _site_call_from_candidates(
        fallback_candidates,
        all_candidates=fallback_candidates,
        ambiguous=False,
        phospho_site_count=len(modified_phospho_sites),
    )


def _parse_ptmprophet_localisation_candidates(
    value: object,
    *,
    modified_phospho_sites: tuple[ModifiedResidue, ...],
    field_name: str,
    row_position: int,
) -> tuple[_LocalisationCandidate, ...]:
    if is_missing(value):
        return ()
    if isinstance(value, int | float) and not isinstance(value, bool):
        probability = _parse_probability(
            value, field_name=field_name, row_position=row_position
        )
        return tuple(
            _LocalisationCandidate(
                residue=site.residue,
                position=site.position,
                probability=probability,
            )
            for site in modified_phospho_sites
        )
    if not isinstance(value, str):
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a PTMProphet "
            "localisation string or probability"
        )
    token = value.strip()
    if token == "":
        return ()

    positioned: list[_LocalisationCandidate] = []
    consumed_spans: set[tuple[int, int]] = set()
    for pattern in _POSITIONED_LOCALISATION_PATTERNS:
        for match in pattern.finditer(token):
            span = match.span()
            if span in consumed_spans:
                continue
            consumed_spans.add(span)
            positioned.append(
                _LocalisationCandidate(
                    residue=match.group(1).upper(),
                    position=int(match.group(2)),
                    probability=_parse_probability(
                        match.group(3),
                        field_name=field_name,
                        row_position=row_position,
                    ),
                )
            )
    if positioned:
        return tuple(_deduplicate_localisation_candidates(positioned))

    residue_only_matches = list(_RESIDUE_ONLY_LOCALISATION_PATTERN.finditer(token))
    if residue_only_matches:
        if len(residue_only_matches) > len(modified_phospho_sites):
            raise PhosPyInputError(
                f"{field_name} row_position={row_position} has "
                "residue-only probabilities but more probability tokens than parsed "
                "modified phosphosites"
            )
        candidates: list[_LocalisationCandidate] = []
        search_start = 0
        for match in residue_only_matches:
            residue = match.group(1).upper()
            mapped_site = _next_modified_site_with_residue(
                modified_phospho_sites,
                residue=residue,
                start=search_start,
                row_position=row_position,
                field_name=field_name,
            )
            candidates.append(
                _LocalisationCandidate(
                    residue=residue,
                    position=mapped_site.position,
                    probability=_parse_probability(
                        match.group(2),
                        field_name=field_name,
                        row_position=row_position,
                    ),
                )
            )
            search_start = modified_phospho_sites.index(mapped_site) + 1
        return tuple(candidates)

    try:
        probability = _parse_probability(
            token,
            field_name=field_name,
            row_position=row_position,
        )
    except PhosPyInputError as exc:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} contains malformed "
            f"PTMProphet localisation string {value!r}; expected tokens like "
            "'S3(0.95)' or 'S3:0.95'"
        ) from exc
    return tuple(
        _LocalisationCandidate(
            residue=site.residue,
            position=site.position,
            probability=probability,
        )
        for site in modified_phospho_sites
    )


def _convert_candidates_to_protein_sites(
    candidates: tuple[_LocalisationCandidate, ...],
    *,
    protein_start: int | None,
    ptmprophet_position_reference: str,
    row_position: int,
) -> tuple[_ProteinSiteCandidate, ...]:
    protein_candidates: list[_ProteinSiteCandidate] = []
    for candidate in candidates:
        if candidate.position is None:
            raise PhosPyInputError(
                "FragPipe PTMProphet localisation candidate is missing a position; "
                f"row_position={row_position}"
            )
        if (
            ptmprophet_position_reference
            == FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN
        ):
            protein_candidates.append(
                _ProteinSiteCandidate(
                    residue=candidate.residue,
                    protein_position=candidate.position,
                    probability=candidate.probability,
                )
            )
            continue
        protein_candidates.append(
            _protein_site_from_peptide_position(
                residue=candidate.residue,
                peptide_position=candidate.position,
                probability=candidate.probability,
                protein_start=protein_start,
                row_position=row_position,
            )
        )
    return tuple(_deduplicate_protein_candidates(protein_candidates))


def _select_localised_sites(
    candidates: tuple[_ProteinSiteCandidate, ...],
    *,
    phospho_site_count: int,
) -> tuple[tuple[_ProteinSiteCandidate, ...], bool]:
    if not candidates:
        return (), False
    if phospho_site_count <= 1 and len(candidates) > 1:
        probabilities = [candidate.probability for candidate in candidates]
        if any(probability is None for probability in probabilities):
            return candidates, True
        numeric_probabilities = tuple(
            probability for probability in probabilities if probability is not None
        )
        max_probability = max(numeric_probabilities)
        top_candidates = tuple(
            candidate
            for candidate in candidates
            if candidate.probability is not None
            and math.isclose(
                float(candidate.probability), max_probability, abs_tol=1e-12
            )
        )
        return top_candidates, len(top_candidates) > 1
    if phospho_site_count > 1 and len(candidates) > phospho_site_count:
        return candidates, True
    return candidates, False


def _site_call_from_candidates(
    selected: tuple[_ProteinSiteCandidate, ...],
    *,
    all_candidates: tuple[_ProteinSiteCandidate, ...],
    ambiguous: bool,
    phospho_site_count: int,
) -> _SiteCall:
    if not selected:
        raise PhosPyInputError("FragPipe importer resolved zero phosphosite candidates")
    ordered_selected = tuple(_deduplicate_protein_candidates(selected))
    site_tokens = tuple(candidate.token for candidate in ordered_selected)
    selected_probabilities = [candidate.probability for candidate in ordered_selected]
    if selected_probabilities and all(
        probability is not None for probability in selected_probabilities
    ):
        numeric_selected_probabilities = tuple(
            probability
            for probability in selected_probabilities
            if probability is not None
        )
        localisation_confidence: object = float(min(numeric_selected_probabilities))
    else:
        localisation_confidence = pd.NA
    return _SiteCall(
        site_tokens=site_tokens,
        peptide_site_string=";".join(site_tokens),
        localisation_confidence=localisation_confidence,
        candidate_sites=";".join(
            candidate.token
            for candidate in _deduplicate_protein_candidates(all_candidates)
        ),
        site_probabilities=_format_site_probabilities(all_candidates),
        ambiguous=bool(ambiguous),
        phospho_site_count=int(phospho_site_count),
    )


def _protein_site_from_peptide_position(
    *,
    residue: str,
    peptide_position: int,
    probability: float | None,
    protein_start: int | None,
    row_position: int,
) -> _ProteinSiteCandidate:
    if protein_start is None:
        raise PhosPyInputError(
            "FragPipe importer requires a protein_start column when PTMProphet "
            "positions are peptide-relative and no explicit site column is mapped; "
            f"row_position={row_position}"
        )
    protein_position = int(protein_start) + int(peptide_position) - 1
    if protein_position < 1:
        raise PhosPyInputError(
            f"FragPipe computed invalid protein position {protein_position}; "
            f"row_position={row_position}"
        )
    return _ProteinSiteCandidate(
        residue=residue.upper(),
        protein_position=protein_position,
        probability=probability,
    )


def _probability_for_site_token(
    token: str,
    *,
    candidates: tuple[_LocalisationCandidate, ...],
    protein_start: int | None,
    ptmprophet_position_reference: str,
    row_position: int,
) -> float | None:
    if not candidates:
        return None
    protein_candidates = _convert_candidates_to_protein_sites(
        candidates,
        protein_start=protein_start,
        ptmprophet_position_reference=ptmprophet_position_reference,
        row_position=row_position,
    )
    for candidate in protein_candidates:
        if candidate.token == token:
            return candidate.probability
    if len(protein_candidates) == 1:
        return protein_candidates[0].probability
    return None


def _resolve_protein_start(
    row: Mapping[str, object],
    *,
    resolved: _ResolvedFragPipeColumns,
    row_position: int,
) -> int | None:
    if resolved.protein_start is None:
        return None
    return _normalise_positive_int(
        row[resolved.protein_start],
        field_name=f"FragPipe {resolved.protein_start}",
        row_position=row_position,
    )


def _parse_probability(value: object, *, field_name: str, row_position: int) -> float:
    if isinstance(value, bool):
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} contains boolean probability"
        )
    if isinstance(value, str):
        try:
            probability = float(value.strip())
        except ValueError as exc:
            raise PhosPyInputError(
                f"{field_name} row_position={row_position} contains non-numeric "
                f"probability {value!r}"
            ) from exc
    elif isinstance(value, (int, float)):
        probability = float(value)
    else:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} contains non-numeric "
            f"probability {value!r}"
        )
    if not math.isfinite(probability) or probability < 0.0 or probability > 1.0:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} probability must be in "
            f"[0.0, 1.0]; value={value!r}"
        )
    return float(probability)


def _deduplicate_localisation_candidates(
    candidates: Sequence[_LocalisationCandidate],
) -> tuple[_LocalisationCandidate, ...]:
    values: dict[tuple[str, int | None], _LocalisationCandidate] = {}
    for candidate in candidates:
        key = (candidate.residue.upper(), candidate.position)
        current = values.get(key)
        if current is None:
            values[key] = candidate
            continue
        if current.probability is None:
            values[key] = candidate
            continue
        if candidate.probability is not None and candidate.probability > (
            current.probability
        ):
            values[key] = candidate
    return tuple(values.values())


def _deduplicate_protein_candidates(
    candidates: Sequence[_ProteinSiteCandidate],
) -> tuple[_ProteinSiteCandidate, ...]:
    values: dict[str, _ProteinSiteCandidate] = {}
    for candidate in candidates:
        current = values.get(candidate.token)
        if current is None:
            values[candidate.token] = candidate
            continue
        if current.probability is None:
            values[candidate.token] = candidate
            continue
        if candidate.probability is not None and candidate.probability > (
            current.probability
        ):
            values[candidate.token] = candidate
    return tuple(values.values())


def _format_site_probabilities(
    candidates: tuple[_ProteinSiteCandidate, ...],
) -> str:
    parts: list[str] = []
    for candidate in _deduplicate_protein_candidates(candidates):
        if candidate.probability is None:
            parts.append(f"{candidate.token}:NA")
            continue
        parts.append(f"{candidate.token}:{float(candidate.probability):.6g}")
    return ";".join(parts)


def _next_modified_site_with_residue(
    sites: tuple[ModifiedResidue, ...],
    *,
    residue: str,
    start: int,
    row_position: int,
    field_name: str,
) -> ModifiedResidue:
    residue = residue.upper()
    for site in sites[start:]:
        if site.residue == residue:
            return site
    raise PhosPyInputError(
        f"{field_name} row_position={row_position} could not align residue-only "
        f"PTMProphet probability for residue {residue!r} to the modified peptide"
    )


def _parse_protein_accession(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    token = first_list_token(value, field_name=field_name, row_position=row_position)
    cleaned = _strip_protein_prefixes(token)
    parts = cleaned.split("|")
    if len(parts) >= 3 and parts[0].strip().lower() in {"sp", "tr"}:
        accession = parts[1].strip()
    else:
        accession = cleaned.split()[0].strip()
    if accession == "":
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} did not contain a protein "
            "accession after parsing"
        )
    return accession


def _strip_protein_prefixes(value: str) -> str:
    token = value.strip()
    changed = True
    while changed:
        changed = False
        upper = token.upper()
        for prefix in (*_CONTAMINANT_PREFIXES, *_DECOY_PREFIXES):
            if upper.startswith(prefix):
                token = token[len(prefix) :].strip()
                changed = True
                break
    return token


def _normalise_positive_int(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> int:
    token = required_text(value, field_name=field_name, row_position=row_position)
    try:
        numeric = float(token)
    except ValueError as exc:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a positive integer; "
            f"value={value!r}"
        ) from exc
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a positive integer; "
            f"value={value!r}"
        )
    integer = int(numeric)
    if integer < 1:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a positive integer; "
            f"value={value!r}"
        )
    return integer


__all__ = [
    "_parse_protein_accession",
    "_resolve_protein_start",
    "_resolve_site_call",
]
