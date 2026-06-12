"""Shared modified-peptide parsing helpers for importer adapters."""

from __future__ import annotations

import re
from dataclasses import dataclass

from phospy.errors.input import PhosPyInputError

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_PHOSPHO_RESIDUES = frozenset("STY")
_MODIFIED_RESIDUE_TOKEN_PATTERN = re.compile(
    r"(?i)(?:^|[^A-Z])(?:p|phospho(?:rylation)?)?\s*([STY])(?:$|[^A-Z])"
)
_MASS_TOKEN_PATTERN = re.compile(r"(?<!\d)79(?:\.9(?:6|7)\d*)?(?!\d)")


@dataclass(frozen=True, slots=True)
class ModifiedResidue:
    """One peptide-relative modified residue parsed from a peptide string."""

    residue: str
    position: int
    modification: str

    def __post_init__(self) -> None:
        residue = str(self.residue).strip().upper()
        if residue not in _AMINO_ACIDS:
            raise PhosPyInputError(
                "modified_residue.residue must be a one-letter amino-acid code"
            )
        if isinstance(self.position, bool) or not isinstance(self.position, int):
            raise PhosPyInputError("modified_residue.position must be an int >= 1")
        if self.position < 1:
            raise PhosPyInputError("modified_residue.position must be an int >= 1")
        modification = str(self.modification).strip()
        if modification == "":
            raise PhosPyInputError(
                "modified_residue.modification must be a non-empty string"
            )
        object.__setattr__(self, "residue", residue)
        object.__setattr__(self, "position", int(self.position))
        object.__setattr__(self, "modification", modification)


@dataclass(frozen=True, slots=True)
class ModifiedPeptideParseResult:
    """Parsed modified peptide with unmodified sequence and phospho positions."""

    sequence: str
    modifications: tuple[ModifiedResidue, ...]

    @property
    def phospho_sites(self) -> tuple[ModifiedResidue, ...]:
        return tuple(
            modification
            for modification in self.modifications
            if is_phosphorylation_label(modification.modification)
        )


def parse_modified_peptide_sequence(
    value: object,
    *,
    field_name: str,
) -> ModifiedPeptideParseResult:
    """Parse common modified-peptide notations without tool-specific semantics.

    The parser intentionally returns peptide-relative positions only. Importers
    remain responsible for translating those positions into protein-scoped
    phosphosite identifiers using their own upstream columns.
    """

    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    source = _strip_peptide_flanks(value.strip())
    if source == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")

    sequence: list[str] = []
    modifications: list[ModifiedResidue] = []
    pending_prefix_modifications: list[str] = []
    index = 0
    while index < len(source):
        char = source[index]
        if char in {"_", "-", " "}:
            index += 1
            continue
        if char in "([{":
            content, index = _read_balanced_group(source, index, field_name=field_name)
            inline_sites = _parse_inline_modified_residues(content)
            if inline_sites:
                for residue, modification in inline_sites:
                    sequence.append(residue)
                    modifications.append(
                        ModifiedResidue(
                            residue=residue,
                            position=len(sequence),
                            modification=modification,
                        )
                    )
                continue
            pending_prefix_modifications.append(content)
            continue
        if char.isalpha():
            if char == "p" and index + 1 < len(source):
                next_residue = source[index + 1].upper()
                if next_residue in _PHOSPHO_RESIDUES:
                    sequence.append(next_residue)
                    modifications.append(
                        ModifiedResidue(
                            residue=next_residue,
                            position=len(sequence),
                            modification="p",
                        )
                    )
                    index += 2
                    continue
            residue = char.upper()
            if residue not in _AMINO_ACIDS:
                raise PhosPyInputError(
                    f"{field_name} contains unsupported amino-acid token {char!r}"
                )
            sequence.append(residue)
            position = len(sequence)
            for modification in pending_prefix_modifications:
                if is_phosphorylation_label(modification) and residue not in (
                    _PHOSPHO_RESIDUES
                ):
                    raise PhosPyInputError(
                        f"{field_name} assigns phosphorylation to non-STY residue "
                        f"{residue!r} at peptide_position={position}"
                    )
                modifications.append(
                    ModifiedResidue(
                        residue=residue,
                        position=position,
                        modification=modification,
                    )
                )
            pending_prefix_modifications.clear()
            index += 1
            while index < len(source) and source[index] in "([{":
                content, index = _read_balanced_group(
                    source,
                    index,
                    field_name=field_name,
                )
                inline_sites = _parse_inline_modified_residues(content)
                if inline_sites:
                    for inline_residue, inline_modification in inline_sites:
                        sequence.append(inline_residue)
                        modifications.append(
                            ModifiedResidue(
                                residue=inline_residue,
                                position=len(sequence),
                                modification=inline_modification,
                            )
                        )
                    continue
                if is_phosphorylation_label(content) and residue not in (
                    _PHOSPHO_RESIDUES
                ):
                    next_residue = _next_residue(source, index)
                    if next_residue in _PHOSPHO_RESIDUES:
                        pending_prefix_modifications.append(content)
                        continue
                    raise PhosPyInputError(
                        f"{field_name} assigns phosphorylation to non-STY residue "
                        f"{residue!r} at peptide_position={position}"
                    )
                modifications.append(
                    ModifiedResidue(
                        residue=residue,
                        position=position,
                        modification=content,
                    )
                )
            continue
        raise PhosPyInputError(f"{field_name} contains unsupported character {char!r}")

    if pending_prefix_modifications:
        raise PhosPyInputError(
            f"{field_name} contains modification annotation without a residue"
        )
    if not sequence:
        raise PhosPyInputError(f"{field_name} must contain an amino-acid sequence")
    return ModifiedPeptideParseResult(
        sequence="".join(sequence),
        modifications=tuple(modifications),
    )


def is_phosphorylation_label(value: object) -> bool:
    """Return whether a modification annotation denotes phosphorylation."""

    if value is None:
        return False
    token = str(value).strip()
    if token == "":
        return False
    normalised = re.sub(r"[\s_\-]+", "", token).lower()
    if normalised in {"p", "ph", "phospho", "phosphorylation"}:
        return True
    if normalised in {"ps", "pt", "py"}:
        return True
    if "phospho" in normalised:
        return True
    if _MASS_TOKEN_PATTERN.search(normalised):
        return True
    return False


def _strip_peptide_flanks(value: str) -> str:
    parts = value.split(".")
    if len(parts) >= 3 and len(parts[0].strip()) == 1 and len(parts[-1].strip()) == 1:
        return ".".join(parts[1:-1]).strip()
    return value.strip()


def _read_balanced_group(
    value: str,
    start: int,
    *,
    field_name: str,
) -> tuple[str, int]:
    opener = value[start]
    closer = {"(": ")", "[": "]", "{": "}"}[opener]
    depth = 1
    index = start + 1
    while index < len(value):
        char = value[index]
        if char == opener:
            depth += 1
        elif char == closer:
            depth -= 1
            if depth == 0:
                content = value[start + 1 : index].strip()
                if content == "":
                    raise PhosPyInputError(
                        f"{field_name} contains an empty modification annotation"
                    )
                return content, index + 1
        index += 1
    raise PhosPyInputError(f"{field_name} contains an unclosed modification annotation")


def _parse_inline_modified_residues(content: str) -> tuple[tuple[str, str], ...]:
    parsed: list[tuple[str, str]] = []
    tokens = [token.strip() for token in re.split(r"\s*[,;]\s*", content)]
    for token in tokens:
        if token == "":
            return ()
        match = _MODIFIED_RESIDUE_TOKEN_PATTERN.search(token)
        if match is None or not is_phosphorylation_label(token):
            return ()
        parsed.append((match.group(1).upper(), token))
    return tuple(parsed)


def _next_residue(value: str, start: int) -> str | None:
    index = start
    while index < len(value):
        char = value[index]
        if char in {"_", "-", " "}:
            index += 1
            continue
        if char.isalpha():
            return char.upper()
        return None
    return None


__all__ = [
    "ModifiedPeptideParseResult",
    "ModifiedResidue",
    "is_phosphorylation_label",
    "parse_modified_peptide_sequence",
]
