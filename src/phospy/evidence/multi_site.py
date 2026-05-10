"""Shared multi-site phosphopeptide parsing and resolution models."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.site_ids import canonicalize_site_components

MULTI_SITE_POLICY_KEEP_JOINT = "keep_joint"
MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING = "exclude_from_sequence_scoring"
MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY = "first_site_compatibility"
MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT = "split_equal_weight"
MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY = "split_with_uncertainty"
MULTI_SITE_POLICY_ERROR = "error"
SUPPORTED_MULTI_SITE_POLICIES: tuple[str, ...] = (
    MULTI_SITE_POLICY_KEEP_JOINT,
    MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY,
    MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
    MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY,
    MULTI_SITE_POLICY_ERROR,
)

_SITE_TOKEN_PATTERN = re.compile(r"^\s*([STYsty])([1-9][0-9]*)\s*$")
_SITE_TOKEN_SPLIT_PATTERN = re.compile(r"\s*[,;]\s*")


@dataclass(frozen=True, slots=True)
class PhosphoSiteToken:
    """One parsed phosphosite token such as ``S1246``."""

    residue: str
    position: int

    def __post_init__(self) -> None:
        residue = str(self.residue).strip().upper()
        if residue not in {"S", "T", "Y"}:
            raise PhosPyInputError(
                "phospho_site_token.residue must be one of: 'S', 'T', 'Y'"
            )
        if isinstance(self.position, bool) or not isinstance(self.position, int):
            raise PhosPyInputError("phospho_site_token.position must be an int >= 1")
        if self.position < 1:
            raise PhosPyInputError("phospho_site_token.position must be an int >= 1")
        object.__setattr__(self, "residue", residue)
        object.__setattr__(self, "position", int(self.position))

    @property
    def token(self) -> str:
        return f"{self.residue}{self.position}"


@dataclass(frozen=True, slots=True)
class MultiSiteObservation:
    """One peptide-row multi-site interpretation before downstream resolution."""

    peptide_row_id: str
    gene_symbol: str
    site_string: str
    tokens: tuple[PhosphoSiteToken, ...]
    joint_site_id: str
    is_multi_site: bool
    declared_multi_site: bool


@dataclass(frozen=True, slots=True)
class MultiSiteHandlingConfig:
    """Config for explicit multi-site behavior by workflow purpose."""

    statistical_modeling_policy: str = MULTI_SITE_POLICY_KEEP_JOINT
    kinase_sequence_scoring_policy: str = (
        MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    )

    def __post_init__(self) -> None:
        _validate_policy(
            self.statistical_modeling_policy,
            field_name="multi_site_handling_config.statistical_modeling_policy",
        )
        _validate_policy(
            self.kinase_sequence_scoring_policy,
            field_name="multi_site_handling_config.kinase_sequence_scoring_policy",
        )


def parse_phospho_site_tokens(
    site_string: object,
    *,
    field_name: str,
) -> tuple[PhosphoSiteToken, ...]:
    if not isinstance(site_string, str):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    stripped = site_string.strip()
    if stripped == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    raw_tokens = _SITE_TOKEN_SPLIT_PATTERN.split(stripped)
    if any(token.strip() == "" for token in raw_tokens):
        raise PhosPyInputError(
            f"{field_name} contains empty site tokens: {site_string!r}"
        )
    parsed: list[PhosphoSiteToken] = []
    for token in raw_tokens:
        match = _SITE_TOKEN_PATTERN.fullmatch(token)
        if match is None:
            raise PhosPyInputError(
                f"{field_name} contains malformed site token {token!r}; expected "
                "tokens like 'S123' separated by ',' or ';'"
            )
        parsed.append(
            PhosphoSiteToken(
                residue=match.group(1).upper(),
                position=int(match.group(2)),
            )
        )
    deduplicated = tuple(dict.fromkeys(item.token for item in parsed))
    if len(deduplicated) != len(parsed):
        raise PhosPyInputError(
            f"{field_name} contains duplicate site tokens: {site_string!r}"
        )
    return tuple(parsed)


def build_multi_site_observation(
    *,
    peptide_row_id: str,
    gene_symbol: str,
    site_string: str,
    declared_multi_site: bool,
    field_name: str,
) -> MultiSiteObservation:
    tokens = parse_phospho_site_tokens(site_string, field_name=field_name)
    is_multi_site = len(tokens) > 1
    if bool(declared_multi_site) != bool(is_multi_site):
        raise PhosPyInputError(
            f"{field_name} parsed {len(tokens)} site token(s) from {site_string!r}, "
            "but multi_site disagrees"
        )
    joint_token = ",".join(item.token for item in tokens)
    joint_site_id = canonicalize_site_components(
        gene_symbol,
        joint_token,
        field_name=f"{field_name}.joint_site_id",
        error_type=PhosPyInputError,
    )
    return MultiSiteObservation(
        peptide_row_id=str(peptide_row_id),
        gene_symbol=str(gene_symbol),
        site_string=str(site_string),
        tokens=tokens,
        joint_site_id=joint_site_id,
        is_multi_site=is_multi_site,
        declared_multi_site=bool(declared_multi_site),
    )


def resolve_observation_site_rows(
    *,
    observation: MultiSiteObservation,
    policy: str,
) -> tuple[dict[str, object], ...]:
    _validate_policy(policy, field_name="multi_site_policy")
    if observation.is_multi_site and policy == MULTI_SITE_POLICY_ERROR:
        raise PhosPyInputError(
            "multi-site observation cannot be resolved under policy='error': "
            f"peptide_row_id={observation.peptide_row_id!r}, "
            f"site_string={observation.site_string!r}"
        )
    if observation.is_multi_site and policy == MULTI_SITE_POLICY_KEEP_JOINT:
        return (
            _resolved_row(
                observation,
                observation.joint_site_id,
                1.0,
                False,
                policy=policy,
            ),
        )
    if (
        observation.is_multi_site
        and policy == MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    ):
        return ()
    if (
        observation.is_multi_site
        and policy == MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY
    ):
        first_site = canonicalize_site_components(
            observation.gene_symbol,
            observation.tokens[0].token,
            field_name="multi_site_policy.first_site_compatibility",
            error_type=PhosPyInputError,
        )
        return (_resolved_row(observation, first_site, 1.0, False, policy=policy),)
    if observation.is_multi_site and policy in {
        MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
        MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY,
    }:
        uncertainty = policy == MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY
        weight = 1.0 / float(len(observation.tokens))
        rows: list[dict[str, object]] = []
        for token in observation.tokens:
            split_site_id = canonicalize_site_components(
                observation.gene_symbol,
                token.token,
                field_name="multi_site_policy.split",
                error_type=PhosPyInputError,
            )
            rows.append(
                _resolved_row(
                    observation,
                    split_site_id,
                    weight,
                    uncertainty,
                    policy=policy,
                )
            )
        return tuple(rows)

    single_site_id = canonicalize_site_components(
        observation.gene_symbol,
        observation.tokens[0].token,
        field_name="multi_site_policy.single_site",
        error_type=PhosPyInputError,
    )
    return (_resolved_row(observation, single_site_id, 1.0, False, policy=policy),)


def resolve_site_mapping_frame(
    *,
    observations: Iterable[MultiSiteObservation],
    policy: str,
) -> pd.DataFrame:
    resolved_rows: list[dict[str, object]] = []
    for observation in observations:
        resolved_rows.extend(
            resolve_observation_site_rows(observation=observation, policy=policy)
        )
    if not resolved_rows:
        return pd.DataFrame(
            columns=(
                "peptide_row_id",
                "site_id",
                "mapping_weight",
                "mapping_uncertainty",
                "multi_site_policy",
                "is_multi_site",
            )
        )
    frame = pd.DataFrame(resolved_rows)
    return frame.loc[
        :,
        [
            "peptide_row_id",
            "site_id",
            "mapping_weight",
            "mapping_uncertainty",
            "multi_site_policy",
            "is_multi_site",
        ],
    ].copy(deep=True)


def _resolved_row(
    observation: MultiSiteObservation,
    site_id: str,
    weight: float,
    uncertainty: bool,
    *,
    policy: str,
) -> dict[str, object]:
    return {
        "peptide_row_id": observation.peptide_row_id,
        "site_id": site_id,
        "mapping_weight": float(weight),
        "mapping_uncertainty": bool(uncertainty),
        "multi_site_policy": str(policy),
        "is_multi_site": bool(observation.is_multi_site),
    }


def _validate_policy(policy: object, *, field_name: str) -> None:
    if not isinstance(policy, str) or policy not in SUPPORTED_MULTI_SITE_POLICIES:
        allowed = ", ".join(repr(item) for item in SUPPORTED_MULTI_SITE_POLICIES)
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")


__all__ = [
    "MULTI_SITE_POLICY_ERROR",
    "MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING",
    "MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY",
    "MULTI_SITE_POLICY_KEEP_JOINT",
    "MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT",
    "MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY",
    "SUPPORTED_MULTI_SITE_POLICIES",
    "MultiSiteHandlingConfig",
    "MultiSiteObservation",
    "PhosphoSiteToken",
    "build_multi_site_observation",
    "parse_phospho_site_tokens",
    "resolve_observation_site_rows",
    "resolve_site_mapping_frame",
]
