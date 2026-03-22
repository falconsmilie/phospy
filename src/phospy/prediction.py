from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .scoring import KinaseScoringResult


@dataclass(slots=True)
class AdaptiveSamplingIterationTrace:
    iteration_index: int
    labels: pd.Series
    probabilities: pd.DataFrame
    positive_weights: pd.Series | None
    negative_weights: pd.Series | None
    sampled_positive_sites: list[str]
    sampled_negative_sites: list[str]


@dataclass(slots=True)
class AdaptiveSamplingEnsembleTrace:
    ensemble_index: int
    initial_negative_sites: list[str]
    iterations: list[AdaptiveSamplingIterationTrace]
    final_prediction_probabilities: pd.DataFrame
    final_top_sites: list[str]


@dataclass(slots=True)
class KinasePredictionDebugTrace:
    kinase: str
    candidate_substrates: list[str]
    negative_pool_sites: list[str]
    ensemble_traces: list[AdaptiveSamplingEnsembleTrace]


@dataclass(slots=True)
class KinasePredictionResult:
    pred_matrix: pd.DataFrame
    substrate_list: dict[str, list[str]]
    debug_traces: dict[str, KinasePredictionDebugTrace] | None = None


@dataclass(slots=True)
class SamplingTraceOverrideEnsemble:
    initial_negative_sites: list[str] | None
    iteration_sample_sites: dict[int, dict[int, list[str]]]

    def get_iteration_sample_sites(
        self, iteration_index: int, class_label: int
    ) -> list[str] | None:
        return self.iteration_sample_sites.get(iteration_index, {}).get(class_label)


@dataclass(slots=True)
class PredictionSamplingTrace:
    ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]]

    @classmethod
    def from_trace_directory(cls, trace_dir: str | Path) -> PredictionSamplingTrace:
        path = Path(trace_dir)
        initial_path = path / "trace_initial_negatives.csv"
        samples_path = path / "trace_iteration_samples.csv"
        if not initial_path.exists() and not samples_path.exists():
            msg = (
                "sampling trace directory must contain trace_initial_negatives.csv "
                "and/or trace_iteration_samples.csv"
            )
            raise ValueError(msg)

        ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]] = {}

        if initial_path.exists():
            initial_df = pd.read_csv(initial_path)
            required_initial_cols = {"kinase", "ensemble", "draw", "site"}
            if not required_initial_cols.issubset(initial_df.columns):
                missing = sorted(required_initial_cols.difference(initial_df.columns))
                raise ValueError(
                    "trace_initial_negatives.csv is missing required columns: "
                    + ", ".join(missing)
                )
            initial_df = initial_df.sort_values(
                ["kinase", "ensemble", "draw"], kind="mergesort"
            )
            for (kinase, ensemble), group in initial_df.groupby(
                ["kinase", "ensemble"], sort=False
            ):
                ensemble_map = ensembles_by_kinase.setdefault(str(kinase), {})
                ensemble_map[int(ensemble)] = SamplingTraceOverrideEnsemble(
                    initial_negative_sites=group.loc[:, "site"].astype(str).tolist(),
                    iteration_sample_sites={},
                )

        if samples_path.exists():
            samples_df = pd.read_csv(samples_path)
            required_sample_cols = {
                "kinase",
                "ensemble",
                "iteration",
                "class_label",
                "draw",
                "site",
            }
            if not required_sample_cols.issubset(samples_df.columns):
                missing = sorted(required_sample_cols.difference(samples_df.columns))
                raise ValueError(
                    "trace_iteration_samples.csv is missing required columns: "
                    + ", ".join(missing)
                )
            samples_df = samples_df.sort_values(
                ["kinase", "ensemble", "iteration", "class_label", "draw"],
                kind="mergesort",
            )
            for (kinase, ensemble, iteration, class_label), group in samples_df.groupby(
                ["kinase", "ensemble", "iteration", "class_label"], sort=False
            ):
                ensemble_map = ensembles_by_kinase.setdefault(str(kinase), {})
                ensemble_override = ensemble_map.setdefault(
                    int(ensemble),
                    SamplingTraceOverrideEnsemble(
                        initial_negative_sites=None,
                        iteration_sample_sites={},
                    ),
                )
                iteration_map = ensemble_override.iteration_sample_sites.setdefault(
                    int(iteration), {}
                )
                iteration_map[int(class_label)] = (
                    group.loc[:, "site"].astype(str).tolist()
                )

        return cls(ensembles_by_kinase=ensembles_by_kinase)

    def get_ensemble_override(
        self, kinase: str, ensemble_index: int
    ) -> SamplingTraceOverrideEnsemble | None:
        return self.ensembles_by_kinase.get(kinase, {}).get(ensemble_index)

    def subset_kinases(self, kinases: list[str] | set[str]) -> PredictionSamplingTrace:
        kinase_set = {str(kinase) for kinase in kinases}
        return PredictionSamplingTrace(
            ensembles_by_kinase={
                kinase: ensemble_map
                for kinase, ensemble_map in self.ensembles_by_kinase.items()
                if kinase in kinase_set
            }
        )


class KinasePredictor:
    """Predict kinase-substrate relationships from phosphosite score matrices.

    This is a narrow native Python port of PhosR's score-to-prediction seam.
    It mirrors the broad structure of ``kinaseSubstratePred()``: candidate
    substrates are selected from the combined score matrix, then an ensemble of
    adaptive SVM models is used to produce a kinase prediction matrix.
    """

    def __init__(self, kernel: str = "rbf") -> None:
        self.kernel = kernel

    def predict(
        self,
        combined_scores: pd.DataFrame,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        capture_debug_trace: bool = False,
        debug_kinases: list[str] | None = None,
        debug_top_n: int = 10,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
    ) -> KinasePredictionResult:
        _validate_positive_int(ensemble_size, name="ensemble_size")
        _validate_positive_int(top, name="top")
        _validate_positive_int(inclusion, name="inclusion")
        _validate_positive_int(n_iterations, name="n_iterations")
        _validate_positive_int(debug_top_n, name="debug_top_n")

        sampling_trace_obj = _coerce_sampling_trace(sampling_trace)

        substrate_list = build_candidate_substrate_list(
            combined_scores,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
        )
        if not substrate_list:
            empty = pd.DataFrame(index=combined_scores.index.copy(), dtype=float)
            return KinasePredictionResult(
                pred_matrix=empty,
                substrate_list={},
                debug_traces={} if capture_debug_trace else None,
            )

        rng = np.random.default_rng(random_state)
        feature_mat = combined_scores.astype(float)
        pred_matrix = pd.DataFrame(
            0.0,
            index=feature_mat.index.copy(),
            columns=list(substrate_list),
        )

        traced_kinases = (
            set(substrate_list)
            if capture_debug_trace and debug_kinases is None
            else set(debug_kinases or [])
        )
        debug_traces: dict[str, KinasePredictionDebugTrace] | None = (
            {} if capture_debug_trace else None
        )

        for kinase, substrates in substrate_list.items():
            positive_train = feature_mat.loc[substrates, :]
            negative_pool = feature_mat.loc[
                ~feature_mat.index.isin(substrates),
                :,
            ]
            if negative_pool.empty:
                msg = (
                    f"No negative pool sites available to train predictor for {kinase}"
                )
                raise ValueError(msg)

            if (
                capture_debug_trace
                and kinase in traced_kinases
                and debug_traces is not None
            ):
                debug_traces[kinase] = KinasePredictionDebugTrace(
                    kinase=kinase,
                    candidate_substrates=list(substrates),
                    negative_pool_sites=negative_pool.index.tolist(),
                    ensemble_traces=[],
                )

            for ensemble_idx in range(ensemble_size):
                ensemble_index = ensemble_idx + 1
                ensemble_override = None
                if sampling_trace_obj is not None:
                    ensemble_override = sampling_trace_obj.get_ensemble_override(
                        kinase=kinase,
                        ensemble_index=ensemble_index,
                    )

                if (
                    ensemble_override is not None
                    and ensemble_override.initial_negative_sites is not None
                ):
                    negative_sites = _validate_override_sites(
                        available_sites=negative_pool.index,
                        sampled_sites=ensemble_override.initial_negative_sites,
                        expected_size=len(positive_train),
                        context=(
                            f"initial negatives for kinase={kinase}, "
                            f"ensemble={ensemble_index}"
                        ),
                    )
                else:
                    negative_indices = rng.choice(
                        negative_pool.index.to_numpy(),
                        size=len(positive_train),
                        replace=len(negative_pool) < len(positive_train),
                    )
                    negative_sites = list(negative_indices.tolist())

                negative_train = negative_pool.loc[negative_sites, :]
                train_mat = pd.concat([positive_train, negative_train], axis=0)
                labels = np.concatenate(
                    [
                        np.repeat(1, len(positive_train)),
                        np.repeat(2, len(negative_train)),
                    ]
                )

                if (
                    capture_debug_trace
                    and kinase in traced_kinases
                    and debug_traces is not None
                ):
                    series, ensemble_trace = _multi_ada_sampling(
                        train_mat=train_mat,
                        test_mat=feature_mat,
                        labels=labels,
                        kernel=self.kernel,
                        n_iterations=n_iterations,
                        rng=rng,
                        capture_trace=True,
                        ensemble_index=ensemble_index,
                        initial_negative_sites=negative_sites,
                        debug_top_n=debug_top_n,
                        sampling_override=ensemble_override,
                    )
                    if ensemble_trace is not None:
                        debug_traces[kinase].ensemble_traces.append(ensemble_trace)
                else:
                    series, _ = _multi_ada_sampling(
                        train_mat=train_mat,
                        test_mat=feature_mat,
                        labels=labels,
                        kernel=self.kernel,
                        n_iterations=n_iterations,
                        rng=rng,
                        capture_trace=False,
                        ensemble_index=ensemble_index,
                        initial_negative_sites=negative_sites,
                        debug_top_n=debug_top_n,
                        sampling_override=ensemble_override,
                    )

                pred_matrix.loc[:, kinase] += series

        pred_matrix /= float(ensemble_size)
        return KinasePredictionResult(
            pred_matrix=pred_matrix,
            substrate_list=substrate_list,
            debug_traces=debug_traces,
        )

    def predict_from_scoring_result(
        self,
        scoring_result: KinaseScoringResult,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        allow_profile_only_fallback: bool = False,
        capture_debug_trace: bool = False,
        debug_kinases: list[str] | None = None,
        debug_top_n: int = 10,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
    ) -> KinasePredictionResult:
        if scoring_result.combined_scores is not None:
            feature_mat = scoring_result.combined_scores
        elif allow_profile_only_fallback:
            feature_mat = scoring_result.profile_scores
        else:
            msg = (
                "scoring_result does not contain combined_scores; pass "
                "allow_profile_only_fallback=True to use profile_scores instead"
            )
            raise ValueError(msg)

        return self.predict(
            combined_scores=feature_mat,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            capture_debug_trace=capture_debug_trace,
            debug_kinases=debug_kinases,
            debug_top_n=debug_top_n,
            sampling_trace=sampling_trace,
        )


def build_candidate_substrate_list(
    combined_scores: pd.DataFrame,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
) -> dict[str, list[str]]:
    """Select candidate kinase substrates from the combined score matrix."""

    _validate_positive_int(top, name="top")
    _validate_positive_int(inclusion, name="inclusion")

    substrate_list: dict[str, list[str]] = {}
    for kinase in combined_scores.columns:
        selected = (
            combined_scores.loc[:, kinase]
            .sort_values(ascending=False, kind="mergesort")
            .head(top)
        )
        sites = selected.loc[selected > score_threshold].index.tolist()
        if len(sites) >= inclusion:
            substrate_list[kinase] = sites
    return substrate_list


def _multi_ada_sampling(
    train_mat: pd.DataFrame,
    test_mat: pd.DataFrame,
    labels: np.ndarray,
    kernel: str,
    n_iterations: int,
    rng: np.random.Generator,
    capture_trace: bool,
    ensemble_index: int,
    initial_negative_sites: list[str],
    debug_top_n: int,
    sampling_override: SamplingTraceOverrideEnsemble | None,
) -> tuple[pd.Series, AdaptiveSamplingEnsembleTrace | None]:
    StandardScaler, SVC = _require_sklearn()

    base_x = train_mat.to_numpy(dtype=float)
    base_y = np.asarray(labels, dtype=int)
    base_index = train_mat.index.copy()
    current_x = base_x
    current_y = base_y
    model = None
    iteration_traces: list[AdaptiveSamplingIterationTrace] = []

    for iteration_index in range(1, n_iterations + 1):
        model = _make_svm(
            StandardScaler=StandardScaler, SVC=SVC, kernel=kernel, rng=rng
        )
        model.fit(current_x, current_y)
        prob_mat = model.predict_proba(base_x)
        prob_df = pd.DataFrame(
            prob_mat,
            index=base_index,
            columns=[str(class_label) for class_label in model.classes_],
        )
        label_series = pd.Series(base_y, index=base_index, dtype=int)

        resampled_x: list[np.ndarray] = []
        resampled_y: list[np.ndarray] = []
        sampled_sites_by_class: dict[int, list[str]] = {}
        weights_by_class: dict[int, pd.Series | None] = {}
        for class_idx, class_label in enumerate(model.classes_):
            class_mask = base_y == class_label
            class_x = base_x[class_mask]
            class_index = base_index[class_mask]
            class_prob = prob_mat[class_mask, class_idx]
            sample_prob = _normalize_probabilities(class_prob)
            weights_by_class[int(class_label)] = (
                pd.Series(sample_prob, index=class_index, dtype=float)
                if sample_prob is not None
                else None
            )
            override_sites = None
            if sampling_override is not None:
                override_sites = sampling_override.get_iteration_sample_sites(
                    iteration_index=iteration_index,
                    class_label=int(class_label),
                )
            if override_sites is not None:
                sampled_idx = _resolve_sampled_site_positions(
                    available_sites=class_index,
                    sampled_sites=override_sites,
                    expected_size=class_x.shape[0],
                    context=(
                        f"iteration={iteration_index}, ensemble={ensemble_index}, "
                        f"class_label={int(class_label)}"
                    ),
                )
            else:
                sampled_idx = rng.choice(
                    class_x.shape[0],
                    size=class_x.shape[0],
                    replace=True,
                    p=sample_prob,
                )
            sampled_sites = class_index[sampled_idx].tolist()
            sampled_sites_by_class[int(class_label)] = sampled_sites
            resampled_x.append(class_x[sampled_idx])
            resampled_y.append(np.repeat(class_label, class_x.shape[0]))

        if capture_trace:
            iteration_traces.append(
                AdaptiveSamplingIterationTrace(
                    iteration_index=iteration_index,
                    labels=label_series,
                    probabilities=prob_df,
                    positive_weights=weights_by_class.get(1),
                    negative_weights=weights_by_class.get(2),
                    sampled_positive_sites=sampled_sites_by_class.get(1, []),
                    sampled_negative_sites=sampled_sites_by_class.get(2, []),
                )
            )

        current_x = np.vstack(resampled_x)
        current_y = np.concatenate(resampled_y)

    if model is None:
        msg = "n_iterations must be at least 1"
        raise ValueError(msg)

    pred = model.predict_proba(test_mat.to_numpy(dtype=float))
    pred_df = pd.DataFrame(
        pred,
        index=test_mat.index.copy(),
        columns=[str(class_label) for class_label in model.classes_],
    )
    positive_idx = np.flatnonzero(model.classes_ == 1)
    if len(positive_idx) != 1:
        msg = f"Expected exactly one positive class labelled 1; found {model.classes_.tolist()}"
        raise ValueError(msg)
    positive_series = pd.Series(
        pred[:, positive_idx[0]], index=test_mat.index.copy(), dtype=float
    )

    ensemble_trace = None
    if capture_trace:
        final_top_sites = (
            positive_series.sort_values(ascending=False)
            .head(debug_top_n)
            .index.tolist()
        )
        ensemble_trace = AdaptiveSamplingEnsembleTrace(
            ensemble_index=ensemble_index,
            initial_negative_sites=list(initial_negative_sites),
            iterations=iteration_traces,
            final_prediction_probabilities=pred_df,
            final_top_sites=final_top_sites,
        )

    return positive_series, ensemble_trace


def _make_svm(
    *, StandardScaler: type, SVC: type, kernel: str, rng: np.random.Generator
):
    from sklearn.pipeline import make_pipeline

    return make_pipeline(
        StandardScaler(),
        SVC(
            kernel=kernel,
            probability=True,
            random_state=int(rng.integers(0, 2**31 - 1)),
        ),
    )


def _normalize_probabilities(values: np.ndarray) -> np.ndarray | None:
    total = float(np.nansum(values))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return values / total


def _coerce_sampling_trace(
    sampling_trace: PredictionSamplingTrace | str | Path | None,
) -> PredictionSamplingTrace | None:
    if sampling_trace is None:
        return None
    if isinstance(sampling_trace, PredictionSamplingTrace):
        return sampling_trace
    return PredictionSamplingTrace.from_trace_directory(sampling_trace)


def _resolve_sampled_site_positions(
    *,
    available_sites: pd.Index,
    sampled_sites: list[str],
    expected_size: int,
    context: str,
) -> np.ndarray:
    sampled_site_list = _validate_override_sites(
        available_sites=available_sites,
        sampled_sites=sampled_sites,
        expected_size=expected_size,
        context=context,
    )
    position_lookup: dict[str, int] = {}
    for position, site in enumerate(available_sites.astype(str).tolist()):
        position_lookup.setdefault(site, position)
    return np.asarray([position_lookup[site] for site in sampled_site_list], dtype=int)


def _validate_override_sites(
    *,
    available_sites: pd.Index,
    sampled_sites: list[str],
    expected_size: int,
    context: str,
) -> list[str]:
    sampled_site_list = [str(site) for site in sampled_sites]
    if len(sampled_site_list) != expected_size:
        msg = (
            f"Sampling override for {context} has {len(sampled_site_list)} rows; "
            f"expected {expected_size}"
        )
        raise ValueError(msg)
    available_site_set = set(available_sites.astype(str).tolist())
    invalid_sites = [
        site for site in sampled_site_list if site not in available_site_set
    ]
    if invalid_sites:
        unique_invalid_sites = sorted(set(invalid_sites))
        msg = (
            f"Sampling override for {context} contains sites outside the "
            f"available training rows: {', '.join(unique_invalid_sites)}"
        )
        raise ValueError(msg)
    return sampled_site_list


def _require_sklearn() -> tuple[type, type]:
    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
    except ImportError as exc:  # pragma: no cover - environment-dependent
        msg = (
            "KinasePredictor requires scikit-learn. Install the package with "
            "the 'ml' extra, for example: pip install .[ml]"
        )
        raise ImportError(msg) from exc
    return StandardScaler, SVC


def _validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise ValueError(f"{name} must be at least 1")
