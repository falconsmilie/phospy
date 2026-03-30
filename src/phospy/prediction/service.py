from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..scoring import KinaseScoringResult
from ..types import PredictionSvmMode, PredictionTraceFormat, PredictionTraceLevel
from ..validation.errors import (
    InputCompatibilityError,
    PredictionConfigurationError,
)
from .models import KinasePredictionDebugTrace, KinasePredictionResult
from .sampling import (
    coerce_sampling_trace,
    make_prediction_random_generators,
    multi_ada_sampling,
    validate_override_sites,
)
from .traces import PredictionSamplingTrace, TraceSink, create_trace_sink
from .validation import (
    validate_positive_int,
    validate_svm_mode,
    validate_trace_format,
    validate_trace_level,
)


class KinasePredictor:
    """Predict kinase-substrate relationships from phosphosite score matrices.

    This is a narrow native Python port of PhosR's score-to-prediction seam.
    It mirrors the broad structure of ``kinaseSubstratePred()``: candidate
    substrates are selected from the combined score matrix, then an ensemble of
    adaptive SVM models is used to produce a kinase prediction matrix.

    Use ``svm_mode='r_parity'`` when you want settings that more closely match
    PhosR's e1071-based learner seam. The default mode preserves the standard
    scikit-learn behaviour.
    """

    def __init__(
        self,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.kernel = kernel
        self.svm_mode = validate_svm_mode(svm_mode)

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
        svm_mode: PredictionSvmMode | None = None,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
        trace_level: PredictionTraceLevel | None = None,
        trace_sink: TraceSink | str | Path | None = None,
        trace_sink_format: PredictionTraceFormat = "csv",
    ) -> KinasePredictionResult:
        validate_positive_int(ensemble_size, name="ensemble_size")
        validate_positive_int(top, name="top")
        validate_positive_int(inclusion, name="inclusion")
        validate_positive_int(n_iterations, name="n_iterations")
        validate_positive_int(debug_top_n, name="debug_top_n")
        resolved_svm_mode = (
            self.svm_mode if svm_mode is None else validate_svm_mode(svm_mode)
        )
        resolved_trace_level = validate_trace_level(
            "summary"
            if trace_level is None and capture_debug_trace
            else trace_level or "none"
        )
        resolved_trace_format = validate_trace_format(trace_sink_format)
        sampling_trace_obj = coerce_sampling_trace(sampling_trace)
        trace_sink_obj = (
            create_trace_sink(trace_sink, fmt=resolved_trace_format)
            if resolved_trace_level == "full"
            else None
        )

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
                debug_traces={} if resolved_trace_level != "none" else None,
                trace_level=resolved_trace_level,
                trace_sink=trace_sink_obj,
            )

        master_rng = np.random.default_rng(random_state)
        feature_mat = combined_scores.astype(float)
        pred_matrix = pd.DataFrame(
            0.0,
            index=feature_mat.index.copy(),
            columns=list(substrate_list),
        )

        traced_kinases = (
            set(substrate_list)
            if resolved_trace_level != "none" and debug_kinases is None
            else set(debug_kinases or [])
        )
        debug_traces: dict[str, KinasePredictionDebugTrace] | None = (
            {} if resolved_trace_level != "none" else None
        )

        for kinase, substrates in substrate_list.items():
            negative_sampling_rng, resampling_rng = make_prediction_random_generators(
                master_rng
            )
            positive_train = feature_mat.loc[substrates, :]
            negative_pool = feature_mat.loc[
                ~feature_mat.index.isin(substrates),
                :,  # noqa: E203
            ]
            if negative_pool.empty:
                msg = (
                    f"No negative pool sites available to train predictor for {kinase}"
                )
                raise PredictionConfigurationError(msg)

            is_traced_kinase = (
                resolved_trace_level != "none" and kinase in traced_kinases
            )
            if is_traced_kinase and debug_traces is not None:
                debug_traces[kinase] = KinasePredictionDebugTrace(
                    kinase=kinase,
                    candidate_substrates=list(substrates),
                    negative_pool_sites=negative_pool.index.tolist(),
                    ensemble_traces=[],
                )
                if trace_sink_obj is not None:
                    trace_sink_obj.write_rows(
                        "trace_selected_candidates",
                        [
                            {
                                "kinase": kinase,
                                "candidate_rank": rank,
                                "site": site,
                            }
                            for rank, site in enumerate(substrates, start=1)
                        ],
                    )
                    trace_sink_obj.write_rows(
                        "trace_negative_pool",
                        [
                            {
                                "kinase": kinase,
                                "pool_index": pool_index,
                                "site": site,
                            }
                            for pool_index, site in enumerate(
                                negative_pool.index.tolist(), start=1
                            )
                        ],
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
                    negative_sites = validate_override_sites(
                        available_sites=negative_pool.index,
                        sampled_sites=ensemble_override.initial_negative_sites,
                        expected_size=len(positive_train),
                        context=(
                            f"initial negatives for kinase={kinase}, "
                            f"ensemble={ensemble_index}"
                        ),
                    )
                else:
                    negative_indices = negative_sampling_rng.choice(
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

                if is_traced_kinase and trace_sink_obj is not None:
                    trace_sink_obj.write_rows(
                        "trace_initial_negatives",
                        [
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_index,
                                "draw": draw,
                                "site": site,
                            }
                            for draw, site in enumerate(negative_sites, start=1)
                        ],
                    )

                if is_traced_kinase and debug_traces is not None:
                    series, ensemble_trace = multi_ada_sampling(
                        train_mat=train_mat,
                        test_mat=feature_mat,
                        labels=labels,
                        kernel=self.kernel,
                        n_iterations=n_iterations,
                        resampling_rng=resampling_rng,
                        capture_trace=True,
                        trace_level=resolved_trace_level,
                        trace_sink=trace_sink_obj,
                        kinase=kinase,
                        ensemble_index=ensemble_index,
                        initial_negative_sites=negative_sites,
                        debug_top_n=debug_top_n,
                        svm_mode=resolved_svm_mode,
                        sampling_override=ensemble_override,
                    )
                    if ensemble_trace is not None:
                        debug_traces[kinase].ensemble_traces.append(ensemble_trace)
                else:
                    series, _ = multi_ada_sampling(
                        train_mat=train_mat,
                        test_mat=feature_mat,
                        labels=labels,
                        kernel=self.kernel,
                        n_iterations=n_iterations,
                        resampling_rng=resampling_rng,
                        capture_trace=False,
                        trace_level="none",
                        trace_sink=None,
                        kinase=kinase,
                        ensemble_index=ensemble_index,
                        initial_negative_sites=negative_sites,
                        debug_top_n=debug_top_n,
                        svm_mode=resolved_svm_mode,
                        sampling_override=ensemble_override,
                    )

                pred_matrix.loc[:, kinase] += series

        pred_matrix /= float(ensemble_size)
        return KinasePredictionResult(
            pred_matrix=pred_matrix,
            substrate_list=substrate_list,
            debug_traces=debug_traces,
            trace_level=resolved_trace_level,
            trace_sink=trace_sink_obj,
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
        svm_mode: PredictionSvmMode | None = None,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
        trace_level: PredictionTraceLevel | None = None,
        trace_sink: TraceSink | str | Path | None = None,
        trace_sink_format: PredictionTraceFormat = "csv",
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
            raise InputCompatibilityError(msg)

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
            svm_mode=svm_mode,
            sampling_trace=sampling_trace,
            trace_level=trace_level,
            trace_sink=trace_sink,
            trace_sink_format=trace_sink_format,
        )


def build_candidate_substrate_list(
    combined_scores: pd.DataFrame,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
) -> dict[str, list[str]]:
    """Select candidate kinase substrates from the combined score matrix."""

    validate_positive_int(top, name="top")
    validate_positive_int(inclusion, name="inclusion")

    substrate_list: dict[str, list[str]] = {}
    for kinase in combined_scores.columns:
        selected = combined_scores.loc[:, kinase].nlargest(top)
        sites = selected.loc[selected > score_threshold].index.tolist()
        if len(sites) >= inclusion:
            substrate_list[kinase] = sites
    return substrate_list


__all__ = ["KinasePredictor", "build_candidate_substrate_list"]
