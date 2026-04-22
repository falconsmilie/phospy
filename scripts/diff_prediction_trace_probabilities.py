#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare R and Python prediction traces at the learner seam by diffing "
            "per-iteration probabilities, decision values, and class-specific resampling weights."
        )
    )
    parser.add_argument(
        "--r-trace-dir",
        default="tests/fixtures/rewrite_parity/r_reference_l6/prediction_trace",
        help="Directory containing R-exported prediction trace CSVs.",
    )
    parser.add_argument(
        "--python-trace-dir",
        default="tests/fixtures/rewrite_parity/python_reference_l6/prediction_trace",
        help="Directory containing Python-exported prediction trace CSVs.",
    )
    parser.add_argument(
        "--kinases",
        default="MAPK9,IRAK1,TBK1,LCK",
        help="Comma-separated kinase list to compare.",
    )
    parser.add_argument(
        "--top-rows",
        type=int,
        default=10,
        help="Number of largest-difference rows to print for each comparison.",
    )
    return parser.parse_args()


def parse_csv_values(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def read_trace_table(trace_dir: Path, filename: str) -> pd.DataFrame:
    path = trace_dir / filename
    if not path.exists():
        raise FileNotFoundError(f"Trace file not found: {path}")
    return pd.read_csv(path)


def _validate_requested_kinases(
    *,
    frame: pd.DataFrame,
    kinases: list[str],
    source_name: str,
) -> pd.DataFrame:
    available = sorted(frame["kinase"].astype(str).drop_duplicates().tolist())
    requested = [str(kinase) for kinase in kinases]
    missing = [kinase for kinase in requested if kinase not in available]
    if missing:
        msg = (
            f"Requested kinases are missing from {source_name}: {', '.join(missing)}. "
            f"Available kinases: {', '.join(available) or 'none'}"
        )
        raise ValueError(msg)
    return frame.loc[frame["kinase"].astype(str).isin(requested)].copy()


def normalize_probability_table(
    frame: pd.DataFrame, kinases: list[str], source_name: str
) -> pd.DataFrame:
    filtered = _validate_requested_kinases(
        frame=frame, kinases=kinases, source_name=source_name
    )
    filtered["kinase"] = filtered["kinase"].astype(str)
    filtered["ensemble"] = filtered["ensemble"].astype(int)
    filtered["iteration"] = filtered["iteration"].astype(int)
    filtered["site"] = filtered["site"].astype(str)
    filtered["label"] = filtered["label"].astype(str)
    return filtered.sort_values(
        ["kinase", "ensemble", "iteration", "site"], kind="mergesort"
    ).reset_index(drop=True)


def normalize_weight_table(
    frame: pd.DataFrame, kinases: list[str], source_name: str
) -> pd.DataFrame:
    filtered = _validate_requested_kinases(
        frame=frame, kinases=kinases, source_name=source_name
    )
    filtered["kinase"] = filtered["kinase"].astype(str)
    filtered["ensemble"] = filtered["ensemble"].astype(int)
    filtered["iteration"] = filtered["iteration"].astype(int)
    filtered["class_label"] = filtered["class_label"].astype(str)
    filtered["site"] = filtered["site"].astype(str)
    if "normalized_weight" not in filtered.columns and "raw_weight" in filtered.columns:
        totals = filtered.groupby(["kinase", "ensemble", "iteration", "class_label"])[
            "raw_weight"
        ].transform("sum")
        filtered["normalized_weight"] = filtered["raw_weight"] / totals
    return filtered.sort_values(
        ["kinase", "ensemble", "iteration", "class_label", "site"], kind="mergesort"
    ).reset_index(drop=True)


def normalize_probability_parameter_table(
    frame: pd.DataFrame, kinases: list[str], source_name: str
) -> pd.DataFrame:
    filtered = _validate_requested_kinases(
        frame=frame, kinases=kinases, source_name=source_name
    )
    filtered["kinase"] = filtered["kinase"].astype(str)
    filtered["ensemble"] = filtered["ensemble"].astype(int)
    filtered["iteration"] = filtered["iteration"].astype(int)
    filtered["class_pair"] = filtered["class_pair"].astype(str)
    filtered["probA"] = filtered["probA"].astype(float)
    filtered["probB"] = filtered["probB"].astype(float)
    return filtered.sort_values(
        ["kinase", "ensemble", "iteration", "class_pair"], kind="mergesort"
    ).reset_index(drop=True)


def summarize_probability_diff(merged: pd.DataFrame) -> pd.DataFrame:
    summary_rows: list[dict[str, object]] = []
    for kinase, group in merged.groupby("kinase", sort=True):
        diff_1 = (group["prob_class_1_py"] - group["prob_class_1_r"]).abs()
        diff_2 = (group["prob_class_2_py"] - group["prob_class_2_r"]).abs()
        summary_rows.append(
            {
                "kinase": kinase,
                "rows": int(len(group)),
                "prob_class_1_mae": float(diff_1.mean()),
                "prob_class_1_max_abs": float(diff_1.max()),
                "prob_class_2_mae": float(diff_2.mean()),
                "prob_class_2_max_abs": float(diff_2.max()),
            }
        )
    return pd.DataFrame(summary_rows)


def normalize_decision_table(
    frame: pd.DataFrame, kinases: list[str], source_name: str
) -> pd.DataFrame:
    filtered = _validate_requested_kinases(
        frame=frame, kinases=kinases, source_name=source_name
    )
    filtered["kinase"] = filtered["kinase"].astype(str)
    filtered["ensemble"] = filtered["ensemble"].astype(int)
    filtered["iteration"] = filtered["iteration"].astype(int)
    filtered["site"] = filtered["site"].astype(str)
    filtered["label"] = filtered["label"].astype(str)
    return filtered.sort_values(
        ["kinase", "ensemble", "iteration", "site"], kind="mergesort"
    ).reset_index(drop=True)


def summarize_decision_diff(merged: pd.DataFrame) -> pd.DataFrame:
    summary_rows: list[dict[str, object]] = []
    for kinase, group in merged.groupby("kinase", sort=True):
        diff = (
            group["decision_value_class_1_py"] - group["decision_value_class_1_r"]
        ).abs()
        summary_rows.append(
            {
                "kinase": kinase,
                "rows": int(len(group)),
                "decision_value_corr": float(
                    group["decision_value_class_1_py"].corr(
                        group["decision_value_class_1_r"], method="pearson"
                    )
                ),
                "decision_value_mae": float(diff.mean()),
                "decision_value_max_abs": float(diff.max()),
            }
        )
    return pd.DataFrame(summary_rows)


def summarize_weight_diff(merged: pd.DataFrame) -> pd.DataFrame:
    summary_rows: list[dict[str, object]] = []
    for (kinase, class_label), group in merged.groupby(
        ["kinase", "class_label"], sort=True
    ):
        diff = (group["normalized_weight_py"] - group["normalized_weight_r"]).abs()
        summary_rows.append(
            {
                "kinase": kinase,
                "class_label": class_label,
                "rows": int(len(group)),
                "weight_mae": float(diff.mean()),
                "weight_max_abs": float(diff.max()),
            }
        )
    return pd.DataFrame(summary_rows)


def summarize_probability_parameter_diff(merged: pd.DataFrame) -> pd.DataFrame:
    summary_rows: list[dict[str, object]] = []
    for kinase, group in merged.groupby("kinase", sort=True):
        diff_a = (group["probA_py"] - group["probA_r"]).abs()
        diff_b = (group["probB_py"] - group["probB_r"]).abs()
        summary_rows.append(
            {
                "kinase": kinase,
                "rows": int(len(group)),
                "probA_mae": float(diff_a.mean()),
                "probA_max_abs": float(diff_a.max()),
                "probB_mae": float(diff_b.mean()),
                "probB_max_abs": float(diff_b.max()),
            }
        )
    return pd.DataFrame(summary_rows)


def main() -> None:
    args = parse_args()
    kinases = parse_csv_values(args.kinases)
    r_trace_dir = Path(args.r_trace_dir)
    python_trace_dir = Path(args.python_trace_dir)

    r_prob = normalize_probability_table(
        read_trace_table(r_trace_dir, "trace_iteration_probabilities.csv"),
        kinases,
        source_name=f"R trace directory ({r_trace_dir})",
    )
    py_prob = normalize_probability_table(
        read_trace_table(python_trace_dir, "trace_iteration_probabilities.csv"),
        kinases,
        source_name=f"Python trace directory ({python_trace_dir})",
    )
    merged_prob = py_prob.merge(
        r_prob,
        on=["kinase", "ensemble", "iteration", "site", "label"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )
    if merged_prob.empty:
        raise ValueError(
            "No overlapping probability rows were found after merging the requested kinases. "
            "Check that both traces were generated from the same kinase list and fixture state."
        )
    prob_summary = summarize_probability_diff(merged_prob)

    r_decision = normalize_decision_table(
        read_trace_table(r_trace_dir, "trace_iteration_decision_values.csv"),
        kinases,
        source_name=f"R trace directory ({r_trace_dir})",
    )
    py_decision = normalize_decision_table(
        read_trace_table(python_trace_dir, "trace_iteration_decision_values.csv"),
        kinases,
        source_name=f"Python trace directory ({python_trace_dir})",
    )
    merged_decision = py_decision.merge(
        r_decision,
        on=["kinase", "ensemble", "iteration", "site", "label"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )
    if merged_decision.empty:
        raise ValueError(
            "No overlapping decision-value rows were found after merging the requested kinases. "
            "Check that both traces were generated from the same kinase list and fixture state."
        )
    decision_summary = summarize_decision_diff(merged_decision)

    r_weights = normalize_weight_table(
        read_trace_table(r_trace_dir, "trace_iteration_resampling_weights.csv"),
        kinases,
        source_name=f"R trace directory ({r_trace_dir})",
    )
    py_weights = normalize_weight_table(
        read_trace_table(python_trace_dir, "trace_iteration_resampling_weights.csv"),
        kinases,
        source_name=f"Python trace directory ({python_trace_dir})",
    )
    merged_weights = py_weights.merge(
        r_weights,
        on=["kinase", "ensemble", "iteration", "class_label", "site"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )
    if merged_weights.empty:
        raise ValueError(
            "No overlapping resampling-weight rows were found after merging the requested kinases. "
            "Check that both traces were generated from the same kinase list and fixture state."
        )
    weight_summary = summarize_weight_diff(merged_weights)

    r_probability_parameters = normalize_probability_parameter_table(
        read_trace_table(r_trace_dir, "trace_iteration_probability_parameters.csv"),
        kinases,
        source_name=f"R trace directory ({r_trace_dir})",
    )
    py_probability_parameters = normalize_probability_parameter_table(
        read_trace_table(
            python_trace_dir, "trace_iteration_probability_parameters.csv"
        ),
        kinases,
        source_name=f"Python trace directory ({python_trace_dir})",
    )
    merged_probability_parameters = py_probability_parameters.merge(
        r_probability_parameters,
        on=["kinase", "ensemble", "iteration", "class_pair"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )
    if merged_probability_parameters.empty:
        raise ValueError(
            "No overlapping probability-parameter rows were found after merging the requested kinases. "
            "Check that both traces were generated from the same kinase list and fixture state."
        )
    probability_parameter_summary = summarize_probability_parameter_diff(
        merged_probability_parameters
    )

    print("Probability diff summary")
    print(prob_summary.to_string(index=False))
    print()

    top_prob = merged_prob.assign(
        prob_class_1_abs_diff=(
            merged_prob["prob_class_1_py"] - merged_prob["prob_class_1_r"]
        ).abs(),
        prob_class_2_abs_diff=(
            merged_prob["prob_class_2_py"] - merged_prob["prob_class_2_r"]
        ).abs(),
    ).sort_values(["prob_class_1_abs_diff", "prob_class_2_abs_diff"], ascending=False)
    print(f"Top {args.top_rows} probability deltas")
    print(
        top_prob.loc[
            :,
            [
                "kinase",
                "ensemble",
                "iteration",
                "site",
                "label",
                "prob_class_1_py",
                "prob_class_1_r",
                "prob_class_1_abs_diff",
                "prob_class_2_py",
                "prob_class_2_r",
                "prob_class_2_abs_diff",
            ],
        ]
        .head(args.top_rows)
        .to_string(index=False)
    )
    print()

    print("Decision-value diff summary")
    print(decision_summary.to_string(index=False))
    print()

    top_decisions = merged_decision.assign(
        decision_value_abs_diff=(
            merged_decision["decision_value_class_1_py"]
            - merged_decision["decision_value_class_1_r"]
        ).abs()
    ).sort_values("decision_value_abs_diff", ascending=False)
    print(f"Top {args.top_rows} decision-value deltas")
    print(
        top_decisions.loc[
            :,
            [
                "kinase",
                "ensemble",
                "iteration",
                "site",
                "label",
                "decision_value_class_1_py",
                "decision_value_class_1_r",
                "decision_value_abs_diff",
            ],
        ]
        .head(args.top_rows)
        .to_string(index=False)
    )
    print()

    print("Probability-parameter diff summary")
    print(probability_parameter_summary.to_string(index=False))
    print()

    top_probability_parameters = merged_probability_parameters.assign(
        probA_abs_diff=(
            merged_probability_parameters["probA_py"]
            - merged_probability_parameters["probA_r"]
        ).abs(),
        probB_abs_diff=(
            merged_probability_parameters["probB_py"]
            - merged_probability_parameters["probB_r"]
        ).abs(),
    ).sort_values(["probA_abs_diff", "probB_abs_diff"], ascending=False)
    print(f"Top {args.top_rows} probability-parameter deltas")
    print(
        top_probability_parameters.loc[
            :,
            [
                "kinase",
                "ensemble",
                "iteration",
                "class_pair",
                "probA_py",
                "probA_r",
                "probA_abs_diff",
                "probB_py",
                "probB_r",
                "probB_abs_diff",
            ],
        ]
        .head(args.top_rows)
        .to_string(index=False)
    )
    print()

    print("Resampling-weight diff summary")
    print(weight_summary.to_string(index=False))
    print()

    top_weights = merged_weights.assign(
        weight_abs_diff=(
            merged_weights["normalized_weight_py"]
            - merged_weights["normalized_weight_r"]
        ).abs()
    ).sort_values("weight_abs_diff", ascending=False)
    print(f"Top {args.top_rows} resampling-weight deltas")
    print(
        top_weights.loc[
            :,
            [
                "kinase",
                "ensemble",
                "iteration",
                "class_label",
                "site",
                "normalized_weight_py",
                "normalized_weight_r",
                "weight_abs_diff",
            ],
        ]
        .head(args.top_rows)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
