#!/usr/bin/env python3
"""Independent revalidation of the frozen B47 HbA1c/FPG endpoint data.

Produces ordinary and modified Hartung–Knapp results, multiplicity-adjusted
inference, study weights, prediction intervals and leave-one-out analyses.
It does not alter the extraction master or the previous meta-analysis release.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import optimize, stats


def bh(p: pd.Series) -> np.ndarray:
    x = pd.to_numeric(p, errors="coerce").to_numpy(float)
    ans = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    order = np.argsort(x[ok])
    ranked = x[ok][order]
    q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.clip(q, 0, 1)
    ans[np.flatnonzero(ok)] = restored
    return ans


def reml(y: np.ndarray, v: np.ndarray) -> float:
    if len(y) < 2:
        return 0.0

    def nll(tau2: float) -> float:
        w = 1.0 / (v + max(float(tau2), 0.0))
        mu = float(np.sum(w * y) / np.sum(w))
        return 0.5 * (np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))

    upper = max(1.0, float(np.var(y, ddof=1) * 10 + np.max(v)))
    res = optimize.minimize_scalar(nll, bounds=(0, upper), method="bounded", options={"xatol": 1e-14})
    return 0.0 if res.x < 1e-8 else float(res.x)


def pool(y: np.ndarray, v: np.ndarray, modified: bool) -> dict[str, float | int | str | bool]:
    k = len(y)
    if k == 1:
        se = math.sqrt(v[0])
        return {
            "k": 1, "estimate": y[0], "se": se, "ci_lb": y[0] - 1.96 * se,
            "ci_ub": y[0] + 1.96 * se, "p_value": 2 * stats.norm.sf(abs(y[0] / se)),
            "tau2": 0.0, "I2": 0.0, "Q": 0.0, "Q_p": np.nan, "q_scale": np.nan,
            "test": "single_study", "exploratory": False, "prediction_lb": np.nan,
            "prediction_ub": np.nan,
        }
    tau2 = reml(y, v)
    w = 1.0 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    q_scale = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    used_scale = max(1.0, q_scale) if modified else q_scale
    se = math.sqrt(used_scale / np.sum(w))
    crit = stats.t.ppf(0.975, k - 1)
    p = float(2 * stats.t.sf(abs(mu / se), k - 1)) if se > 0 else np.nan
    wf = 1.0 / v
    muf = float(np.sum(wf * y) / np.sum(wf))
    Q = float(np.sum(wf * (y - muf) ** 2))
    I2 = max(0.0, 100 * (Q - (k - 1)) / Q) if Q > 0 else 0.0
    pred_lb = pred_ub = np.nan
    if k >= 3:
        pred_crit = stats.t.ppf(0.975, k - 2)
        pred_se = math.sqrt(tau2 + se**2)
        pred_lb, pred_ub = mu - pred_crit * pred_se, mu + pred_crit * pred_se
    return {
        "k": k, "estimate": mu, "se": se, "ci_lb": mu - crit * se,
        "ci_ub": mu + crit * se, "p_value": p, "tau2": tau2, "I2": I2,
        "Q": Q, "Q_p": float(stats.chi2.sf(Q, k - 1)), "q_scale": q_scale,
        "test": "modified_HK" if modified else "ordinary_HK", "exploratory": k == 2,
        "prediction_lb": pred_lb, "prediction_ub": pred_ub,
    }


def input_audit(d: pd.DataFrame) -> pd.DataFrame:
    recalc_y = pd.to_numeric(d["mean_trt"], errors="coerce") - pd.to_numeric(d["mean_ctrl"], errors="coerce")
    recalc_v = (
        pd.to_numeric(d["sd_trt"], errors="coerce") ** 2 / pd.to_numeric(d["n_trt"], errors="coerce")
        + pd.to_numeric(d["sd_ctrl"], errors="coerce") ** 2 / pd.to_numeric(d["n_ctrl"], errors="coerce")
    )
    return pd.DataFrame({
        "trial_id": d["trial_id"], "herb": d["matched_b47_herb"], "outcome": d["outcome"],
        "yi_stored": d["yi_md"], "yi_recalculated": recalc_y,
        "abs_delta_yi": np.abs(pd.to_numeric(d["yi_md"]) - recalc_y),
        "vi_stored": d["vi"], "vi_recalculated": recalc_v,
        "abs_delta_vi": np.abs(pd.to_numeric(d["vi"]) - recalc_v),
        "valid_n": (pd.to_numeric(d["n_trt"]) > 1) & (pd.to_numeric(d["n_ctrl"]) > 1),
        "positive_sd": (pd.to_numeric(d["sd_trt"]) > 0) & (pd.to_numeric(d["sd_ctrl"]) > 0),
    })


def trial_characteristics(d: pd.DataFrame, master_path: Path) -> pd.DataFrame:
    master = pd.read_csv(master_path, encoding="utf-8-sig")
    ids = set(d["trial_id"].astype(str))
    m = master[master["trial_id"].astype(str).isin(ids)].copy()
    aggregate = d.groupby("trial_id", as_index=False).agg(
        herb=("matched_b47_herb", "first"), outcomes=("outcome", lambda x: "|".join(sorted(set(x)))),
        analysis_n=("n_trt", "max"), analysis_n_control=("n_ctrl", "max"),
        time_point_weeks=("time_point_weeks", "max"), extracted_rows=("outcome", "size"),
    )
    cols = ["trial_id", "year", "title", "doi", "matched_b47_herb", "intervention_detail", "plant_part",
            "dose", "duration_weeks", "design", "country", "disease", "n_total", "n_trt", "n_ctrl",
            "comparator", "cointervention", "pdf_file", "notes"]
    return aggregate.merge(m[cols], on="trial_id", how="left", suffixes=("_meta", "_master"), validate="one_to_one")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rebuild the frozen RCT meta-analysis from structured endpoint data."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to meta_input_clean.csv."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory for regenerated RCT meta-analysis outputs."
    )
    parser.add_argument(
        "--study-master",
        type=Path,
        default=None,
        help="Optional study_master_b47.csv used only to regenerate trial characteristics."
    )
    parser.add_argument(
        "--old-results",
        type=Path,
        default=None,
        help="Optional previous meta_results_summary.csv used only for implementation comparison."
    )
    args = parser.parse_args()

    input_path = args.input
    out = args.output
    study_master = args.study_master
    old_results = args.old_results

    if not input_path.is_file():
        parser.error(f"Input file not found: {input_path}")
    if study_master is not None and not study_master.is_file():
        parser.error(f"Study master not found: {study_master}")
    if old_results is not None and not old_results.is_file():
        parser.error(f"Old results file not found: {old_results}")

    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(input_path)
    required = {"trial_id", "matched_b47_herb", "outcome", "yi_md", "vi", "n_trt", "n_ctrl"}
    if missing := required - set(d.columns):
        raise RuntimeError(f"Missing input columns: {sorted(missing)}")
    if d["trial_id"].nunique() != 19 or len(d) != 33:
        raise RuntimeError(f"Frozen input mismatch: expected 19 trials/33 rows, got {d.trial_id.nunique()}/{len(d)}")

    audit = input_audit(d)
    audit.to_csv(out / "meta_input_arithmetic_audit.tsv", sep="\t", index=False)
    if audit["abs_delta_yi"].max() > 1e-10 or audit["abs_delta_vi"].max() > 1e-8:
        raise RuntimeError("Stored effect arithmetic did not reproduce")
    if study_master is not None:
        trial_characteristics(d, study_master).to_csv(
            out / "trial_characteristics_reconstructed.tsv",
            sep="\t",
            index=False,
        )

    rows, weights, loo_rows = [], [], []
    for (herb, outcome), g in d.groupby(["matched_b47_herb", "outcome"], sort=True):
        y, v = g["yi_md"].to_numpy(float), g["vi"].to_numpy(float)
        for modified in (False, True):
            z = pool(y, v, modified)
            rows.append({"herb": herb, "outcome": outcome, "n_total": int((g.n_trt + g.n_ctrl).sum()),
                         "trial_ids": "|".join(g.trial_id.astype(str)), **z})
        tau2 = reml(y, v) if len(y) > 1 else 0.0
        w = 1 / (v + tau2)
        for tid, yi, vi, wi in zip(g.trial_id, y, v, w):
            weights.append({"herb": herb, "outcome": outcome, "trial_id": tid, "yi": yi,
                            "se": math.sqrt(vi), "random_weight_fraction": wi / w.sum(), "tau2_pool": tau2})
        if len(g) >= 3:
            for tid in g.trial_id:
                keep = g.trial_id != tid
                z = pool(y[keep], v[keep], True)
                loo_rows.append({"herb": herb, "outcome": outcome, "omitted_trial": tid, **z})

    res = pd.DataFrame(rows)
    pooled = (res["k"] >= 2) & (res["test"] == "modified_HK")
    res["bh_q_across_herb_outcome_pools"] = np.nan
    res.loc[pooled, "bh_q_across_herb_outcome_pools"] = bh(res.loc[pooled, "p_value"])
    res["publication_bias_assessment"] = np.where(res["k"] >= 10, "eligible", "not_assessed_k_lt_10")
    res.to_csv(out / "meta_results_revalidated.tsv", sep="\t", index=False)
    pd.DataFrame(weights).to_csv(out / "study_random_effect_weights.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows).to_csv(out / "leave_one_out.tsv", sep="\t", index=False)

    if old_results is not None:
        old = pd.read_csv(old_results)
        current = res[res.test == "ordinary_HK"].rename(
            columns={"herb": "matched_b47_herb"}
        )
        cmp = current.merge(
            old,
            on=["matched_b47_herb", "outcome"],
            suffixes=("_rebuilt", "_old"),
            validate="one_to_one",
        )
        for col in ["estimate", "ci_lb", "ci_ub", "p_value", "I2", "tau2", "Q"]:
            cmp[f"delta_{col}"] = (
                pd.to_numeric(cmp[f"{col}_rebuilt"])
                - pd.to_numeric(cmp[f"{col}_old"])
            )
        cmp.to_csv(
            out / "comparison_with_previous_meta.tsv",
            sep="\t",
            index=False,
        )

    spec = {
        "input": input_path.name, "frozen_trial_count": 19, "frozen_endpoint_rows": 33,
        "effect": "raw mean difference; intervention minus control",
        "units": {"HbA1c": "%", "FPG": "mg/dL"},
        "tau2": "REML", "primary_small_sample_inference": "modified Hartung-Knapp, scale=max(1,q)",
        "secondary_inference": "ordinary Hartung-Knapp", "multiplicity": "BH across all herb-by-outcome pools with k>=2",
        "k2": "exploratory", "k1": "narrative only", "publication_bias": "not assessed where k<10",
        "prediction_interval": "t(k-2), reported only for k>=3", "python": platform.python_version(),
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
    }
    (out / "analysis_spec.json").write_text(json.dumps(spec, ensure_ascii=False, indent=2), encoding="utf-8")
    primary = res[(res.test == "modified_HK") & (res.k >= 2)].sort_values("p_value")
    print(primary[["herb", "outcome", "k", "estimate", "ci_lb", "ci_ub", "p_value", "bh_q_across_herb_outcome_pools", "I2"]].to_string(index=False))


if __name__ == "__main__":
    main()
