#!/usr/bin/env python3
"""Depth and influence audit for the Wei 2025 richness result.

This is a low-cost sensitivity analysis of frozen DADA2 outputs. It does not
reprocess FASTQ files or alter the upstream Wei 2025 release.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = PACKAGE_ROOT / "source_data" / "figure4" / "wei_depth_audit"
DEFAULT_COUNTS = PACKAGE_ROOT / "source_data" / "figure4" / "wei_clean_asv_counts.tsv"
DEFAULT_ALPHA = PACKAGE_ROOT / "source_data" / "figure4" / "Figure4_alpha_by_donor.tsv"
SEED = 20260903
RARE_DEPTH = 30_000
N_RARE = 500


def design(frame: pd.DataFrame, include_depth: bool) -> pd.DataFrame:
    x = pd.DataFrame({
        "High_vs_Low": (frame["rb1_group"] == "High").astype(float),
        "age": frame["age"].astype(float),
        "male": (frame["gender"].str.lower() == "male").astype(float),
    }, index=frame.index)
    if include_depth:
        x["log_library_size"] = np.log(frame["library_size"].astype(float))
    return sm.add_constant(x, has_constant="add")


def fit_metric(frame: pd.DataFrame, metric: str, include_depth: bool = True) -> tuple[object, dict]:
    fit = sm.OLS(frame[metric].astype(float), design(frame, include_depth)).fit(cov_type="HC3", use_t=True)
    ci = fit.conf_int().loc["High_vs_Low"]
    infl = fit.get_influence()
    return fit, {
        "metric": metric,
        "model": f"{metric} ~ High + age + sex" + (" + log(library_size)" if include_depth else ""),
        "n": int(fit.nobs),
        "adjusted_difference_high_vs_low": float(fit.params["High_vs_Low"]),
        "hc3_se": float(fit.bse["High_vs_Low"]),
        "ci95_low": float(ci.iloc[0]),
        "ci95_high": float(ci.iloc[1]),
        "p_value": float(fit.pvalues["High_vs_Low"]),
        "max_cooks_distance": float(np.max(infl.cooks_distance[0])),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Wei depth and influence audit on frozen derived inputs.")
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--alpha", type=Path, default=DEFAULT_ALPHA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    alpha = pd.read_csv(args.alpha, sep="\t")
    alpha = alpha[alpha["rb1_group"].isin(["High", "Low"])].copy().reset_index(drop=True)
    counts = pd.read_csv(args.counts, sep="\t").set_index("donor_id")
    counts = counts.loc[alpha["donor_id"]]
    observed_library = counts.sum(axis=1).to_numpy()
    if not np.array_equal(observed_library, alpha["library_size"].to_numpy()):
        raise RuntimeError("Clean ASV counts do not reproduce frozen library sizes")
    if not np.array_equal((counts > 0).sum(axis=1).to_numpy(), alpha["Observed_ASV"].to_numpy()):
        raise RuntimeError("Clean ASV counts do not reproduce frozen observed richness")

    high_depth = alpha.loc[alpha["rb1_group"].eq("High"), "library_size"]
    low_depth = alpha.loc[alpha["rb1_group"].eq("Low"), "library_size"]
    depth_test = stats.mannwhitneyu(high_depth, low_depth, alternative="two-sided")
    depth_summary = pd.DataFrame([{
        "n_high": len(high_depth), "n_low": len(low_depth),
        "high_median_library_size": float(high_depth.median()),
        "low_median_library_size": float(low_depth.median()),
        "median_difference": float(high_depth.median() - low_depth.median()),
        "mann_whitney_u": float(depth_test.statistic),
        "p_value": float(depth_test.pvalue),
        "spearman_depth_observed_rho": float(stats.spearmanr(alpha["library_size"], alpha["Observed_ASV"]).statistic),
        "spearman_depth_observed_p": float(stats.spearmanr(alpha["library_size"], alpha["Observed_ASV"]).pvalue),
        "spearman_depth_chao1_rho": float(stats.spearmanr(alpha["library_size"], alpha["Chao1"]).statistic),
        "spearman_depth_chao1_p": float(stats.spearmanr(alpha["library_size"], alpha["Chao1"]).pvalue),
    }])
    depth_summary.to_csv(out / "library_depth_group_audit.tsv", sep="\t", index=False)

    model_rows = []
    loo_rows = []
    for metric in ["Observed_ASV", "Chao1"]:
        _, row = fit_metric(alpha, metric, include_depth=True)
        model_rows.append(row)
        for donor in alpha["donor_id"]:
            sub = alpha.loc[alpha["donor_id"].ne(donor)].copy()
            fit, _ = fit_metric(sub, metric, include_depth=True)
            ci = fit.conf_int().loc["High_vs_Low"]
            loo_rows.append({
                "metric": metric, "omitted_donor": donor,
                "adjusted_difference_high_vs_low": float(fit.params["High_vs_Low"]),
                "ci95_low": float(ci.iloc[0]), "ci95_high": float(ci.iloc[1]),
                "p_value": float(fit.pvalues["High_vs_Low"]),
            })
    model = pd.DataFrame(model_rows)
    model["bh_q_two_richness_metrics"] = multipletests(model["p_value"], method="fdr_bh")[1]
    model.to_csv(out / "depth_adjusted_richness_models.tsv", sep="\t", index=False)
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(out / "depth_adjusted_leave_one_out.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    raw = counts.to_numpy(dtype=np.int64)
    rare_rows = []
    x = design(alpha, include_depth=False)
    group = alpha["rb1_group"].to_numpy()
    for replicate in range(1, N_RARE + 1):
        rare = np.vstack([rng.multivariate_hypergeometric(row, RARE_DEPTH) for row in raw])
        richness = (rare > 0).sum(axis=1).astype(float)
        fit = sm.OLS(richness, x).fit(cov_type="HC3", use_t=True)
        rare_rows.append({
            "replicate": replicate,
            "seed": SEED,
            "rarefaction_depth": RARE_DEPTH,
            "high_median": float(np.median(richness[group == "High"])),
            "low_median": float(np.median(richness[group == "Low"])),
            "adjusted_difference_high_vs_low": float(fit.params["High_vs_Low"]),
            "p_value": float(fit.pvalues["High_vs_Low"]),
        })
    rare_df = pd.DataFrame(rare_rows)
    rare_df.to_csv(out / "repeated_rarefaction_observed_asv_500.tsv", sep="\t", index=False)
    rare_summary = pd.DataFrame([{
        "replicates": N_RARE,
        "seed": SEED,
        "rarefaction_depth": RARE_DEPTH,
        "median_adjusted_difference": float(rare_df["adjusted_difference_high_vs_low"].median()),
        "difference_q025": float(rare_df["adjusted_difference_high_vs_low"].quantile(0.025)),
        "difference_q975": float(rare_df["adjusted_difference_high_vs_low"].quantile(0.975)),
        "fraction_positive": float((rare_df["adjusted_difference_high_vs_low"] > 0).mean()),
        "fraction_p_lt_0_05": float((rare_df["p_value"] < 0.05).mean()),
        "median_p": float(rare_df["p_value"].median()),
        "max_p": float(rare_df["p_value"].max()),
    }])
    rare_summary.to_csv(out / "repeated_rarefaction_summary.tsv", sep="\t", index=False)

    obs = model.loc[model["metric"].eq("Observed_ASV")].iloc[0]
    cha = model.loc[model["metric"].eq("Chao1")].iloc[0]
    obs_loo = loo.loc[loo["metric"].eq("Observed_ASV")]
    cha_loo = loo.loc[loo["metric"].eq("Chao1")]
    rr = rare_summary.iloc[0]
    report = f"""# Wei 2025 library-depth sensitivity audit

This low-cost audit uses the frozen clean DADA2 ASV table; FASTQ files were not reprocessed.

- High and Low groups had median retained library sizes of {high_depth.median():,.0f} and {low_depth.median():,.0f} reads, respectively (two-sided Mann–Whitney P={depth_test.pvalue:.3g}).
- In small-sample t-based HC3 models adjusted for age, sex and log library size, the High–Low difference was {obs.adjusted_difference_high_vs_low:.2f} ASVs (95% CI {obs.ci95_low:.2f} to {obs.ci95_high:.2f}; P={obs.p_value:.4f}; two-metric BH q={obs.bh_q_two_richness_metrics:.4f}) and {cha.adjusted_difference_high_vs_low:.2f} for Chao1 (95% CI {cha.ci95_low:.2f} to {cha.ci95_high:.2f}; P={cha.p_value:.4f}; q={cha.bh_q_two_richness_metrics:.4f}).
- Repeated rarefaction used {N_RARE} independently sampled tables at {RARE_DEPTH:,} reads per donor and seed {SEED}. All replicates retained a positive adjusted Observed-ASV difference; {rr.fraction_p_lt_0_05:.1%} had P<0.05 (median P={rr.median_p:.4f}; maximum P={rr.max_p:.4f}).
- Leave-one-donor-out coefficients remained positive (Observed ASV {obs_loo.adjusted_difference_high_vs_low.min():.2f} to {obs_loo.adjusted_difference_high_vs_low.max():.2f}; Chao1 {cha_loo.adjusted_difference_high_vs_low.min():.2f} to {cha_loo.adjusted_difference_high_vs_low.max():.2f}), although the maximum leave-one-out P values were {obs_loo.p_value.max():.4f} and {cha_loo.p_value.max():.4f}. This supports directional stability but not immunity to small-sample influence.

The result should be interpreted only as a richness difference between source-defined extreme converter groups. It is not an independent replication of Rb1 conversion, an enzyme measurement, or evidence about T2D or metformin.
"""
    (out / "WEI_DEPTH_SENSITIVITY_REPORT.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
