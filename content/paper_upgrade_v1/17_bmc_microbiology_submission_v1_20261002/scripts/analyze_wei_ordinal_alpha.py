#!/usr/bin/env python3
"""Prespecified secondary Low–Middle–High trend analysis for alpha diversity.

The source analysis plan encodes Rb1 group as Low=0, Middle=1 and High=2.
This script fits age- and sex-adjusted linear models with HC3 standard errors
for all four alpha-diversity metrics and applies BH correction across the four
ordinal coefficients.  It does not impute continuous conversion measurements.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "source_data" / "figure4" / "Figure4_alpha_by_donor.tsv"
OUTPUT = ROOT / "source_data" / "figure4" / "Figure4_alpha_ordinal_all50.tsv"
METRICS = ["Observed_ASV", "Chao1", "Shannon", "Simpson"]


def main() -> None:
    data = pd.read_csv(INPUT, sep="\t")
    data["rb1_ordinal"] = data["rb1_group"].map({"Low": 0, "Middle": 1, "High": 2})
    if data["rb1_ordinal"].isna().any() or len(data) != 50:
        raise ValueError("Expected 50 donors with Low/Middle/High rb1_group labels")

    rows: list[dict[str, float | int | str]] = []
    for metric in METRICS:
        fit = smf.ols(f"{metric} ~ rb1_ordinal + age + C(gender)", data=data).fit(cov_type="HC3")
        ci = fit.conf_int().loc["rb1_ordinal"]
        rows.append(
            {
                "metric": metric,
                "n": int(fit.nobs),
                "ordinal_coding": "Low=0; Middle=1; High=2",
                "adjustment": "age + sex",
                "variance": "HC3",
                "beta_per_category": float(fit.params["rb1_ordinal"]),
                "robust_se": float(fit.bse["rb1_ordinal"]),
                "ci95_low": float(ci.iloc[0]),
                "ci95_high": float(ci.iloc[1]),
                "p_value": float(fit.pvalues["rb1_ordinal"]),
            }
        )
    out = pd.DataFrame(rows)
    out["q_BH_four_metrics"] = multipletests(out["p_value"], method="fdr_bh")[1]
    out.to_csv(OUTPUT, sep="\t", index=False)
    print(OUTPUT)


if __name__ == "__main__":
    main()
