#!/usr/bin/env python3
"""Post-run invariant tests for metabolism-aware v2 outputs."""
from pathlib import Path
import json
import sys

import pandas as pd


OUT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    routes = pd.read_csv(OUT / "01_routes/transformation_routes_audited_v2.tsv", sep="\t")
    paths = pd.read_csv(OUT / "01_routes/herb_transformation_paths_v2.tsv", sep="\t")
    ranking = pd.read_csv(OUT / "02_rankings/KG_metabolism_aware_v2_two_axis.tsv", sep="\t")
    alpha = pd.read_csv(OUT / "02_rankings/KG_metabolism_aware_v2_alpha_sensitivity.tsv", sep="\t")
    checks = json.loads((OUT / "03_validation/validation_checks.json").read_text(encoding="utf-8"))

    require(len(routes) == 95, "expected 95 audited terminal routes")
    require(routes["route_id"].nunique() == 95, "route IDs must be unique")
    require((routes["route_type"] == "microbial").sum() == 75, "expected 75 microbial routes")
    require((routes["route_type"] == "host").sum() == 20, "expected 20 host routes")
    require(routes["mapping_conflict"].sum() == 8, "expected 8 audited mapping conflicts")
    require(paths.loc[paths["mapping_conflict"], "inclusive_path_weight"].eq(0).all(), "conflicts must have zero weight")
    require(len(ranking) == 104 and ranking["herb_canonical"].nunique() == 104, "all 104 materials must remain")
    require(ranking[["direct_adjusted_z", "primary_transform_score_adjusted_z"]].notna().all().all(), "scores cannot be missing")
    require(alpha["alpha"].nunique() == 5 and len(alpha) == 520, "alpha grid must contain 5 x 104 rows")
    require(bool(checks["frozen_reported_order_matches_frozen_score"]), "frozen baseline order check failed")
    require(bool(checks["alpha_zero_matches_direct_adjusted_comparator"]), "alpha zero reproduction failed")
    require(not bool(checks["clinical_data_used_for_tuning"]), "clinical tuning is prohibited")
    print("PASS: all metabolism-aware v2 output invariants satisfied")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise
