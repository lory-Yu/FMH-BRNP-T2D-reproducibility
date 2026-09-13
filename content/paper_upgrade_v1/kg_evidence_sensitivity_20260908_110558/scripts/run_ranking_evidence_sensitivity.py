#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_KG_DIR = SCRIPT_DIR.parent
DEFAULT_V2_DIR = DEFAULT_KG_DIR.parent / "metabolism_aware_v2_20260907_013732"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    sd = np.std(x, ddof=1)
    return np.zeros(len(x)) if sd <= 0 else (x - np.mean(x)) / sd


def rank_desc(x: np.ndarray) -> np.ndarray:
    return pd.Series(x).rank(method="min", ascending=False).astype(int).to_numpy()


def aggregate(paths: pd.DataFrame, herbs: pd.DataFrame, scenario: str, deduplicate_structure: bool = False) -> tuple[pd.DataFrame, dict]:
    p = paths.copy()
    if deduplicate_structure and len(p):
        p["chem_route_key"] = p[["herb_canonical", "parent_final_parent_inchikey", "product_final_parent_inchikey", "target_gene", "microbe_name", "route_type"]].fillna("").agg("|".join, axis=1)
        p = p.sort_values(["primary_path_weight", "route_id"], ascending=[False, True]).drop_duplicates("chem_route_key")
    agg = p.groupby("herb_canonical").agg(
        transform_score_raw=("primary_path_weight", "sum"),
        n_paths=("route_id", "size"), n_routes=("route_id", "nunique"),
        n_route_substrates=("parent_cid", "nunique"), n_products=("product_cid", "nunique"),
        n_targets=("target_gene", "nunique"),
    ).reset_index()
    out = herbs.merge(agg, on="herb_canonical", how="left")
    for col in ["transform_score_raw", "n_paths", "n_routes", "n_route_substrates", "n_products", "n_targets"]:
        out[col] = out[col].fillna(0)
    y = np.log1p(out["transform_score_raw"].to_numpy(float))
    X = np.column_stack([np.ones(len(out)), np.log1p(out["n_compounds"].to_numpy(float)), np.log1p(out["n_route_substrates"].to_numpy(float))])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    out["transform_adjusted_z"] = zscore(resid)
    out["transform_rank"] = rank_desc(out["transform_adjusted_z"].to_numpy())
    out["scenario"] = scenario
    reg = {"scenario": scenario, "intercept": float(beta[0]), "beta_log1p_n_compounds": float(beta[1]), "beta_log1p_n_route_substrates": float(beta[2]), "residual_sd": float(np.std(resid, ddof=1)), "n_paths_used": int(len(p)), "n_routes_used": int(p["route_id"].nunique())}
    return out, reg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reproduce ranking evidence-sensitivity scenarios.")
    parser.add_argument("--v2-dir", type=Path, default=DEFAULT_V2_DIR)
    parser.add_argument("--qc-dir", type=Path, default=DEFAULT_KG_DIR / "01_chemistry")
    parser.add_argument("--output", type=Path, default=DEFAULT_KG_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    v2_dir = args.v2_dir.resolve()
    qc_dir = args.qc_dir.resolve()
    out = args.output.resolve()
    ranking_out = out / "02_ranking"
    logs_out = out / "logs"
    ranking_out.mkdir(parents=True, exist_ok=True)
    logs_out.mkdir(parents=True, exist_ok=True)

    paths_path = v2_dir / "01_routes/herb_transformation_paths_v2.tsv"
    ranks_path = v2_dir / "02_rankings/KG_metabolism_aware_v2_two_axis.tsv"
    local_route_qc_path = qc_dir / "07_route_entity_crosswalk.tsv"
    final_route_qc_path = qc_dir / "11_route_entity_crosswalk_final.tsv"

    paths = pd.read_csv(paths_path, sep="\t", dtype={"parent_cid": str, "product_cid": str}, low_memory=False)
    deposited = pd.read_csv(ranks_path, sep="\t", low_memory=False)
    herbs = deposited[["herb_canonical", "n_compounds", "primary_transform_rank", "primary_transform_score_adjusted_z"]].copy()
    local = pd.read_csv(local_route_qc_path, sep="\t", dtype={"parent_cid": str, "product_cid": str}, low_memory=False)[["route_id", "both_entities_strict_usable"]]
    final = pd.read_csv(final_route_qc_path, sep="\t", dtype={"parent_cid": str, "product_cid": str}, low_memory=False)[["route_id", "both_entities_final_strict_usable", "parent_final_parent_inchikey", "product_final_parent_inchikey"]]
    p = paths.merge(local, on="route_id", how="left", validate="many_to_one").merge(final, on="route_id", how="left", validate="many_to_one")
    base_eligible = p["primary_path_eligible"].fillna(False).astype(bool)
    gins_text = p[["substrate_name", "product_name"]].fillna("").agg(" ".join, axis=1).str.contains(r"ginsenoside|compound\s*k|protopanaxadiol|\brb1\b|\brd\b|\bf2\b", case=False, regex=True)

    scenarios = [
        ("original_primary", base_eligible, False),
        ("local_consensus_strict", base_eligible & p["both_entities_strict_usable"].fillna(False).astype(bool), False),
        ("pubchem_resolved_strict", base_eligible & p["both_entities_final_strict_usable"].fillna(False).astype(bool), False),
        ("pubchem_resolved_structure_deduplicated", base_eligible & p["both_entities_final_strict_usable"].fillna(False).astype(bool), True),
        ("bbBgl_dispute_removed", base_eligible & ~gins_text, False),
    ]
    outputs, regressions = [], []
    for name, mask, dedup in scenarios:
        result, reg = aggregate(p[mask].copy(), herbs[["herb_canonical", "n_compounds"]], name, dedup)
        outputs.append(result)
        regressions.append(reg)
    long = pd.concat(outputs, ignore_index=True)
    long.to_csv(ranking_out / "01_ranking_scenarios_long.tsv", sep="\t", index=False)

    original = long[long["scenario"] == "original_primary"].set_index("herb_canonical")
    deposited_idx = deposited.set_index("herb_canonical")
    reproduces = original["transform_rank"].astype(int).equals(deposited_idx.loc[original.index, "primary_transform_rank"].astype(int))
    top_original = set(original.nsmallest(20, "transform_rank").index)
    summary_rows = []
    shift_rows = []
    for scenario, g in long.groupby("scenario", sort=False):
        x = g.set_index("herb_canonical").loc[original.index]
        rho = spearmanr(original["transform_rank"], x["transform_rank"]).statistic
        top = set(x.nsmallest(20, "transform_rank").index)
        summary_rows.append({
            "scenario": scenario, "n_paths_used": int(x["n_paths"].sum()), "n_routes_used": int(x["n_routes"].max()) if False else int(p.loc[dict((n,m) for n,m,_ in scenarios)[scenario], "route_id"].nunique()),
            "n_herbs_with_paths": int(x["n_paths"].gt(0).sum()), "spearman_rank_vs_original": float(rho),
            "top20_overlap_vs_original": len(top & top_original), "top20_jaccard_vs_original": len(top & top_original) / len(top | top_original),
            "max_absolute_rank_shift": int((x["transform_rank"] - original["transform_rank"]).abs().max()),
            "bbBgl_or_ginsenoside_paths_removed": int((base_eligible & gins_text).sum()) if scenario == "bbBgl_dispute_removed" else 0,
        })
        delta = (x["transform_rank"] - original["transform_rank"]).rename("rank_delta_vs_original")
        temp = x[["transform_rank", "transform_score_raw", "n_paths", "n_routes"]].join(original["transform_rank"].rename("original_rank")).join(delta)
        temp["scenario"] = scenario
        shift_rows.append(temp.reset_index())
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(ranking_out / "02_ranking_scenario_summary.tsv", sep="\t", index=False)
    shifts = pd.concat(shift_rows, ignore_index=True)
    shifts.sort_values(["scenario", "rank_delta_vs_original"], key=lambda s: s.abs() if s.name == "rank_delta_vs_original" else s, ascending=[True, False]).to_csv(ranking_out / "03_herb_rank_shifts.tsv", sep="\t", index=False)

    top_table = long[long["transform_rank"] <= 20].sort_values(["scenario", "transform_rank"])
    top_table.to_csv(ranking_out / "04_top20_by_scenario.tsv", sep="\t", index=False)
    pd.DataFrame(regressions).to_csv(ranking_out / "05_regression_diagnostics.tsv", sep="\t", index=False)

    report = {
        "original_primary_reproduces_deposited_v2_rank": bool(reproduces),
        "v2_path_rows": int(len(paths)),
        "primary_path_rows": int(base_eligible.sum()),
        "ginsenoside_or_BbBgl_path_rows_in_v2": int((base_eligible & gins_text).sum()),
        "new_BglX_Rb1_to_Rd_edge_enters_existing_v2_rank": False,
        "reason_new_BglX_does_not_enter_rank": "The v2 transformation route table contains no Rb1, Rd, F2 or Compound K route and no independently curated Rd-to-T2D-target terminal route.",
        "scenario_summary": summary.to_dict(orient="records"),
        "input_sha256": {
            "herb_transformation_paths_v2.tsv": sha256(paths_path),
            "KG_metabolism_aware_v2_two_axis.tsv": sha256(ranks_path),
            "07_route_entity_crosswalk.tsv": sha256(local_route_qc_path),
            "11_route_entity_crosswalk_final.tsv": sha256(final_route_qc_path),
        },
    }
    (logs_out / "ranking_sensitivity_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
