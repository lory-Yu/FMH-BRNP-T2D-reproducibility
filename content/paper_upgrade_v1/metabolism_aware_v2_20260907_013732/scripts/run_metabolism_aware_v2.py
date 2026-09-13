#!/usr/bin/env python3
"""Evidence-bounded metabolism-aware FMH prioritisation v2.

Writes only below the supplied output directory. The frozen v1 ranking and all
source data are read-only inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


BASELINE: Path
HERB_COMPOUND: Path
MICRO: Path
MICRO_RAW: Path
HOST: Path
EDGE_AUDIT: Path
ENTITY_AUDIT: Path

ALPHAS = [0.0, 0.25, 0.50, 0.75, 1.0]
TOP_N = 20
N_NULL = 1000
SEED = 20260907


def cid_series(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce").astype("Int64").astype(str).replace("<NA>", "")


def split_tokens(values: pd.Series) -> list[str]:
    out: set[str] = set()
    for value in values.dropna().astype(str):
        for token in value.replace(",", ";").split(";"):
            token = token.strip()
            if token and token.lower() not in {"nan", "none"}:
                out.add(token)
    return sorted(out)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    sd = float(np.std(arr, ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        return np.zeros(len(arr), dtype=float)
    return (arr - float(np.mean(arr))) / sd


def minmax(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    lo, hi = float(np.min(arr)), float(np.max(arr))
    if hi - lo <= 1e-12:
        return np.zeros(len(arr), dtype=float)
    return (arr - lo) / (hi - lo)


def rank_desc(values: pd.Series | np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="min", ascending=False).astype(int).to_numpy()


def degree_adjust(
    score: pd.Series | np.ndarray,
    n_compounds: pd.Series | np.ndarray,
    n_opportunities: pd.Series | np.ndarray,
    degree_2_label: str = "n_opportunities",
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    y = np.asarray(score, dtype=float)
    x = np.column_stack(
        [
            np.ones(len(y)),
            np.log1p(np.asarray(n_compounds, dtype=float)),
            np.log1p(np.asarray(n_opportunities, dtype=float)),
        ]
    )
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    resid = y - fitted
    denom = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - np.sum(resid**2) / denom) if denom > 0 else float("nan")
    return resid, zscore(resid), {
        "intercept": float(beta[0]),
        "beta_log1p_n_compounds": float(beta[1]),
        f"beta_log1p_{degree_2_label}": float(beta[2]),
        "r_squared": r2,
        "residual_sd": float(np.std(resid, ddof=1)),
    }


def build_occurrences() -> tuple[pd.DataFrame, pd.DataFrame]:
    hc = pd.read_csv(HERB_COMPOUND, low_memory=False)
    hc = hc[hc["fmh_official"] == True].copy()  # noqa: E712
    hc["herb_canonical"] = hc["herb_canonical"].astype(str).str.strip()
    hc["parent_cid"] = cid_series(hc["pubchem_cid"])
    hc = hc[hc["parent_cid"] != ""].copy()
    hc["lotus_verified"] = hc["lotus_verified"].fillna(False).astype(bool)
    hc["source_db"] = hc["source_db"].fillna("").astype(str)
    hc["sources_all"] = hc["sources_all"].fillna("").astype(str)
    hc["inchikey_available"] = hc["inchikey"].fillna("").astype(str).str.len() > 0

    rows = []
    for (herb, cid), g in hc.groupby(["herb_canonical", "parent_cid"], sort=False):
        sources = split_tokens(pd.concat([g["source_db"], g["sources_all"]], ignore_index=True))
        rows.append(
            {
                "herb_canonical": herb,
                "parent_cid": cid,
                "compound_names": "; ".join(sorted(set(g["compound_name"].dropna().astype(str)))),
                "lotus_verified_any": bool(g["lotus_verified"].any()),
                "inchikey_available_any": bool(g["inchikey_available"].any()),
                "occurrence_database_count": len(sources),
                "occurrence_databases": "; ".join(sources),
                "occurrence_evidence_class": (
                    "database_cross_verified" if g["lotus_verified"].any()
                    else "database_integrated_unverified"
                ),
                "occurrence_weight": 1.0 if g["lotus_verified"].any() else 0.5,
                "occurrence_primary_publication_retained": False,
            }
        )
    occ = pd.DataFrame(rows)
    herb_degree = (
        hc.groupby("herb_canonical")["dedupe_key"].nunique().rename("n_compounds").reset_index()
    )
    return occ, herb_degree


def build_microbial_routes() -> pd.DataFrame:
    terminal = pd.read_csv(MICRO, low_memory=False)
    raw = pd.read_csv(MICRO_RAW, low_memory=False)
    terminal["parent_cid"] = cid_series(terminal["substrate_pubchem_cid"])
    terminal["product_cid"] = cid_series(terminal["metabolite_pubchem_cid_mg"])

    grouped = []
    for key, g in raw.groupby("dedup_key", sort=False):
        upstream_models = split_tokens(g["human_mouse"])
        target_models = split_tokens(g["human_mouse_mg"])
        trans_pmids = split_tokens(g["pmid"])
        target_pmids = split_tokens(g["pmid_mg"])
        sources = split_tokens(pd.concat([g["source"], g["source_mg"]], ignore_index=True))
        grouped.append(
            {
                "dedup_key": key,
                "upstream_models": ";".join(upstream_models),
                "target_models": ";".join(target_models),
                "transformation_pmids": ";".join(trans_pmids),
                "target_pmids": ";".join(target_pmids),
                "all_source_ids": ";".join(sorted(set(trans_pmids + target_pmids))),
                "source_databases": ";".join(sources),
                "raw_model_rows": len(g),
            }
        )
    support = pd.DataFrame(grouped)
    d = terminal.merge(support, on="dedup_key", how="left", validate="one_to_one")
    d["route_id"] = "MICRO-" + d["link_id"].astype(int).astype(str).str.zfill(3)
    d["route_type"] = "microbial"
    d["product_name"] = d["metabolite_name"]
    d["is_self_loop"] = d["substrate_equals_metabolite"].fillna(False).astype(bool)

    def classify(row: pd.Series) -> tuple[str, float]:
        models = set(str(row["upstream_models"]).split(";")) - {"", "nan"}
        if models == {"human"}:
            return "human_only", 1.0
        if "human" in models and "mouse" in models:
            return "mixed_human_mouse", 0.70
        if models == {"mouse"}:
            return "mouse_only", 0.35
        return "model_unspecified", 0.20

    classes = d.apply(classify, axis=1)
    d["model_evidence_class"] = [x[0] for x in classes]
    d["model_weight"] = [x[1] for x in classes]
    d["association_weight"] = np.where(d["associative_mode_mg"] == "causally", 1.0, 0.50)
    d["route_base_weight"] = d["model_weight"] * d["association_weight"]
    d["mapping_conflict"] = False
    d["route_audit_status"] = "primary-paper experiment details pending"
    d["human_database_route_eligible"] = (
        (d["model_evidence_class"] == "human_only")
        & (d["target_models"] == "human")
        & d["transformation_pmids"].fillna("").ne("")
        & d["target_pmids"].fillna("").ne("")
        & (d["associative_mode_mg"] == "causally")
        & (~d["is_self_loop"])
    )
    d["route_exclusion_reason"] = np.select(
        [
            d["is_self_loop"],
            d["model_evidence_class"].eq("mixed_human_mouse"),
            d["model_evidence_class"].eq("mouse_only"),
            d["associative_mode_mg"].ne("causally"),
        ],
        [
            "substrate-product self-loop",
            "mixed human/mouse upstream evidence",
            "mouse-only upstream evidence",
            "target association deposited as correlational",
        ],
        default="eligible for audited human database subset",
    )
    keep = [
        "route_id", "route_type", "dedup_key", "parent_cid", "substrate_name",
        "product_cid", "product_name", "target_gene", "microbe_name",
        "upstream_models", "target_models", "transformation_pmids", "target_pmids",
        "all_source_ids", "source_databases", "raw_model_rows", "alteration",
        "associative_mode_mg", "model_evidence_class", "route_base_weight",
        "mapping_conflict", "is_self_loop", "route_audit_status",
        "human_database_route_eligible", "route_exclusion_reason",
    ]
    return d[keep].copy()


def build_host_routes() -> pd.DataFrame:
    host = pd.read_csv(HOST, low_memory=False)
    audit = pd.read_csv(EDGE_AUDIT, sep="\t", low_memory=False)
    transform_audit = audit[
        audit["relation_type"].eq("putative host enzymatic transformation (logical sub-edge)")
    ][["source_link_id", "verification_status", "audit_note"]].copy()
    host["source_link_id"] = "HOST-" + host["link_id"].astype(int).astype(str).str.zfill(3)
    d = host.merge(transform_audit, on="source_link_id", how="left", validate="one_to_one")
    d["route_id"] = d["source_link_id"]
    d["route_type"] = "host"
    d["parent_cid"] = cid_series(d["substrate_cid"])
    d["product_cid"] = cid_series(d["product_cid"])
    d["product_name"] = d["product"]
    d["mapping_conflict"] = d["verification_status"].eq("mapping conflict")
    d["is_self_loop"] = d["substrate_equals_metabolite"].fillna(False).astype(bool)
    d["route_base_weight"] = np.where(d["mapping_conflict"], 0.0, 0.20)
    d["model_evidence_class"] = np.where(
        d["mapping_conflict"], "excluded_mapping_conflict", "host_cascade_unverified"
    )
    d["human_database_route_eligible"] = False
    d["upstream_models"] = "host context not retained"
    d["target_models"] = "human target panel assumed"
    d["transformation_pmids"] = ""
    d["target_pmids"] = ""
    d["all_source_ids"] = ""
    d["source_databases"] = "BRENDA/HMDB/database bridge; terminal version incomplete"
    d["raw_model_rows"] = d["raw_record_count"]
    d["route_audit_status"] = d["verification_status"].fillna("not audited")
    d["route_exclusion_reason"] = np.where(
        d["mapping_conflict"],
        "audited EC/reaction mapping conflict; zero weight",
        "primary reaction and product-target provenance incomplete; sensitivity only",
    )
    d = d.rename(
        columns={
            "substrate_name": "substrate_name",
            "enzyme_name": "microbe_name",
            "action_type": "alteration",
        }
    )
    d["associative_mode_mg"] = "not retained"
    keep = [
        "route_id", "route_type", "dedup_key", "parent_cid", "substrate_name",
        "product_cid", "product_name", "target_gene", "microbe_name",
        "upstream_models", "target_models", "transformation_pmids", "target_pmids",
        "all_source_ids", "source_databases", "raw_model_rows", "alteration",
        "associative_mode_mg", "model_evidence_class", "route_base_weight",
        "mapping_conflict", "is_self_loop", "route_audit_status",
        "human_database_route_eligible", "route_exclusion_reason",
    ]
    return d[keep].copy()


def build_paths(occ: pd.DataFrame, routes: pd.DataFrame) -> pd.DataFrame:
    substrate_degree = occ.groupby("parent_cid")["herb_canonical"].nunique().to_dict()
    product_degree = routes.groupby("product_cid")["target_gene"].nunique().to_dict()
    paths = occ.merge(routes, on="parent_cid", how="inner", validate="many_to_many")
    paths["substrate_herb_degree"] = paths["parent_cid"].map(substrate_degree).astype(int)
    paths["product_target_degree"] = paths["product_cid"].map(product_degree).astype(int)
    paths["specificity_weight"] = 1.0 / np.sqrt(
        np.log1p(paths["substrate_herb_degree"].astype(float))
        * np.log1p(paths["product_target_degree"].astype(float))
    )
    paths["inclusive_path_weight"] = (
        paths["occurrence_weight"] * paths["route_base_weight"] * paths["specificity_weight"]
    )
    paths["primary_path_eligible"] = (
        paths["human_database_route_eligible"]
        & paths["lotus_verified_any"]
        & (~paths["mapping_conflict"])
    )
    paths["primary_path_weight"] = np.where(
        paths["primary_path_eligible"], paths["inclusive_path_weight"], 0.0
    )
    paths["path_exclusion_reason"] = np.where(
        ~paths["lotus_verified_any"],
        "herb-compound occurrence lacks LOTUS cross-verification",
        paths["route_exclusion_reason"],
    )
    return paths.sort_values(["herb_canonical", "route_id"]).reset_index(drop=True)


def aggregate_score(
    paths: pd.DataFrame,
    herbs: pd.DataFrame,
    weight_col: str,
    eligible_mask: pd.Series | None = None,
) -> tuple[pd.DataFrame, dict[str, float]]:
    p = paths if eligible_mask is None else paths[eligible_mask].copy()
    agg = (
        p.groupby("herb_canonical")
        .agg(
            transform_score_raw=(weight_col, "sum"),
            n_routes=("route_id", "nunique"),
            n_route_substrates=("parent_cid", "nunique"),
            n_products=("product_cid", "nunique"),
            n_transform_targets=("target_gene", "nunique"),
        )
        .reset_index()
    )
    out = herbs.merge(agg, on="herb_canonical", how="left")
    for c in ["transform_score_raw", "n_routes", "n_route_substrates", "n_products", "n_transform_targets"]:
        out[c] = out[c].fillna(0)
    out["transform_score_log"] = np.log1p(out["transform_score_raw"].astype(float))
    resid, adj_z, reg = degree_adjust(
        out["transform_score_log"], out["n_compounds"], out["n_route_substrates"]
    )
    out["transform_degree_residual"] = resid
    out["transform_score_adjusted_z"] = adj_z
    out["transform_score_adjusted_01"] = minmax(adj_z)
    out["transform_rank"] = rank_desc(adj_z)
    return out, reg


def count_sources(series: pd.Series) -> int:
    return len(split_tokens(series))


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    """Render a compact Markdown table without optional third-party packages."""
    cols = [str(c) for c in df.columns]
    lines = ["|" + "|".join(cols) + "|", "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in df.itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, float):
                text = f"{value:.6g}"
            else:
                text = str(value)
            values.append(text.replace("|", "\\|"))
        lines.append("|" + "|".join(values) + "|")
    return "\n".join(lines)


def make_herb_features(
    baseline: pd.DataFrame,
    herb_degree: pd.DataFrame,
    paths: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    base = baseline.merge(herb_degree, on="herb_canonical", how="left", validate="one_to_one")
    base["frozen_rank"] = base["kg_rank"].astype(int)
    base["frozen_direct_score"] = base["composite_score"].astype(float)
    if base["n_compounds"].isna().any():
        missing = base.loc[base["n_compounds"].isna(), "herb_canonical"].tolist()
        raise ValueError(f"Missing compound degree for herbs: {missing}")

    direct_resid, direct_z, direct_reg = degree_adjust(
        base["composite_score"], base["n_compounds"], base["n_t2d_targets"],
        degree_2_label="n_t2d_targets",
    )
    base["direct_degree_residual"] = direct_resid
    base["direct_adjusted_z"] = direct_z
    base["direct_adjusted_rank"] = rank_desc(direct_z)

    primary, primary_reg = aggregate_score(
        paths, base[["herb_canonical", "n_compounds"]], "primary_path_weight",
        paths["primary_path_eligible"],
    )
    inclusive, inclusive_reg = aggregate_score(
        paths, base[["herb_canonical", "n_compounds"]], "inclusive_path_weight",
        ~paths["mapping_conflict"],
    )
    primary = primary.rename(columns={c: f"primary_{c}" for c in primary.columns if c not in {"herb_canonical", "n_compounds"}})
    inclusive = inclusive.rename(columns={c: f"inclusive_{c}" for c in inclusive.columns if c not in {"herb_canonical", "n_compounds"}})
    out = base.merge(primary.drop(columns="n_compounds"), on="herb_canonical", how="left")
    out = out.merge(inclusive.drop(columns="n_compounds"), on="herb_canonical", how="left")

    descriptive = []
    for herb in out["herb_canonical"]:
        g = paths[paths["herb_canonical"] == herb]
        descriptive.append(
            {
                "herb_canonical": herb,
                "n_all_paths": len(g),
                "n_unique_routes_all": g["route_id"].nunique(),
                "n_human_only_microbial_routes": g.loc[g["model_evidence_class"].eq("human_only"), "route_id"].nunique(),
                "n_mixed_microbial_routes": g.loc[g["model_evidence_class"].eq("mixed_human_mouse"), "route_id"].nunique(),
                "n_mouse_only_microbial_routes": g.loc[g["model_evidence_class"].eq("mouse_only"), "route_id"].nunique(),
                "n_host_unverified_routes": g.loc[g["model_evidence_class"].eq("host_cascade_unverified"), "route_id"].nunique(),
                "n_mapping_conflict_routes": g.loc[g["mapping_conflict"], "route_id"].nunique(),
                "n_source_identifiers": count_sources(g["all_source_ids"]),
            }
        )
    out = out.merge(pd.DataFrame(descriptive), on="herb_canonical", how="left")

    out["evidence_tier"] = np.select(
        [
            out["primary_n_routes"] > 0,
            out["n_mixed_microbial_routes"] > 0,
            (out["n_mouse_only_microbial_routes"] > 0) | (out["n_host_unverified_routes"] > 0),
        ],
        [1, 2, 3],
        default=4,
    ).astype(int)
    out["evidence_tier_label"] = out["evidence_tier"].map(
        {
            1: "audited human database subset available",
            2: "mixed human/mouse microbial evidence only",
            3: "mouse-only or provenance-incomplete host evidence only",
            4: "no non-conflicting CID-connected route",
        }
    )
    out = out.sort_values(
        ["evidence_tier", "primary_transform_score_adjusted_z", "inclusive_transform_score_adjusted_z", "direct_adjusted_z"],
        ascending=[True, False, False, False],
    ).reset_index(drop=True)
    out["tiered_research_priority_rank"] = np.arange(1, len(out) + 1)
    return out, {"direct": direct_reg, "primary_transform": primary_reg, "inclusive_transform": inclusive_reg}


def alpha_sensitivity(features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    metrics = []
    direct = features["direct_adjusted_z"].to_numpy(float)
    transform = features["primary_transform_score_adjusted_z"].to_numpy(float)
    direct_rank = rank_desc(direct)
    direct_top = set(features.loc[direct_rank <= TOP_N, "herb_canonical"])
    for alpha in ALPHAS:
        score = (1.0 - alpha) * direct + alpha * transform
        rank = rank_desc(score)
        top = set(features.loc[rank <= TOP_N, "herb_canonical"])
        rho = spearmanr(direct_rank, rank)
        for herb, value, rk in zip(features["herb_canonical"], score, rank):
            rows.append(
                {
                    "alpha": alpha,
                    "herb_canonical": herb,
                    "exploratory_composite_score": value,
                    "exploratory_rank": int(rk),
                }
            )
        metrics.append(
            {
                "alpha": alpha,
                "spearman_vs_direct_adjusted": float(rho.statistic),
                "top20_overlap_with_direct_adjusted": len(top & direct_top),
                "top20_jaccard_with_direct_adjusted": len(top & direct_top) / len(top | direct_top),
                "alpha_selected_using_clinical_data": False,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(metrics)


def scenario_sensitivity(paths: pd.DataFrame, base: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scenarios = {
        "primary_human_lotus_causal": paths["primary_path_eligible"],
        "human_any_occurrence_causal": (
            paths["model_evidence_class"].eq("human_only")
            & paths["associative_mode_mg"].eq("causally")
            & (~paths["mapping_conflict"])
        ),
        "human_plus_mixed_lotus_causal": (
            paths["model_evidence_class"].isin(["human_only", "mixed_human_mouse"])
            & paths["associative_mode_mg"].eq("causally")
            & paths["lotus_verified_any"]
            & (~paths["mapping_conflict"])
        ),
        "all_microbial_nonconflict": paths["route_type"].eq("microbial") & (~paths["mapping_conflict"]),
        "inclusive_nonconflict": ~paths["mapping_conflict"],
    }
    long_rows = []
    primary_rank: np.ndarray | None = None
    primary_top: set[str] | None = None
    metrics = []
    for name, mask in scenarios.items():
        ranked, _ = aggregate_score(
            paths,
            base[["herb_canonical", "n_compounds"]],
            "inclusive_path_weight",
            mask,
        )
        rank = ranked["transform_rank"].to_numpy(int)
        top = set(ranked.loc[rank <= TOP_N, "herb_canonical"])
        if primary_rank is None:
            primary_rank = rank.copy()
            primary_top = set(top)
        rho = spearmanr(primary_rank, rank)
        metrics.append(
            {
                "scenario": name,
                "n_eligible_paths": int(mask.sum()),
                "n_herbs_with_routes": int(ranked["n_routes"].gt(0).sum()),
                "spearman_vs_primary": float(rho.statistic),
                "top20_overlap_vs_primary": len(top & primary_top),
            }
        )
        for _, row in ranked.iterrows():
            long_rows.append(
                {
                    "scenario": name,
                    "herb_canonical": row["herb_canonical"],
                    "transform_score_raw": row["transform_score_raw"],
                    "transform_score_adjusted_z": row["transform_score_adjusted_z"],
                    "transform_rank": int(row["transform_rank"]),
                    "n_routes": int(row["n_routes"]),
                }
            )
    return pd.DataFrame(long_rows), pd.DataFrame(metrics)


def degrees(left: np.ndarray, right: np.ndarray) -> tuple[dict[str, int], dict[str, int]]:
    dl: defaultdict[str, int] = defaultdict(int)
    dr: defaultdict[str, int] = defaultdict(int)
    for a, b in zip(left.tolist(), right.tolist()):
        dl[str(a)] += 1
        dr[str(b)] += 1
    return dict(dl), dict(dr)


def swap_bipartite(
    left: np.ndarray, right: np.ndarray, n_swaps: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, int]:
    left = left.copy()
    right = right.copy()
    edges = set(zip(left.tolist(), right.tolist()))
    accepted = 0
    for _ in range(n_swaps):
        i, j = rng.integers(0, len(left), size=2)
        if i == j:
            continue
        a, x = left[i], right[i]
        b, y = left[j], right[j]
        if a == b or x == y or (a, y) in edges or (b, x) in edges:
            continue
        edges.remove((a, x))
        edges.remove((b, y))
        edges.add((a, y))
        edges.add((b, x))
        right[i], right[j] = y, x
        accepted += 1
    return left, right, accepted


def topology_null(
    paths: pd.DataFrame,
    occurrences: pd.DataFrame,
    features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible_routes = paths.loc[paths["human_database_route_eligible"], [
        "route_id", "parent_cid", "product_cid", "target_gene", "route_base_weight",
    ]].drop_duplicates()
    occ = occurrences[occurrences["lotus_verified_any"]][["herb_canonical", "parent_cid"]].drop_duplicates()
    route_cids = set(eligible_routes["parent_cid"])
    occ = occ[occ["parent_cid"].isin(route_cids)].copy()
    left0 = occ["herb_canonical"].to_numpy(object)
    right0 = occ["parent_cid"].to_numpy(object)
    dl0, dr0 = degrees(left0, right0)
    route_by_cid = {cid: g.copy() for cid, g in eligible_routes.groupby("parent_cid")}
    product_degree = eligible_routes.groupby("product_cid")["target_gene"].nunique().to_dict()
    direct_map = features.set_index("herb_canonical")["direct_adjusted_z"].to_dict()
    n_comp_map = features.set_index("herb_canonical")["n_compounds"].to_dict()
    herbs = features["herb_canonical"].tolist()

    def score_edges(left: np.ndarray, right: np.ndarray) -> pd.DataFrame:
        rows = []
        substrate_deg = defaultdict(int)
        for h, c in zip(left.tolist(), right.tolist()):
            substrate_deg[c] += 1
        for herb, cid in zip(left.tolist(), right.tolist()):
            rg = route_by_cid.get(cid)
            if rg is None:
                continue
            for _, rr in rg.iterrows():
                specificity = 1.0 / math.sqrt(
                    math.log1p(substrate_deg[cid]) * math.log1p(product_degree[rr["product_cid"]])
                )
                rows.append((herb, cid, rr["route_id"], float(rr["route_base_weight"]) * specificity))
        q = pd.DataFrame(rows, columns=["herb_canonical", "parent_cid", "route_id", "weight"])
        if q.empty:
            agg = pd.DataFrame({"herb_canonical": herbs, "raw": 0.0, "opp": 0})
        else:
            agg = q.groupby("herb_canonical").agg(raw=("weight", "sum"), opp=("parent_cid", "nunique")).reset_index()
            agg = pd.DataFrame({"herb_canonical": herbs}).merge(agg, on="herb_canonical", how="left").fillna(0)
        y = np.log1p(agg["raw"].to_numpy(float))
        _, adj, _ = degree_adjust(
            y,
            np.array([n_comp_map[h] for h in herbs], dtype=float),
            agg["opp"].to_numpy(float),
            degree_2_label="n_route_substrates",
        )
        agg["adjusted"] = adj
        agg["direct"] = [direct_map[h] for h in herbs]
        return agg

    observed = score_edges(left0, right0)
    obs_rho = float(spearmanr(observed["direct"], observed["adjusted"]).statistic)
    obs_top = set(observed.nlargest(TOP_N, "adjusted")["herb_canonical"])
    direct_top = set(observed.nlargest(TOP_N, "direct")["herb_canonical"])
    obs_overlap = len(obs_top & direct_top)

    rng = np.random.default_rng(SEED)
    rows = []
    for i in range(N_NULL):
        l, r, accepted = swap_bipartite(left0, right0, 10 * len(left0), rng)
        dl1, dr1 = degrees(l, r)
        if dl1 != dl0 or dr1 != dr0:
            raise RuntimeError(f"degree preservation failed at null iteration {i}")
        q = score_edges(l, r)
        rho = float(spearmanr(q["direct"], q["adjusted"]).statistic)
        top = set(q.nlargest(TOP_N, "adjusted")["herb_canonical"])
        rows.append(
            {
                "iteration": i,
                "spearman_direct_vs_transform": rho,
                "top20_overlap_direct_vs_transform": len(top & direct_top),
                "accepted_swaps": accepted,
                "left_degree_preserved": True,
                "right_degree_preserved": True,
            }
        )
    null = pd.DataFrame(rows)
    p_rho = (np.sum(np.abs(null["spearman_direct_vs_transform"]) >= abs(obs_rho)) + 1) / (N_NULL + 1)
    p_overlap = (np.sum(null["top20_overlap_direct_vs_transform"] >= obs_overlap) + 1) / (N_NULL + 1)
    summary = pd.DataFrame(
        [
            {
                "metric": "spearman_direct_vs_transform",
                "observed": obs_rho,
                "null_mean": null["spearman_direct_vs_transform"].mean(),
                "null_sd": null["spearman_direct_vs_transform"].std(ddof=1),
                "empirical_p": p_rho,
                "tail": "two-sided absolute",
            },
            {
                "metric": "top20_overlap_direct_vs_transform",
                "observed": obs_overlap,
                "null_mean": null["top20_overlap_direct_vs_transform"].mean(),
                "null_sd": null["top20_overlap_direct_vs_transform"].std(ddof=1),
                "empirical_p": p_overlap,
                "tail": "greater or equal",
            },
        ]
    )
    return null, summary


def correlations(features: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("frozen_rank", "n_compounds"),
        ("frozen_rank", "n_t2d_targets"),
        ("direct_adjusted_rank", "n_compounds"),
        ("direct_adjusted_rank", "n_t2d_targets"),
        ("primary_transform_rank", "n_compounds"),
        ("primary_transform_rank", "primary_n_route_substrates"),
        ("tiered_research_priority_rank", "n_compounds"),
    ]
    rows = []
    for a, b in pairs:
        stat = spearmanr(features[a], features[b])
        rows.append({"variable_1": a, "variable_2": b, "spearman_rho": stat.statistic, "p_value_descriptive": stat.pvalue})
    return pd.DataFrame(rows)


def write_report(
    out: Path,
    baseline: pd.DataFrame,
    routes: pd.DataFrame,
    occurrences: pd.DataFrame,
    paths: pd.DataFrame,
    features: pd.DataFrame,
    alpha_metrics: pd.DataFrame,
    scenario_metrics: pd.DataFrame,
    null_summary: pd.DataFrame,
    regressions: dict[str, dict[str, float]],
    checks: dict[str, object],
) -> None:
    tier_counts = features["evidence_tier_label"].value_counts().to_dict()
    top = features.nsmallest(20, "tiered_research_priority_rank")
    lines = [
        "# Metabolism-aware v2 execution report",
        "",
        "## Outcome",
        "",
        "The v2 pipeline was implemented and executed without changing the frozen v1 ranking. "
        "It now reads CID-connected microbial and host transformation routes, separates evidence classes, "
        "downweights broad topology and emits two-axis and sensitivity rankings.",
        "",
        "This is an evidence-bounded hypothesis-prioritisation framework. It is not a validated predictor "
        "of biochemical conversion, clinical efficacy or patient-specific response.",
        "",
        "## Input and route accounting",
        "",
        f"- Frozen FMH universe: **{len(baseline)}** materials.",
        f"- Audited terminal routes read: **{len(routes)}** ({(routes.route_type == 'microbial').sum()} microbial; {(routes.route_type == 'host').sum()} host).",
        f"- Official herb-compound CID occurrences: **{len(occurrences)}** unique herb-CID pairs.",
        f"- Constructed herb-route paths: **{len(paths)}**.",
        f"- Audited human-database-subset eligible paths: **{int(paths.primary_path_eligible.sum())}**.",
        f"- Mapping-conflict paths retained with zero weight: **{int(paths.mapping_conflict.sum())}**.",
        "",
        "## Evidence-tier distribution",
        "",
    ]
    for label, n in tier_counts.items():
        lines.append(f"- {label}: **{n}** materials.")
    lines += [
        "",
        "The first tier does not mean that the complete route has been adjudicated as a human mechanism. "
        "It means that the deposited microbial record is human-only, source identifiers are retained, the "
        "target-side record is human and causal as deposited, and the herb-compound occurrence is LOTUS-cross-verified.",
        "",
        "## Evidence-tier-first top 20 research priorities",
        "",
        "|rank|material|tier|primary routes|primary products|primary targets|direct adjusted rank|",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in top.iterrows():
        lines.append(
            f"|{int(r.tiered_research_priority_rank)}|{r.herb_canonical}|{int(r.evidence_tier)}|"
            f"{int(r.primary_n_routes)}|{int(r.primary_n_products)}|{int(r.primary_n_transform_targets)}|"
            f"{int(r.direct_adjusted_rank)}|"
        )
    lines += [
        "",
        "## Alpha sensitivity",
        "",
        dataframe_to_markdown(alpha_metrics),
        "",
        "No alpha was selected on the basis of RCT results.",
        "",
        "## Evidence-scenario sensitivity",
        "",
        dataframe_to_markdown(scenario_metrics),
        "",
        "## Degree-preserving topology null",
        "",
        dataframe_to_markdown(null_summary),
        "",
        "The topology null tests overlap between the new transformation evidence axis and the direct-edge "
        "axis. It is not a test of treatment efficacy or biochemical prediction accuracy.",
        "",
        "## Regression diagnostics",
        "",
        "```json",
        json.dumps(regressions, ensure_ascii=False, indent=2),
        "```",
        "",
        "## Validation checks",
        "",
        "```json",
        json.dumps(checks, ensure_ascii=False, indent=2),
        "```",
        "",
        "## Defensible conclusion",
        "",
        "The implementation gap has been repaired at the computational level: transformation routes now "
        "enter a separately reported evidence axis and prespecified sensitivity rerankings. The analysis "
        "does not establish predictive superiority because no independent transformation ground truth was "
        "available and critical route provenance remains incomplete.",
        "",
        "## Claims that remain prohibited",
        "",
        "- The v2 rank predicts clinical efficacy.",
        "- The ranked materials increase Rb1-to-compound-K conversion in T2D.",
        "- Docking, homology or database routes demonstrate enzyme activity.",
        "- Lack of an eligible route is evidence of biological absence.",
        "- The ranking is a patient-specific recommendation.",
        "",
        "## Remaining evidence requirement",
        "",
        "A fully verified primary route score still requires source-level adjudication of herb occurrence, "
        "transformation assays and product-target evidence. Independent donor-level Rb1/Rd/F2/compound-K "
        "measurements or enzyme assays are required for biochemical validation.",
    ]
    (out / "04_reports/METABOLISM_AWARE_V2_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    global BASELINE, HERB_COMPOUND, MICRO, MICRO_RAW, HOST, EDGE_AUDIT, ENTITY_AUDIT

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing the seven required metabolism-aware v2 input files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for regenerated metabolism-aware v2 results.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()

    BASELINE = input_dir / "KG_all_FMH_ranked.csv"
    HERB_COMPOUND = input_dir / "2_3_herb_compound_clean.csv"
    MICRO = input_dir / "Supplementary_Table_Xb_gut_microbial_links.csv"
    MICRO_RAW = input_dir / "Supplementary_Table_Xb_raw_microbial_329.csv"
    HOST = input_dir / "Supplementary_Table_Xc_host_enzymatic_links.csv"
    EDGE_AUDIT = input_dir / "KG_edge_evidence.tsv"
    ENTITY_AUDIT = input_dir / "KG_entity_and_mapping_audit.tsv"

    required_inputs = [
        BASELINE,
        HERB_COMPOUND,
        MICRO,
        MICRO_RAW,
        HOST,
        EDGE_AUDIT,
        ENTITY_AUDIT,
    ]

    missing = [str(x) for x in required_inputs if not x.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing required input files:\\n  " + "\\n  ".join(missing)
        )
    out = args.output.resolve()
    for sub in ["01_routes", "02_rankings", "03_validation", "04_reports", "logs"]:
        (out / sub).mkdir(parents=True, exist_ok=True)

    baseline = pd.read_csv(BASELINE, low_memory=False)
    baseline["herb_canonical"] = baseline["herb_canonical"].astype(str).str.strip()
    if len(baseline) != 104 or baseline["herb_canonical"].nunique() != 104:
        raise ValueError("Frozen baseline is not the expected 104-unique-material universe")
    expected = baseline.sort_values(["composite_score", "herb_canonical"], ascending=[False, True])["herb_canonical"].tolist()
    reported = baseline.sort_values("kg_rank")["herb_canonical"].tolist()
    frozen_order_matches_score = expected == reported

    occurrences, herb_degree = build_occurrences()
    microbial = build_microbial_routes()
    host = build_host_routes()
    routes = pd.concat([microbial, host], ignore_index=True)
    if routes["route_id"].duplicated().any():
        raise ValueError("Duplicate route IDs")
    paths = build_paths(occurrences, routes)

    features, regressions = make_herb_features(baseline, herb_degree, paths)
    alpha_long, alpha_metrics = alpha_sensitivity(features)
    scenario_long, scenario_metrics = scenario_sensitivity(paths, features)
    null, null_summary = topology_null(paths, occurrences, features)
    corr = correlations(features)

    alpha0 = alpha_long[alpha_long["alpha"] == 0].sort_values("herb_canonical")
    direct = features.sort_values("herb_canonical")
    alpha0_matches = np.array_equal(
        alpha0["exploratory_rank"].to_numpy(int), direct["direct_adjusted_rank"].to_numpy(int)
    )
    conflict_zero = bool(paths.loc[paths["mapping_conflict"], "inclusive_path_weight"].eq(0).all())
    checks = {
        "frozen_baseline_rows": len(baseline),
        "frozen_baseline_unique_materials": baseline["herb_canonical"].nunique(),
        "frozen_reported_order_matches_frozen_score": frozen_order_matches_score,
        "route_ids_unique": not routes["route_id"].duplicated().any(),
        "microbial_terminal_routes": int((routes["route_type"] == "microbial").sum()),
        "host_terminal_routes": int((routes["route_type"] == "host").sum()),
        "mapping_conflict_routes": int(routes["mapping_conflict"].sum()),
        "mapping_conflict_paths_have_zero_weight": conflict_zero,
        "alpha_zero_matches_direct_adjusted_comparator": alpha0_matches,
        "all_104_materials_retained": len(features) == 104 and features["herb_canonical"].nunique() == 104,
        "clinical_data_used_for_tuning": False,
    }
    if not all([
        frozen_order_matches_score,
        conflict_zero,
        alpha0_matches,
        checks["all_104_materials_retained"],
    ]):
        raise RuntimeError(f"Critical validation check failed: {checks}")

    routes.to_csv(out / "01_routes/transformation_routes_audited_v2.tsv", sep="\t", index=False)
    occurrences.to_csv(out / "01_routes/herb_compound_occurrences_v2.tsv", sep="\t", index=False)
    paths.to_csv(out / "01_routes/herb_transformation_paths_v2.tsv", sep="\t", index=False)
    paths[~paths["primary_path_eligible"]].to_csv(
        out / "01_routes/KG_v2_excluded_or_sensitivity_paths.tsv", sep="\t", index=False
    )

    features.to_csv(out / "02_rankings/KG_metabolism_aware_v2_two_axis.tsv", sep="\t", index=False)
    features.sort_values("tiered_research_priority_rank").to_csv(
        out / "02_rankings/KG_metabolism_aware_v2_tiered_ranking.tsv", sep="\t", index=False
    )
    alpha_long.to_csv(out / "02_rankings/KG_metabolism_aware_v2_alpha_sensitivity.tsv", sep="\t", index=False)
    scenario_long.to_csv(out / "02_rankings/KG_metabolism_aware_v2_evidence_scenarios.tsv", sep="\t", index=False)

    alpha_metrics.to_csv(out / "03_validation/alpha_sensitivity_summary.tsv", sep="\t", index=False)
    scenario_metrics.to_csv(out / "03_validation/evidence_scenario_summary.tsv", sep="\t", index=False)
    null.to_csv(out / "03_validation/topology_null_1000.tsv", sep="\t", index=False)
    null_summary.to_csv(out / "03_validation/topology_null_summary.tsv", sep="\t", index=False)
    corr.to_csv(out / "03_validation/degree_and_rank_correlations.tsv", sep="\t", index=False)
    (out / "03_validation/regression_diagnostics.json").write_text(
        json.dumps(regressions, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (out / "03_validation/validation_checks.json").write_text(
        json.dumps(checks, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    write_report(
        out, baseline, routes, occurrences, paths, features, alpha_metrics,
        scenario_metrics, null_summary, regressions, checks,
    )

    inputs = [BASELINE, HERB_COMPOUND, MICRO, MICRO_RAW, HOST, EDGE_AUDIT, ENTITY_AUDIT]
    manifest = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": __import__("scipy").__version__,
        "seed": SEED,
        "null_iterations": N_NULL,
        "inputs": [{"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} for p in inputs],
        "checks": checks,
    }
    (out / "logs/run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    output_rows = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "OUTPUT_SHA256.tsv":
            output_rows.append({"relative_path": str(p.relative_to(out)), "bytes": p.stat().st_size, "sha256": sha256(p)})
    pd.DataFrame(output_rows).to_csv(out / "04_reports/OUTPUT_SHA256.tsv", sep="\t", index=False)
    print(json.dumps({"output": str(out), "checks": checks}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
