#!/usr/bin/env python3
"""Build the gutMGene–AGORA2 reaction and flux-feasibility matrix.

The analysis is evidence-bounded.  It does not add missing reactions and it does
not treat host target genes as microbial model genes.  Flux tests use an
artificial substrate-only-carbon medium and therefore test model feasibility,
not in-vivo activity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import sys
import time
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import cobra
import networkx as nx
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import read_sbml_model


CHAIN = None
T2D_LINKS = None
PUBCHEM = None
SBML_DIR = None
BRIDGE = None
DEFAULT_OUT = Path("outputs")
TOL = 1e-8


# NCBI Taxonomy-verified legacy names used by AGORA2 v2.01, plus two explicit
# source spelling corrections.  Broad family/phylum/no-rank records are not
# forced onto strain models.
TAXON_ALIASES: dict[str, dict[str, str]] = {
    "Acidipropionibacterium acidipropionici": {
        "prefix": "Propionibacterium_acidipropionici_",
        "basis": "NCBI TaxID 1748 synonym",
    },
    "Agathobacter rectalis": {
        "prefix": "Eubacterium_rectale_",
        "basis": "NCBI TaxID 39491 synonym",
    },
    "Anaerobutyricum hallii": {
        "prefix": "Eubacterium_hallii_",
        "basis": "NCBI TaxID 39488 synonym",
    },
    "Enterocloster asparagiformis": {
        "prefix": "Clostridium_asparagiforme_",
        "basis": "NCBI TaxID 333367 synonym",
    },
    "Enterocloster bolteae": {
        "prefix": "Clostridium_bolteae_",
        "basis": "NCBI TaxID 208479 synonym",
    },
    "Enterocloster citroniae": {
        "prefix": "Clostridium_citroniae_",
        "basis": "NCBI TaxID 358743 synonym",
    },
    "Enterocloster clostridioformis": {
        "prefix": "Clostridium_clostridioforme_",
        "basis": "NCBI TaxID 1531 synonym",
    },
    "Holdemanella biformis": {
        "prefix": "Eubacterium_biforme_",
        "basis": "legacy AGORA2 name for NCBI TaxID 1735",
    },
    "Hungatella hathewayi": {
        "prefix": "Clostridium_hathewayi_",
        "basis": "NCBI TaxID 154046 synonym",
    },
    "Lacticaseibacillus paracasei": {
        "prefix": "Lactobacillus_paracasei_",
        "basis": "NCBI TaxID 1597 synonym",
    },
    "Lactiplantibacillus plantarum": {
        "prefix": "Lactobacillus_plantarum_",
        "basis": "NCBI TaxID 1590 synonym",
    },
    "Ligilactobacillus murinus": {
        "prefix": "Lactobacillus_murinus_",
        "basis": "NCBI TaxID 1622 synonym",
    },
    "Ligilactobacillus salivarius": {
        "prefix": "Lactobacillus_salivarius_",
        "basis": "NCBI TaxID 1624 synonym",
    },
    "Limosilactobacillus mucosae": {
        "prefix": "Lactobacillus_mucosae_",
        "basis": "NCBI TaxID 97478 synonym",
    },
    "[Clostridium] scindens": {
        "prefix": "Clostridium_scindens_",
        "basis": "NCBI TaxID 29347 synonym",
    },
    "[Ruminococcus] torques": {
        "prefix": "Ruminococcus_torques_",
        "basis": "NCBI TaxID 33039 synonym",
    },
    "Lactobacilli reuteri": {
        "prefix": "Lactobacillus_reuteri_",
        "basis": "source spelling correction",
    },
    "[Clostridium] nexile 1": {
        "prefix": "Clostridium_nexile_",
        "basis": "source bracket/suffix normalization",
    },
    "Bififidobacterium thetaiotaomicron": {
        "prefix": "Bacteroides_thetaiotaomicron_",
        "basis": "source spelling/genus correction; strain DSM 2079",
    },
}

CURRENCY_BASE_IDS = {
    "h",
    "h2o",
    "atp",
    "adp",
    "amp",
    "pi",
    "ppi",
    "nad",
    "nadh",
    "nadp",
    "nadph",
    "coa",
    "accoa",
    "fad",
    "fadh2",
    "q8",
    "q8h2",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def norm_name(value: object) -> str:
    text = clean_text(value).lower()
    text = text.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    return re.sub(r"[^a-z0-9]+", "", text)


def as_tokens(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if clean_text(item)]
    text = clean_text(value)
    if not text or text.lower() in {"nan", "none"}:
        return []
    return [text]


def parse_formula(formula: object) -> dict[str, int]:
    text = clean_text(formula)
    counts: dict[str, int] = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", text):
        counts[element] = counts.get(element, 0) + int(count or 1)
    return counts


def base_metabolite_id(met_id: str) -> str:
    return re.sub(r"\[[^\]]+\]$", "", met_id)


def build_taxon_mapping(chain: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    bridge = pd.read_csv(BRIDGE, encoding="utf-8-sig", low_memory=False)
    taxa = (
        chain[["microbe_name", "microbe_ncbi_id", "microbe_rank"]]
        .drop_duplicates("microbe_name")
        .sort_values("microbe_name")
    )
    bridge_by_name = bridge.set_index("gutmgene_name", drop=False)
    all_files = sorted(path.name for path in SBML_DIR.glob("*.xml"))
    mapping_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []

    for row in taxa.itertuples(index=False):
        name = str(row.microbe_name)
        model_files: list[str] = []
        match_level = "none"
        match_basis = "no reliable AGORA2 strain mapping"
        analysis_eligible = False
        if name in bridge_by_name.index:
            hit = bridge_by_name.loc[name]
            if isinstance(hit, pd.DataFrame):
                hit = hit.iloc[0]
            if clean_text(hit.get("agora2_files", "")):
                model_files = [
                    part.strip()
                    for part in str(hit["agora2_files"]).split(";")
                    if part.strip()
                ]
                match_level = str(hit["match_level"])
                match_basis = "original AGORA2 filename bridge"
                analysis_eligible = match_level in {"species", "genus"}
        if not model_files and name in TAXON_ALIASES:
            alias = TAXON_ALIASES[name]
            model_files = [f for f in all_files if f.startswith(alias["prefix"])]
            if model_files:
                match_level = "species_alias"
                match_basis = alias["basis"]
                analysis_eligible = True

        existing = [f for f in model_files if (SBML_DIR / f).is_file()]
        missing = sorted(set(model_files) - set(existing))
        mapping_rows.append(
            {
                "microbe_name": name,
                "microbe_ncbi_id": row.microbe_ncbi_id,
                "microbe_rank": row.microbe_rank,
                "model_match_level": match_level,
                "model_match_basis": match_basis,
                "model_count": len(existing),
                "analysis_eligible": analysis_eligible and bool(existing),
                "missing_model_file_count": len(missing),
                "model_files_pipe": "|".join(existing),
            }
        )
        if existing:
            for filename in existing:
                edge_rows.append(
                    {
                        "microbe_name": name,
                        "microbe_ncbi_id": row.microbe_ncbi_id,
                        "microbe_rank": row.microbe_rank,
                        "model_match_level": match_level,
                        "model_match_basis": match_basis,
                        "analysis_eligible": analysis_eligible,
                        "agora2_model_file": filename,
                        "agora2_model_path": str(SBML_DIR / filename),
                    }
                )
        else:
            edge_rows.append(
                {
                    "microbe_name": name,
                    "microbe_ncbi_id": row.microbe_ncbi_id,
                    "microbe_rank": row.microbe_rank,
                    "model_match_level": match_level,
                    "model_match_basis": match_basis,
                    "analysis_eligible": False,
                    "agora2_model_file": "",
                    "agora2_model_path": "",
                }
            )

    summary = pd.DataFrame(mapping_rows)
    edges = pd.DataFrame(edge_rows)
    if len(summary) != 147 or summary["microbe_name"].nunique() != 147:
        raise ValueError("Taxon mapping must contain exactly 147 gutMGene classifications")
    return summary, edges


def build_routes(chain: pd.DataFrame) -> pd.DataFrame:
    chain = chain.copy()
    chain["substrate_pubchem_cid"] = pd.to_numeric(
        chain["substrate_pubchem_cid"], errors="coerce"
    ).astype("Int64")
    chain["metabolite_pubchem_cid"] = pd.to_numeric(
        chain["metabolite_pubchem_cid"], errors="coerce"
    ).astype("Int64")
    keys = [
        "microbe_name",
        "substrate_pubchem_cid",
        "substrate_name",
        "metabolite_pubchem_cid",
        "metabolite_name",
    ]

    def join_unique(values: pd.Series) -> str:
        return "|".join(sorted({clean_text(v) for v in values if clean_text(v)}))

    routes = (
        chain.groupby(keys, dropna=False, sort=True)
        .agg(
            path_ids=("path_id", join_unique),
            target_genes=("target_gene", join_unique),
            transformation_pmids=("pmid", join_unique),
            target_pmids=("pmid_mg", join_unique),
            current_fmh_herbs=("current_fmh_herbs", join_unique),
            raw_chain_rows=("path_id", "size"),
            any_both_edges_causal=("both_edges_causal", "max"),
            any_host_species_concordant=("host_species_concordant", "max"),
        )
        .reset_index()
    )
    routes.insert(0, "taxon_route_id", [f"TR{i:04d}" for i in range(1, len(routes) + 1)])
    t2d = pd.read_csv(T2D_LINKS, low_memory=False)
    t2d_keys = {
        (int(r.substrate_pubchem_cid), int(r.metabolite_pubchem_cid_mg), str(r.target_gene))
        for r in t2d.itertuples(index=False)
        if pd.notna(r.substrate_pubchem_cid) and pd.notna(r.metabolite_pubchem_cid_mg)
    }
    routes["t2d_core_target_genes"] = routes.apply(
        lambda row: "|".join(
            sorted(
                gene
                for gene in clean_text(row["target_genes"]).split("|")
                if gene
                and pd.notna(row["substrate_pubchem_cid"])
                and pd.notna(row["metabolite_pubchem_cid"])
                and (
                    int(row["substrate_pubchem_cid"]),
                    int(row["metabolite_pubchem_cid"]),
                    gene,
                )
                in t2d_keys
            )
        ),
        axis=1,
    )
    routes["in_t2d_core_75"] = routes["t2d_core_target_genes"].str.len().gt(0)
    return routes


def mapping_candidates(model: cobra.Model, cid: object, name: object, inchi: object) -> list[dict[str, Any]]:
    if pd.isna(cid):
        return []
    target_cid = str(int(cid))
    target_name = norm_name(name)
    target_inchi = clean_text(inchi)
    candidates: dict[str, dict[str, Any]] = {}
    for met in model.metabolites:
        annotations = met.annotation or {}
        model_cids = as_tokens(annotations.get("pubchem.compound"))
        model_inchi = clean_text(annotations.get("inchi"))
        model_name = norm_name(met.name)
        methods: list[str] = []
        if target_cid in model_cids:
            methods.append("exact_pubchem_cid")
        if target_inchi and model_inchi and target_inchi == model_inchi:
            methods.append("exact_inchi")
        if target_name and model_name and target_name == model_name:
            methods.append("exact_normalized_name")
        if not methods:
            continue
        priority = min(
            {"exact_pubchem_cid": 0, "exact_inchi": 1, "exact_normalized_name": 2}[m]
            for m in methods
        )
        candidates[met.id] = {
            "id": met.id,
            "name": met.name,
            "formula": met.formula or "",
            "charge": met.charge,
            "compartment": met.compartment,
            "pubchem_cids": "|".join(sorted(model_cids)),
            "inchi": model_inchi,
            "methods": "|".join(sorted(methods)),
            "priority": priority,
        }
    return sorted(
        candidates.values(),
        key=lambda x: (x["priority"], 0 if x["compartment"] == "e" else 1, x["id"]),
    )


def direct_reactions(
    model: cobra.Model, substrate_ids: set[str], product_ids: set[str]
) -> list[cobra.Reaction]:
    hits: list[cobra.Reaction] = []
    for reaction in model.reactions:
        if reaction.boundary:
            continue
        sub_coeffs = [
            reaction.metabolites.get(model.metabolites.get_by_id(mid), 0.0)
            for mid in substrate_ids
        ]
        prod_coeffs = [
            reaction.metabolites.get(model.metabolites.get_by_id(mid), 0.0)
            for mid in product_ids
        ]
        forward = reaction.upper_bound > TOL and any(c < 0 for c in sub_coeffs) and any(
            c > 0 for c in prod_coeffs
        )
        reverse = reaction.lower_bound < -TOL and any(c > 0 for c in sub_coeffs) and any(
            c < 0 for c in prod_coeffs
        )
        if forward or reverse:
            hits.append(reaction)
    return hits


def build_reaction_graph(model: cobra.Model) -> nx.DiGraph:
    graph = nx.DiGraph()
    for reaction in model.reactions:
        if reaction.boundary:
            continue
        reactants = [met for met, coef in reaction.metabolites.items() if coef < 0]
        products = [met for met, coef in reaction.metabolites.items() if coef > 0]
        if reaction.upper_bound > TOL:
            for met in reactants:
                if base_metabolite_id(met.id) not in CURRENCY_BASE_IDS:
                    graph.add_edge(f"m:{met.id}", f"r:{reaction.id}")
            for met in products:
                if base_metabolite_id(met.id) not in CURRENCY_BASE_IDS:
                    graph.add_edge(f"r:{reaction.id}", f"m:{met.id}")
        if reaction.lower_bound < -TOL:
            for met in products:
                if base_metabolite_id(met.id) not in CURRENCY_BASE_IDS:
                    graph.add_edge(f"m:{met.id}", f"r:{reaction.id}")
            for met in reactants:
                if base_metabolite_id(met.id) not in CURRENCY_BASE_IDS:
                    graph.add_edge(f"r:{reaction.id}", f"m:{met.id}")
    return graph


def shortest_reaction_path(
    graph: nx.DiGraph, substrate_ids: list[str], product_ids: list[str]
) -> tuple[list[str], list[str]]:
    best_nodes: list[str] = []
    for substrate in substrate_ids:
        source = f"m:{substrate}"
        if source not in graph:
            continue
        for product in product_ids:
            target = f"m:{product}"
            if target not in graph:
                continue
            try:
                path = nx.shortest_path(graph, source=source, target=target)
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                continue
            if not best_nodes or len(path) < len(best_nodes):
                best_nodes = path
    reactions = [node[2:] for node in best_nodes if node.startswith("r:")]
    metabolites = [node[2:] for node in best_nodes if node.startswith("m:")]
    return reactions, metabolites


def exchange_for_metabolite(model: cobra.Model, met: cobra.Metabolite) -> cobra.Reaction | None:
    for reaction in model.exchanges:
        if met in reaction.metabolites:
            return reaction
    return None


def configure_substrate_only_carbon(
    model: cobra.Model, substrate_boundary: cobra.Reaction, allow_oxygen: bool
) -> None:
    for exchange in model.exchanges:
        if exchange.lower_bound >= 0:
            continue
        met = next(iter(exchange.metabolites), None)
        formula = parse_formula(met.formula if met is not None else "")
        unknown = not formula
        contains_carbon = formula.get("C", 0) > 0
        is_oxygen = met is not None and base_metabolite_id(met.id).lower() in {"o2", "oxygen"}
        if contains_carbon or unknown or (is_oxygen and not allow_oxygen):
            exchange.lower_bound = 0.0
    substrate_boundary.lower_bound = -10.0
    if substrate_boundary.upper_bound < 0:
        substrate_boundary.upper_bound = 0.0


def flux_test(
    model: cobra.Model,
    substrate_met: cobra.Metabolite,
    product_met: cobra.Metabolite,
    allow_oxygen: bool,
) -> dict[str, object]:
    original_objective = model.objective
    substrate_exchange = exchange_for_metabolite(model, substrate_met)
    product_exchange = exchange_for_metabolite(model, product_met)
    exchange_ready = substrate_exchange is not None and product_exchange is not None
    with model:
        if substrate_exchange is None:
            substrate_boundary = model.add_boundary(
                substrate_met, type="sink", reaction_id=f"SK_FMH_SUB_{substrate_met.id}"
            )
        else:
            substrate_boundary = substrate_exchange
        if product_exchange is None:
            product_boundary = model.add_boundary(
                product_met, type="demand", reaction_id=f"DM_FMH_PROD_{product_met.id}"
            )
        else:
            product_boundary = product_exchange

        configure_substrate_only_carbon(model, substrate_boundary, allow_oxygen)
        product_boundary.lower_bound = 0.0
        product_boundary.upper_bound = max(1000.0, product_boundary.upper_bound)
        model.objective = product_boundary
        model.objective_direction = "max"

        substrate_boundary.lower_bound = 0.0
        without = model.optimize(raise_error=False)
        without_value = (
            float(without.objective_value)
            if without.status == "optimal" and without.objective_value is not None
            else math.nan
        )

        substrate_boundary.lower_bound = -10.0
        with_substrate = model.optimize(raise_error=False)
        with_value = (
            float(with_substrate.objective_value)
            if with_substrate.status == "optimal" and with_substrate.objective_value is not None
            else math.nan
        )
        substrate_flux = (
            float(with_substrate.fluxes.get(substrate_boundary.id, math.nan))
            if with_substrate.status == "optimal"
            else math.nan
        )

        growth_value = math.nan
        fva_min = math.nan
        fva_max = math.nan
        try:
            model.objective = original_objective
            model.objective_direction = "max"
            growth = model.optimize(raise_error=False)
            if growth.status == "optimal" and growth.objective_value is not None:
                growth_value = float(growth.objective_value)
                if exchange_ready and growth_value > TOL:
                    fva = flux_variability_analysis(
                        model,
                        reaction_list=[product_boundary],
                        fraction_of_optimum=0.9,
                        loopless=False,
                        processes=1,
                    )
                    fva_min = float(fva.loc[product_boundary.id, "minimum"])
                    fva_max = float(fva.loc[product_boundary.id, "maximum"])
        except Exception:
            pass

    delta = (
        with_value - without_value
        if math.isfinite(with_value) and math.isfinite(without_value)
        else math.nan
    )
    return {
        "exchange_ready": exchange_ready,
        "substrate_boundary_id": substrate_boundary.id,
        "product_boundary_id": product_boundary.id,
        "product_max_with_substrate": with_value,
        "product_max_without_substrate": without_value,
        "substrate_attributable_product_delta": delta,
        "substrate_boundary_flux_at_product_max": substrate_flux,
        "substrate_dependent_flux_feasible": bool(math.isfinite(delta) and delta > TOL),
        "growth_max_substrate_only": growth_value,
        "product_fva_min_at_90pct_growth": fva_min,
        "product_fva_max_at_90pct_growth": fva_max,
    }


def serialize_candidates(candidates: list[dict[str, Any]]) -> dict[str, str]:
    if not candidates:
        return {
            "ids": "",
            "names": "",
            "formulas": "",
            "charges": "",
            "compartments": "",
            "pubchem_cids": "",
            "inchis": "",
            "methods": "",
        }
    return {
        "ids": "|".join(str(c["id"]) for c in candidates),
        "names": "|".join(str(c["name"]) for c in candidates),
        "formulas": "|".join(str(c["formula"]) for c in candidates),
        "charges": "|".join(str(c["charge"]) for c in candidates),
        "compartments": "|".join(str(c["compartment"]) for c in candidates),
        "pubchem_cids": "||".join(str(c["pubchem_cids"]) for c in candidates),
        "inchis": "||".join(str(c["inchi"]) for c in candidates),
        "methods": "||".join(str(c["methods"]) for c in candidates),
    }


def read_agora2_model(path: Path) -> tuple[cobra.Model, str]:
    """Read an AGORA2 model, repairing only invalid legacy text encoding."""
    try:
        return read_sbml_model(str(path)), "native_utf8"
    except Exception as original_exc:
        raw = path.read_bytes()
        try:
            raw.decode("utf-8")
        except UnicodeDecodeError:
            text = raw.decode("cp1252", errors="replace")
            text = "".join(ch for ch in text if ch in "\t\n\r" or ord(ch) >= 0x20)
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".xml", encoding="utf-8", delete=True
            ) as handle:
                handle.write(text)
                handle.flush()
                return read_sbml_model(handle.name), "encoding_sanitized_copy"
        raise original_exc


def process_model(task: tuple[str, list[dict[str, Any]], dict[int, dict[str, Any]]]) -> dict[str, Any]:
    filename, route_records, compounds = task
    path = SBML_DIR / filename
    started = time.time()
    try:
        model, load_mode = read_agora2_model(path)
    except Exception as exc:
        return {
            "model_file": filename,
            "model_qc": {
                "agora2_model_file": filename,
                "load_status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "runtime_seconds": time.time() - started,
            },
            "rows": [
                {
                    **route,
                    "agora2_model_file": filename,
                    "model_load_status": "error",
                    "model_error": f"{type(exc).__name__}: {exc}",
                }
                for route in route_records
            ],
        }

    graph = build_reaction_graph(model)
    default_growth = model.slim_optimize(error_value=math.nan)
    model_qc = {
        "agora2_model_file": filename,
        "agora2_model_id": model.id,
        "load_status": "loaded",
        "load_mode": load_mode,
        "reaction_count": len(model.reactions),
        "metabolite_count": len(model.metabolites),
        "gene_count": len(model.genes),
        "exchange_count": len(model.exchanges),
        "default_objective_value": float(default_growth) if default_growth is not None else math.nan,
        "runtime_seconds": 0.0,
        "error": "",
    }
    rows: list[dict[str, Any]] = []
    for route in route_records:
        substrate_cid = int(route["substrate_pubchem_cid"])
        product_cid = (
            int(route["metabolite_pubchem_cid"])
            if pd.notna(route["metabolite_pubchem_cid"])
            else None
        )
        substrate_reference = compounds[substrate_cid]
        product_reference = compounds.get(product_cid, {}) if product_cid is not None else {}
        substrate_candidates = mapping_candidates(
            model,
            substrate_cid,
            route["substrate_name"],
            substrate_reference.get("pubchem_inchi", ""),
        )
        product_candidates = (
            mapping_candidates(
                model,
                product_cid,
                route["metabolite_name"],
                product_reference.get("pubchem_inchi", ""),
            )
            if product_cid is not None
            else []
        )
        substrate_serialized = serialize_candidates(substrate_candidates)
        product_serialized = serialize_candidates(product_candidates)
        selected_substrate = (
            model.metabolites.get_by_id(substrate_candidates[0]["id"])
            if substrate_candidates
            else None
        )
        selected_product = (
            model.metabolites.get_by_id(product_candidates[0]["id"])
            if product_candidates
            else None
        )

        direct: list[cobra.Reaction] = []
        path_reactions: list[str] = []
        path_metabolites: list[str] = []
        anaerobic: dict[str, object] = {}
        oxygen: dict[str, object] = {}
        if substrate_candidates and product_candidates:
            direct = direct_reactions(
                model,
                {c["id"] for c in substrate_candidates},
                {c["id"] for c in product_candidates},
            )
            path_reactions, path_metabolites = shortest_reaction_path(
                graph,
                [c["id"] for c in substrate_candidates],
                [c["id"] for c in product_candidates],
            )
            try:
                anaerobic = flux_test(model, selected_substrate, selected_product, False)
            except Exception as exc:
                anaerobic = {"flux_error": f"{type(exc).__name__}: {exc}"}
            try:
                oxygen = flux_test(model, selected_substrate, selected_product, True)
            except Exception as exc:
                oxygen = {"flux_error": f"{type(exc).__name__}: {exc}"}

        rows.append(
            {
                **route,
                "agora2_model_file": filename,
                "agora2_model_id": model.id,
                "model_load_status": "loaded",
                "model_error": "",
                "substrate_mapping_candidate_count": len(substrate_candidates),
                "substrate_metabolite_ids": substrate_serialized["ids"],
                "substrate_model_names": substrate_serialized["names"],
                "substrate_model_formulas": substrate_serialized["formulas"],
                "substrate_model_charges": substrate_serialized["charges"],
                "substrate_model_compartments": substrate_serialized["compartments"],
                "substrate_model_pubchem_cids": substrate_serialized["pubchem_cids"],
                "substrate_model_inchis": substrate_serialized["inchis"],
                "substrate_mapping_methods": substrate_serialized["methods"],
                "product_mapping_candidate_count": len(product_candidates),
                "product_metabolite_ids": product_serialized["ids"],
                "product_model_names": product_serialized["names"],
                "product_model_formulas": product_serialized["formulas"],
                "product_model_charges": product_serialized["charges"],
                "product_model_compartments": product_serialized["compartments"],
                "product_model_pubchem_cids": product_serialized["pubchem_cids"],
                "product_model_inchis": product_serialized["inchis"],
                "product_mapping_methods": product_serialized["methods"],
                "direct_reaction_count": len(direct),
                "direct_reaction_ids": "|".join(r.id for r in direct),
                "direct_reaction_names": "|".join(r.name for r in direct),
                "direct_reaction_equations": "||".join(r.reaction for r in direct),
                "direct_reaction_ec_codes": "||".join(
                    "|".join(as_tokens(r.annotation.get("ec-code"))) for r in direct
                ),
                "topological_path_exists": bool(path_reactions),
                "topological_path_reaction_count": len(path_reactions),
                "topological_path_reaction_ids": "|".join(path_reactions),
                "topological_path_metabolite_ids": "|".join(path_metabolites),
                **{f"anaerobic_{k}": v for k, v in anaerobic.items()},
                **{f"oxygen_permissive_{k}": v for k, v in oxygen.items()},
            }
        )

    model_qc["runtime_seconds"] = time.time() - started
    return {"model_file": filename, "model_qc": model_qc, "rows": rows}


def unmapped_rows(
    routes: pd.DataFrame, taxon_summary: pd.DataFrame, taxon_edges: pd.DataFrame
) -> list[dict[str, Any]]:
    mapped_taxa = set(
        taxon_edges.loc[taxon_edges["agora2_model_file"].astype(str).str.len().gt(0), "microbe_name"]
    )
    status = taxon_summary.set_index("microbe_name").to_dict(orient="index")
    rows = []
    for route in routes.loc[~routes["microbe_name"].isin(mapped_taxa)].to_dict(orient="records"):
        meta = status[route["microbe_name"]]
        rows.append(
            {
                **route,
                "agora2_model_file": "",
                "agora2_model_id": "",
                "model_match_level": meta["model_match_level"],
                "model_match_basis": meta["model_match_basis"],
                "taxonomy_analysis_eligible": False,
                "model_load_status": "not_testable_no_model",
                "model_error": "",
                "substrate_mapping_candidate_count": 0,
                "product_mapping_candidate_count": 0,
                "direct_reaction_count": 0,
                "topological_path_exists": False,
            }
        )
    return rows


def main() -> None:
    global CHAIN, T2D_LINKS, PUBCHEM, SBML_DIR, BRIDGE
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain", type=Path, required=True, help="Checksum-identified two-hop route table")
    parser.add_argument("--t2d-links", type=Path, required=True, help="T2D target-link table")
    parser.add_argument("--pubchem", type=Path, required=True, help="Compound structure table")
    parser.add_argument("--models", type=Path, required=True, help="AGORA2 SBML directory")
    parser.add_argument("--bridge", type=Path, required=True, help="Taxon/model bridge table")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-models", type=int, default=0)
    args = parser.parse_args()
    CHAIN, T2D_LINKS, PUBCHEM = args.chain, args.t2d_links, args.pubchem
    SBML_DIR, BRIDGE = args.models, args.bridge
    args.output.mkdir(parents=True, exist_ok=True)

    for required in [CHAIN, T2D_LINKS, PUBCHEM, BRIDGE, SBML_DIR]:
        if not required.exists():
            raise FileNotFoundError(required)

    chain = pd.read_csv(CHAIN, low_memory=False)
    if len(chain) != 8063 or chain["microbe_name"].nunique() != 147:
        raise ValueError("Expected frozen 8,063-row / 147-taxon gutMGene chain")
    compounds_df = pd.read_csv(PUBCHEM, low_memory=False)
    compounds_df["cid"] = pd.to_numeric(compounds_df["cid"], errors="raise").astype(int)
    if compounds_df["cid"].nunique() != 90:
        raise ValueError("Expected 90 unique PubChem structures")
    compounds = compounds_df.set_index("cid").to_dict(orient="index")

    taxon_summary, taxon_edges = build_taxon_mapping(chain)
    routes = build_routes(chain)
    taxon_summary.to_csv(args.output / "01_taxon_model_mapping_147.csv", index=False)
    taxon_edges.to_csv(args.output / "02_taxon_model_edges.csv", index=False)
    routes.to_csv(args.output / "03_unique_taxon_routes.csv", index=False)

    mapped_edges = taxon_edges.loc[
        taxon_edges["agora2_model_file"].astype(str).str.len().gt(0)
    ].copy()
    expanded = routes.merge(
        mapped_edges[
            [
                "microbe_name",
                "model_match_level",
                "model_match_basis",
                "analysis_eligible",
                "agora2_model_file",
            ]
        ],
        on="microbe_name",
        how="inner",
        validate="many_to_many",
    ).rename(columns={"analysis_eligible": "taxonomy_analysis_eligible"})
    if args.max_models:
        selected = sorted(expanded["agora2_model_file"].unique())[: args.max_models]
        expanded = expanded.loc[expanded["agora2_model_file"].isin(selected)].copy()

    by_model: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in expanded.to_dict(orient="records"):
        by_model[row["agora2_model_file"]].append(row)
    tasks = [(filename, records, compounds) for filename, records in sorted(by_model.items())]
    print(
        f"Processing {len(tasks)} AGORA2 models and {len(expanded)} taxon–strain–route rows "
        f"with {args.workers} workers",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    started = time.time()
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(process_model, task): task[0] for task in tasks}
        for done, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            rows.extend(result["rows"])
            qc_rows.append(result["model_qc"])
            if done % 50 == 0 or done == len(tasks):
                elapsed = time.time() - started
                print(
                    f"Completed {done}/{len(tasks)} models; {len(rows)} rows; "
                    f"elapsed {elapsed / 60:.1f} min",
                    flush=True,
                )

    if not args.max_models:
        rows.extend(unmapped_rows(routes, taxon_summary, taxon_edges))
    matrix = pd.DataFrame(rows)
    matrix.insert(0, "matrix_row_id", [f"MX{i:06d}" for i in range(1, len(matrix) + 1)])
    matrix.to_csv(args.output / "04_route_strain_flux_matrix_raw.csv", index=False)
    model_qc = pd.DataFrame(qc_rows).sort_values("agora2_model_file")
    model_qc.to_csv(args.output / "05_model_load_qc.csv", index=False)

    manifest = {
        "created_at_utc": pd.Timestamp.utcnow().isoformat(),
        "python": platform.python_version(),
        "cobra": cobra.__version__,
        "networkx": nx.__version__,
        "pandas": pd.__version__,
        "input_sha256": {
            str(CHAIN): sha256(CHAIN),
            str(T2D_LINKS): sha256(T2D_LINKS),
            str(PUBCHEM): sha256(PUBCHEM),
            str(BRIDGE): sha256(BRIDGE),
        },
        "agora2_model_directory": str(SBML_DIR),
        "agora2_xml_count": len(list(SBML_DIR.glob("*.xml"))),
        "chain_rows": len(chain),
        "taxa": chain["microbe_name"].nunique(),
        "unique_cids": compounds_df["cid"].nunique(),
        "unique_taxon_routes": len(routes),
        "mapped_taxa": int(taxon_summary["analysis_eligible"].sum()),
        "unmapped_taxa": int((~taxon_summary["analysis_eligible"]).sum()),
        "mapped_model_files": int(mapped_edges["agora2_model_file"].nunique()),
        "matrix_rows": len(matrix),
        "media_definition": (
            "all carbon-containing and unknown-formula exchange uptake closed; "
            "target substrate uptake capped at 10; non-carbon exchanges retain model bounds; "
            "primary anaerobic and oxygen-permissive sensitivity scenarios"
        ),
        "limitations": [
            "topological paths ignore multi-reactant AND logic",
            "flux feasibility is conditional on an artificial substrate-only-carbon medium",
            "missing reactions were not gap-filled",
            "host target genes were not interpreted as microbial genes",
        ],
    }
    (args.output / "00_run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
