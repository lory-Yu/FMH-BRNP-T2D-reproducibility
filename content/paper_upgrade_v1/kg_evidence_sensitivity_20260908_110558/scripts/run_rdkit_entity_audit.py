#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors, rdMolHash
from rdkit.Chem.MolStandardize import rdMolStandardize


RDLogger.DisableLog("rdApp.*")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cid_text(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64").astype(str).replace("<NA>", "")


def fetch_pubchem(cids: list[str], pubchem_raw: Path) -> dict[str, dict]:
    pubchem_raw.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict] = {}
    props = "Title,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight,Charge"
    for i in range(0, len(cids), 40):
        batch = cids[i : i + 40]
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(batch)}/property/{props}/JSON"
        response = requests.get(url, timeout=120, headers={"User-Agent": "FMH-entity-audit/1.0"})
        response.raise_for_status()
        path = pubchem_raw / f"missing_cid_properties_{i//40+1:02d}.json"
        path.write_bytes(response.content)
        for rec in response.json().get("PropertyTable", {}).get("Properties", []):
            records[str(rec["CID"])] = rec
    return records


def analyze_structure(smiles: str) -> dict:
    result = {
        "parse_ok": False, "sanitization_status": "missing", "rdkit_exact_isomeric_smiles": "",
        "rdkit_exact_inchikey": "", "connectivity_smiles": "", "parent_isomeric_smiles": "",
        "parent_inchikey": "", "parent_connectivity_smiles": "", "tautomer_parent_smiles": "",
        "n_fragments": None, "formal_charge": None, "stereo_centres_specified": None,
        "stereo_centres_unassigned": None, "molecular_formula": "", "exact_mass": None,
        "standardization_error": "",
    }
    if not smiles or str(smiles).lower() == "nan":
        return result
    try:
        mol = Chem.MolFromSmiles(str(smiles), sanitize=True)
        if mol is None:
            result["sanitization_status"] = "parse_failed"
            return result
        result["parse_ok"] = True
        result["sanitization_status"] = "sanitized"
        result["rdkit_exact_isomeric_smiles"] = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        result["rdkit_exact_inchikey"] = Chem.MolToInchiKey(mol)
        no_stereo = Chem.Mol(mol)
        Chem.RemoveStereochemistry(no_stereo)
        result["connectivity_smiles"] = Chem.MolToSmiles(no_stereo, canonical=True, isomericSmiles=False)
        result["n_fragments"] = len(Chem.GetMolFrags(mol))
        result["formal_charge"] = Chem.GetFormalCharge(mol)
        centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        result["stereo_centres_specified"] = sum(tag != "?" for _, tag in centres)
        result["stereo_centres_unassigned"] = sum(tag == "?" for _, tag in centres)
        result["molecular_formula"] = rdMolDescriptors.CalcMolFormula(mol)
        result["exact_mass"] = Descriptors.ExactMolWt(mol)

        parent = rdMolStandardize.FragmentParent(mol)
        parent = rdMolStandardize.Uncharger().uncharge(parent)
        result["parent_isomeric_smiles"] = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        result["parent_inchikey"] = Chem.MolToInchiKey(parent)
        parent_no_stereo = Chem.Mol(parent)
        Chem.RemoveStereochemistry(parent_no_stereo)
        result["parent_connectivity_smiles"] = Chem.MolToSmiles(parent_no_stereo, canonical=True, isomericSmiles=False)
        # MolHash is substantially faster and more stable for an audit key than
        # enumerating every possible tautomer for >10,000 large natural products.
        result["tautomer_parent_smiles"] = rdMolHash.MolHash(
            parent, rdMolHash.HashFunction.HetAtomTautomer
        )
    except Exception as exc:
        result["standardization_error"] = f"{type(exc).__name__}: {exc}"
        if result["parse_ok"]:
            result["sanitization_status"] = "parsed_standardization_failed"
    return result


def joined_unique(values: pd.Series, limit: int = 50) -> str:
    vals = sorted({str(x).strip() for x in values.dropna() if str(x).strip() and str(x).lower() != "nan"})
    return "; ".join(vals[:limit]) + (f"; ...(+{len(vals)-limit})" if len(vals) > limit else "")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit chemical identities used by the metabolism-aware route analysis.")
    parser.add_argument("--occ", type=Path, required=True, help="Herb-compound occurrence TSV (local derived input).")
    parser.add_argument("--routes", type=Path, required=True, help="Audited transformation-route TSV.")
    parser.add_argument("--herb-chem", type=Path, required=True, help="Locally obtained herb-compound integration table; not redistributed.")
    parser.add_argument("--route-structures", type=Path, required=True, help="Previously audited route-structure CSV.")
    parser.add_argument("--pubchem-raw", type=Path, required=True, help="Directory for a local PubChem response cache.")
    parser.add_argument("--output", type=Path, required=True, help="Output root; 01_chemistry and logs are created below it.")
    parser.add_argument("--fetch-missing", action="store_true", help="Retrieve missing public PubChem records through PUG REST.")
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    occ_path = args.occ.resolve()
    routes_path = args.routes.resolve()
    herb_chem_path = args.herb_chem.resolve()
    route_structures_path = args.route_structures.resolve()
    pubchem_raw = args.pubchem_raw.resolve()
    out = args.output.resolve()
    chemistry_out = out / "01_chemistry"
    logs_out = out / "logs"
    chemistry_out.mkdir(parents=True, exist_ok=True)
    logs_out.mkdir(parents=True, exist_ok=True)

    occ = pd.read_csv(occ_path, sep="\t", dtype={"parent_cid": str}, low_memory=False)
    official_cids = set(occ["parent_cid"].astype(str))
    routes = pd.read_csv(routes_path, sep="\t", dtype={"parent_cid": str, "product_cid": str})
    route_cids = set(routes["parent_cid"].dropna().astype(str)) | set(routes["product_cid"].dropna().astype(str))
    target_cids = official_cids | route_cids
    hc = pd.read_csv(herb_chem_path, low_memory=False)
    hc = hc[hc["fmh_official"] == True].copy()  # noqa: E712
    hc["cid"] = cid_text(hc["pubchem_cid"])
    hc = hc[hc["cid"].isin(target_cids)].copy()
    hc["smiles"] = hc["smiles"].fillna("").astype(str).str.strip()
    hc["inchikey"] = hc["inchikey"].fillna("").astype(str).str.strip()

    source_rows = hc[["cid", "smiles", "inchikey", "compound_name", "herb_canonical", "source_db", "lotus_verified"]].copy()
    route_structures = pd.read_csv(route_structures_path, dtype={"cid": str}, low_memory=False)
    route_structures["cid"] = cid_text(route_structures["cid"])
    route_structures["smiles"] = route_structures["isomeric_smiles"].fillna(route_structures["canonical_smiles"]).fillna("")
    route_supplement = pd.DataFrame({
        "cid": route_structures["cid"], "smiles": route_structures["smiles"],
        "inchikey": route_structures["pubchem_inchikey"].fillna(""),
        "compound_name": route_structures["pubchem_title"].fillna(route_structures["gutmgene_name"]).fillna(""),
        "herb_canonical": "", "source_db": "previously audited PubChem route structure",
        "lotus_verified": False,
    })
    route_supplement = route_supplement[route_supplement["cid"].isin(target_cids)]
    source_rows = pd.concat([source_rows, route_supplement], ignore_index=True)
    source_rows["smiles"] = source_rows["smiles"].fillna("").astype(str).str.strip()
    present = set(source_rows.loc[source_rows["smiles"].ne(""), "cid"])
    missing = sorted(target_cids - present, key=int)
    pd.DataFrame({"pubchem_cid": missing}).to_csv(chemistry_out / "00_missing_structure_cids_before_pubchem.tsv", sep="\t", index=False)
    pubchem: dict[str, dict] = {}
    if args.fetch_missing and missing:
        pubchem = fetch_pubchem(missing, pubchem_raw)
    else:
        for path in sorted(pubchem_raw.glob("missing_cid_properties_*.json")):
            for rec in json.loads(path.read_text()).get("PropertyTable", {}).get("Properties", []):
                pubchem[str(rec["CID"])] = rec

    fetched_rows = []
    for cid in missing:
        rec = pubchem.get(cid, {})
        smiles = rec.get("SMILES", rec.get("IsomericSMILES", rec.get("ConnectivitySMILES", rec.get("CanonicalSMILES", ""))))
        fetched_rows.append({"cid": cid, "smiles": smiles or "", "inchikey": rec.get("InChIKey", ""), "compound_name": rec.get("Title", ""), "herb_canonical": "", "source_db": "PubChem PUG REST missing-structure fill", "lotus_verified": False})
    source_rows = pd.concat([source_rows, pd.DataFrame(fetched_rows)], ignore_index=True)
    source_rows = source_rows[source_rows["smiles"].ne("")].copy()
    variants = source_rows.drop_duplicates(["cid", "smiles"]).reset_index(drop=True)
    analyses = pd.DataFrame([analyze_structure(x) for x in variants["smiles"]])
    audit = pd.concat([variants, analyses], axis=1)
    audit.to_csv(chemistry_out / "01_structure_variant_audit.tsv", sep="\t", index=False)

    summaries = []
    for cid in sorted(target_cids, key=int):
        g = audit[audit["cid"] == cid]
        base = source_rows[source_rows["cid"] == cid]
        valid = g[g["parse_ok"] == True]  # noqa: E712
        exact_keys = set(valid["rdkit_exact_inchikey"].dropna()) - {""}
        conn = set(valid["connectivity_smiles"].dropna()) - {""}
        parent = set(valid["parent_inchikey"].dropna()) - {""}
        taut = set(valid["tautomer_parent_smiles"].dropna()) - {""}
        supplied_keys = set(base["inchikey"].dropna()) - {""}
        if len(g) == 0:
            status = "unresolved_no_structure"
        elif len(valid) == 0:
            status = "unresolved_parse_failure"
        elif len(conn) > 1:
            status = "conflict_within_cid_connectivity"
        elif len(exact_keys) > 1:
            status = "stereo_or_protonation_variants_within_cid"
        else:
            status = "resolved_consistent"
        summaries.append({
            "pubchem_cid": cid,
            "compound_names": joined_unique(base["compound_name"]),
            "n_herbs": int(base.loc[base["herb_canonical"].ne(""), "herb_canonical"].nunique()),
            "n_source_rows": int(len(base)),
            "n_raw_smiles": int(g["smiles"].nunique()),
            "n_parseable_variants": int(len(valid)),
            "n_exact_structure_keys": len(exact_keys),
            "n_connectivity_keys": len(conn),
            "n_parent_keys": len(parent),
            "n_tautomer_parent_keys": len(taut),
            "n_supplied_inchikeys": len(supplied_keys),
            "has_multifragment_variant": bool(valid["n_fragments"].fillna(0).gt(1).any()),
            "has_charged_variant": bool(valid["formal_charge"].fillna(0).ne(0).any()),
            "has_unassigned_stereocentre": bool(valid["stereo_centres_unassigned"].fillna(0).gt(0).any()),
            "representative_exact_inchikey": sorted(exact_keys)[0] if len(exact_keys) == 1 else "",
            "representative_connectivity_smiles": sorted(conn)[0] if len(conn) == 1 else "",
            "representative_parent_inchikey": sorted(parent)[0] if len(parent) == 1 else "",
            "representative_tautomer_parent_smiles": sorted(taut)[0] if len(taut) == 1 else "",
            "entity_audit_status": status,
            "usable_for_strict_route_mapping": status in {"resolved_consistent", "stereo_or_protonation_variants_within_cid"},
            "interpretation": "structure identity audit only; not evidence of herbal occurrence or biological activity",
        })
    consensus = pd.DataFrame(summaries)

    for col, label in [
        ("representative_exact_inchikey", "exact_structure"),
        ("representative_parent_inchikey", "standardized_parent"),
        ("representative_tautomer_parent_smiles", "tautomer_parent"),
    ]:
        counts = consensus.loc[consensus[col].ne("")].groupby(col)["pubchem_cid"].agg(list)
        duplicate_map = {key: vals for key, vals in counts.items() if len(vals) > 1}
        consensus[f"cross_cid_{label}_group_size"] = consensus[col].map(lambda x: len(duplicate_map.get(x, [])) if x else 0)
        consensus[f"cross_cid_{label}_cids"] = consensus[col].map(lambda x: ";".join(duplicate_map.get(x, [])) if x else "")

    consensus.to_csv(chemistry_out / "02_cid_consensus_audit.tsv", sep="\t", index=False)
    for col, filename in [
        ("representative_exact_inchikey", "03_cross_cid_exact_duplicate_groups.tsv"),
        ("representative_parent_inchikey", "04_cross_cid_parent_duplicate_groups.tsv"),
        ("representative_tautomer_parent_smiles", "05_cross_cid_tautomer_parent_groups.tsv"),
    ]:
        grouped = consensus[consensus[col].ne("")].groupby(col).agg(n_cids=("pubchem_cid", "nunique"), cids=("pubchem_cid", lambda x: ";".join(sorted(x, key=int))), compound_names=("compound_names", joined_unique), n_herbs=("n_herbs", "sum")).reset_index()
        grouped[grouped["n_cids"] > 1].to_csv(chemistry_out / filename, sep="\t", index=False)

    herb_map = occ[["herb_canonical", "parent_cid"]].merge(consensus, left_on="parent_cid", right_on="pubchem_cid", how="left")
    herb_summary = herb_map.groupby("herb_canonical").agg(
        n_unique_cids=("parent_cid", "nunique"),
        n_strict_usable_cids=("usable_for_strict_route_mapping", "sum"),
        n_unresolved_no_structure=("entity_audit_status", lambda x: int((x == "unresolved_no_structure").sum())),
        n_parse_failures=("entity_audit_status", lambda x: int((x == "unresolved_parse_failure").sum())),
        n_within_cid_connectivity_conflicts=("entity_audit_status", lambda x: int((x == "conflict_within_cid_connectivity").sum())),
        n_stereo_or_protonation_variant_cids=("entity_audit_status", lambda x: int((x == "stereo_or_protonation_variants_within_cid").sum())),
        n_cross_cid_exact_duplicate_members=("cross_cid_exact_structure_group_size", lambda x: int((x > 1).sum())),
        n_cross_cid_parent_duplicate_members=("cross_cid_standardized_parent_group_size", lambda x: int((x > 1).sum())),
    ).reset_index()
    herb_summary["strict_usable_fraction"] = herb_summary["n_strict_usable_cids"] / herb_summary["n_unique_cids"]
    herb_summary.to_csv(chemistry_out / "06_herb_entity_qc_summary.tsv", sep="\t", index=False)

    parent_qc = consensus.add_prefix("parent_")
    product_qc = consensus.add_prefix("product_")
    route_qc = routes.merge(parent_qc, left_on="parent_cid", right_on="parent_pubchem_cid", how="left").merge(product_qc, left_on="product_cid", right_on="product_pubchem_cid", how="left")
    route_qc["both_entities_strict_usable"] = route_qc["parent_usable_for_strict_route_mapping"].fillna(False) & route_qc["product_usable_for_strict_route_mapping"].fillna(False)
    gins_pat = r"ginsenoside|compound\s*k|protopanaxadiol|\brb1\b|\brd\b|\bf2\b"
    route_qc["ginsenoside_route_text_match"] = route_qc[["substrate_name", "product_name"]].fillna("").agg(" ".join, axis=1).str.contains(gins_pat, case=False, regex=True)
    route_qc.to_csv(chemistry_out / "07_route_entity_crosswalk.tsv", sep="\t", index=False)

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "rdkit": rdBase.rdkitVersion,
        "official_materials": int(occ["herb_canonical"].nunique()),
        "herb_cid_occurrences": int(len(occ)),
        "unique_herb_parent_cids": int(len(official_cids)),
        "unique_route_parent_product_cids": int(len(route_cids)),
        "unique_cids": int(len(consensus)),
        "missing_structures_before_pubchem": len(missing),
        "pubchem_missing_records_retrieved": len(pubchem),
        "status_counts": dict(Counter(consensus["entity_audit_status"])),
        "exact_cross_cid_duplicate_groups": int((consensus.groupby("representative_exact_inchikey")["pubchem_cid"].nunique() > 1).sum()),
        "parent_cross_cid_duplicate_groups": int((consensus.groupby("representative_parent_inchikey")["pubchem_cid"].nunique() > 1).sum()),
        "route_rows": int(len(routes)),
        "primary_route_rows": int(routes["human_database_route_eligible"].sum()),
        "primary_routes_with_both_entities_strict_usable": int((route_qc["human_database_route_eligible"] & route_qc["both_entities_strict_usable"]).sum()),
        "ginsenoside_route_rows_in_v2": int(route_qc["ginsenoside_route_text_match"].sum()),
        "input_sha256": {
            "herb_compound_occurrences": sha256(occ_path),
            "transformation_routes": sha256(routes_path),
            "herb_compound_integration": sha256(herb_chem_path),
            "route_structures": sha256(route_structures_path),
        },
    }
    (logs_out / "rdkit_entity_audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(parse_args())
