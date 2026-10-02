#!/usr/bin/env python3
"""RDKit validation and evidence-tiering for the AGORA2 flux matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize


DEFAULT_INPUT = Path("04_route_strain_flux_matrix_raw.csv")
PUBCHEM = None
TAXON_MAP = None
MODEL_QC = None
RUN_MANIFEST = None
DEFAULT_OUT = Path("outputs")
TOL = 1e-8
RDLogger.DisableLog("rdApp.*")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def first_double_pipe(value: object) -> str:
    return clean(value).split("||", 1)[0]


def first_pipe(value: object) -> str:
    return clean(value).split("|", 1)[0]


def molecule_from(smiles: object, inchi: object) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(clean(smiles)) if clean(smiles) else None
    if mol is None and clean(inchi):
        mol = Chem.MolFromInchi(clean(inchi))
    return mol


def normalized_keys(mol: Chem.Mol) -> dict[str, str]:
    exact_key = Chem.MolToInchiKey(mol)
    uncharged = rdMolStandardize.Uncharger().uncharge(Chem.Mol(mol))
    parent = rdMolStandardize.FragmentParent(uncharged)
    parent_key = Chem.MolToInchiKey(parent)
    return {
        "exact_inchikey": exact_key,
        "connectivity_key": exact_key.split("-")[0],
        "parent_inchikey": parent_key,
        "parent_connectivity_key": parent_key.split("-")[0],
        "canonical_isomeric_smiles": Chem.MolToSmiles(parent, isomericSmiles=True),
        "canonical_connectivity_smiles": Chem.MolToSmiles(parent, isomericSmiles=False),
    }


def mapping_relation(
    target_smiles: object,
    target_inchi: object,
    model_inchi: object,
    methods: object,
) -> tuple[str, bool, str, str]:
    method_text = clean(methods)
    if not method_text:
        return "not_mapped", False, "", ""
    target = molecule_from(target_smiles, target_inchi)
    model = Chem.MolFromInchi(clean(model_inchi)) if clean(model_inchi) else None
    if target is not None and model is not None:
        target_keys = normalized_keys(target)
        model_keys = normalized_keys(model)
        if target_keys["exact_inchikey"] == model_keys["exact_inchikey"]:
            return (
                "exact_stereochemical_structure",
                True,
                target_keys["parent_connectivity_key"],
                model_keys["parent_connectivity_key"],
            )
        if target_keys["connectivity_key"] == model_keys["connectivity_key"]:
            return (
                "same_connectivity_protonation_or_stereo_diff",
                True,
                target_keys["parent_connectivity_key"],
                model_keys["parent_connectivity_key"],
            )
        if target_keys["parent_connectivity_key"] == model_keys["parent_connectivity_key"]:
            return (
                "same_parent_connectivity_after_standardization",
                True,
                target_keys["parent_connectivity_key"],
                model_keys["parent_connectivity_key"],
            )
        return (
            "structure_conflict",
            False,
            target_keys["parent_connectivity_key"],
            model_keys["parent_connectivity_key"],
        )
    if "exact_pubchem_cid" in method_text:
        return "exact_cid_annotation_structure_unavailable", True, "", ""
    return "name_only_structure_unavailable", False, "", ""


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(False, index=frame.index)
    values = frame[column]
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(math.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def classify_support(row: pd.Series) -> str:
    if clean(row.get("model_load_status")) == "not_testable_no_model":
        return "F7_no_AGORA2_strain_model"
    if clean(row.get("model_load_status")) != "loaded":
        return "F8_model_load_error"
    substrate_count = int(row.get("substrate_mapping_candidate_count", 0) or 0)
    product_count = int(row.get("product_mapping_candidate_count", 0) or 0)
    if substrate_count == 0 and product_count == 0:
        return "F6_both_compounds_absent"
    if substrate_count == 0:
        return "F6_substrate_absent"
    if product_count == 0:
        return "F6_product_absent"
    if not bool(row.get("substrate_structure_verified")) or not bool(
        row.get("product_structure_verified")
    ):
        return "F5_compound_mapping_unverified"
    if bool(row.get("anaerobic_substrate_dependent_flux_feasible")):
        if bool(row.get("anaerobic_exchange_ready")):
            return "F1_anaerobic_exchange_flux_supported"
        return "F2_anaerobic_artificial_sink_flux_only"
    if bool(row.get("oxygen_permissive_substrate_dependent_flux_feasible")):
        if bool(row.get("oxygen_permissive_exchange_ready")):
            return "F2_oxygen_permissive_exchange_flux_only"
        return "F3_oxygen_permissive_artificial_sink_only"
    if int(row.get("direct_reaction_count", 0) or 0) > 0:
        return "F3_direct_reaction_but_flux_blocked"
    if bool(row.get("topological_path_exists")):
        return "F4_topological_path_only"
    return "F4_compounds_present_no_supported_route"


def audit_pubchem(compounds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in compounds.itertuples(index=False):
        mol = molecule_from(row.isomeric_smiles, row.pubchem_inchi)
        if mol is None:
            rows.append({"cid": row.cid, "rdkit_parse_ok": False})
            continue
        keys = normalized_keys(mol)
        rows.append(
            {
                "cid": int(row.cid),
                "rdkit_parse_ok": True,
                "rdkit_isomeric_smiles": keys["canonical_isomeric_smiles"],
                "rdkit_connectivity_smiles": keys["canonical_connectivity_smiles"],
                "rdkit_inchikey": keys["exact_inchikey"],
                "rdkit_parent_connectivity_key": keys["parent_connectivity_key"],
                "rdkit_formula": rdMolDescriptors.CalcMolFormula(mol),
                "pubchem_inchikey_matches_rdkit": clean(row.pubchem_inchikey)
                == keys["exact_inchikey"],
            }
        )
    return compounds.merge(pd.DataFrame(rows), on="cid", how="left", validate="one_to_one")


def pick_column(frame: pd.DataFrame, source: str, splitter=first_pipe) -> pd.Series:
    if source not in frame:
        return pd.Series("", index=frame.index, dtype=object)
    return frame[source].map(splitter)


def mapping_columns(
    frame: pd.DataFrame, prefix: str, compounds_by_cid: dict[int, dict[str, object]]
) -> pd.DataFrame:
    result = pd.DataFrame(index=frame.index)
    result[f"{prefix}_selected_metabolite_id"] = pick_column(
        frame, f"{prefix}_metabolite_ids"
    )
    result[f"{prefix}_selected_model_name"] = pick_column(
        frame, f"{prefix}_model_names"
    )
    result[f"{prefix}_selected_model_formula"] = pick_column(
        frame, f"{prefix}_model_formulas"
    )
    result[f"{prefix}_selected_model_charge"] = pick_column(
        frame, f"{prefix}_model_charges"
    )
    result[f"{prefix}_selected_model_compartment"] = pick_column(
        frame, f"{prefix}_model_compartments"
    )
    result[f"{prefix}_selected_model_pubchem_cids"] = pick_column(
        frame, f"{prefix}_model_pubchem_cids", first_double_pipe
    )
    result[f"{prefix}_selected_model_inchi"] = pick_column(
        frame, f"{prefix}_model_inchis", first_double_pipe
    )
    result[f"{prefix}_selected_mapping_methods"] = pick_column(
        frame, f"{prefix}_mapping_methods", first_double_pipe
    )

    cid_col = "substrate_pubchem_cid" if prefix == "substrate" else "metabolite_pubchem_cid"
    relations = []
    for idx, row in frame.iterrows():
        raw_cid = row.get(cid_col)
        ref: dict[str, object] = {}
        if pd.notna(raw_cid):
            ref = compounds_by_cid.get(int(float(raw_cid)), {})
        relations.append(
            mapping_relation(
                ref.get("isomeric_smiles", ""),
                ref.get("pubchem_inchi", ""),
                result.at[idx, f"{prefix}_selected_model_inchi"],
                result.at[idx, f"{prefix}_selected_mapping_methods"],
            )
        )
    result[
        [
            f"{prefix}_structure_relation",
            f"{prefix}_structure_verified",
            f"{prefix}_reference_parent_connectivity_key",
            f"{prefix}_model_parent_connectivity_key",
        ]
    ] = pd.DataFrame(relations, index=frame.index)
    return result


def summarize_routes(frame: pd.DataFrame) -> pd.DataFrame:
    route_cols = [
        "taxon_route_id",
        "microbe_name",
        "substrate_pubchem_cid",
        "substrate_name",
        "metabolite_pubchem_cid",
        "metabolite_name",
        "target_genes",
        "t2d_core_target_genes",
        "in_t2d_core_75",
        "current_fmh_herbs",
        "raw_chain_rows",
    ]
    rows = []
    for _, group in frame.groupby("taxon_route_id", sort=True, dropna=False):
        base = {c: group.iloc[0].get(c, "") for c in route_cols}
        loaded = group["model_load_status"].eq("loaded")
        both = (
            numeric_series(group, "substrate_mapping_candidate_count").gt(0)
            & numeric_series(group, "product_mapping_candidate_count").gt(0)
        )
        formal = bool_series(group, "formal_flux_eligible")
        rows.append(
            {
                **base,
                "strain_models_considered": int(group["agora2_model_file"].replace("", pd.NA).nunique()),
                "strain_models_loaded": int(loaded.sum()),
                "strain_models_both_compounds_mapped": int((loaded & both).sum()),
                "strain_models_both_structures_verified": int(
                    (
                        loaded
                        & bool_series(group, "substrate_structure_verified")
                        & bool_series(group, "product_structure_verified")
                    ).sum()
                ),
                "strain_models_direct_reaction": int(
                    numeric_series(group, "direct_reaction_count").gt(0).sum()
                ),
                "strain_models_topological_path": int(
                    bool_series(group, "topological_path_exists").sum()
                ),
                "strain_models_formal_anaerobic_flux": int(formal.sum()),
                "any_formal_anaerobic_flux": bool(formal.any()),
                "formal_flux_fraction_loaded": (
                    float(formal.sum() / loaded.sum()) if loaded.sum() else math.nan
                ),
                "best_support_tier": sorted(group["cobrapy_support_tier"].dropna().astype(str))[0]
                if group["cobrapy_support_tier"].notna().any()
                else "",
            }
        )
    return pd.DataFrame(rows)


def build_crosswalk(frame: pd.DataFrame) -> pd.DataFrame:
    blocks = []
    for role, cid_col, name_col in [
        ("substrate", "substrate_pubchem_cid", "substrate_name"),
        ("product", "metabolite_pubchem_cid", "metabolite_name"),
    ]:
        block = pd.DataFrame(
            {
                "role": role,
                "gutmgene_pubchem_cid": frame[cid_col],
                "gutmgene_name": frame[name_col],
                "agora2_metabolite_id": frame[f"{role}_selected_metabolite_id"],
                "agora2_name": frame[f"{role}_selected_model_name"],
                "agora2_formula": frame[f"{role}_selected_model_formula"],
                "agora2_charge": frame[f"{role}_selected_model_charge"],
                "agora2_inchi": frame[f"{role}_selected_model_inchi"],
                "mapping_methods": frame[f"{role}_selected_mapping_methods"],
                "structure_relation": frame[f"{role}_structure_relation"],
                "structure_verified": frame[f"{role}_structure_verified"],
                "microbe_name": frame["microbe_name"],
                "agora2_model_file": frame["agora2_model_file"],
            }
        )
        blocks.append(block)
    long = pd.concat(blocks, ignore_index=True)
    long = long.loc[long["gutmgene_pubchem_cid"].notna()].copy()
    keys = [
        "role",
        "gutmgene_pubchem_cid",
        "gutmgene_name",
        "agora2_metabolite_id",
        "agora2_name",
        "agora2_formula",
        "agora2_charge",
        "agora2_inchi",
        "mapping_methods",
        "structure_relation",
        "structure_verified",
    ]
    return (
        long.groupby(keys, dropna=False, sort=True)
        .agg(
            taxon_count=("microbe_name", "nunique"),
            strain_model_count=("agora2_model_file", "nunique"),
            matrix_row_count=("microbe_name", "size"),
        )
        .reset_index()
    )


def build_compound_coverage(frame: pd.DataFrame, compounds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ref in compounds.itertuples(index=False):
        cid = int(ref.cid)
        sub = frame[pd.to_numeric(frame["substrate_pubchem_cid"], errors="coerce").eq(cid)]
        pro = frame[pd.to_numeric(frame["metabolite_pubchem_cid"], errors="coerce").eq(cid)]
        mapped = pd.concat(
            [
                sub.assign(_mapped=numeric_series(sub, "substrate_mapping_candidate_count").gt(0)),
                pro.assign(_mapped=numeric_series(pro, "product_mapping_candidate_count").gt(0)),
            ],
            ignore_index=True,
        )
        verified_count = int(bool_series(sub, "substrate_structure_verified").sum()) + int(
            bool_series(pro, "product_structure_verified").sum()
        )
        rows.append(
            {
                "cid": cid,
                "preferred_name": ref.gutmgene_name or ref.pubchem_title,
                "present_as_substrate": not sub.empty,
                "present_as_product": not pro.empty,
                "models_with_any_mapping": int(mapped.loc[mapped["_mapped"], "agora2_model_file"].nunique())
                if not mapped.empty
                else 0,
                "taxa_with_any_mapping": int(mapped.loc[mapped["_mapped"], "microbe_name"].nunique())
                if not mapped.empty
                else 0,
                "verified_mapping_rows": verified_count,
                "rdkit_parse_ok": bool(ref.rdkit_parse_ok),
                "pubchem_inchikey": ref.pubchem_inchikey,
                "molecular_formula": ref.molecular_formula,
            }
        )
    return pd.DataFrame(rows)


def write_report(out: Path, metrics: dict[str, object], checks: dict[str, object]) -> None:
    text = f"""# AGORA2–gutMGene–FMH COBRApy feasibility audit

## Scope

This analysis links the audited gutMGene taxon–substrate–product routes to local AGORA2 v2.01 strain models. It is a model-feasibility audit, not evidence that a reaction occurs in vivo. Host target genes remain downstream host annotations and are never treated as microbial enzyme genes.

## Main counts

- Source gutMGene chain rows: **{metrics['source_chain_rows']:,}**
- gutMGene taxa: **{metrics['taxa']:,}**; taxonomically eligible for an AGORA2 strain mapping: **{metrics['mapped_taxa']:,}**
- Unique taxon–substrate–product routes: **{metrics['taxon_routes']:,}**
- Relevant AGORA2 strain models: **{metrics['strain_models']:,}**; successfully parsed: **{metrics['models_loaded']:,}**
- PubChem compounds audited by RDKit: **{metrics['compounds']:,}**
- Taxon–strain–route matrix rows: **{metrics['matrix_rows']:,}**
- Rows with both compounds mapped and RDKit-verified: **{metrics['both_structures_verified']:,}**
- Formal anaerobic exchange-flux positive rows (F1): **{metrics['formal_flux_rows']:,}**
- Routes with at least one F1 strain: **{metrics['formal_flux_routes']:,}**
- T2D-core subset rows/routes with F1 support: **{metrics['t2d_formal_flux_rows']:,} / {metrics['t2d_formal_flux_routes']:,}**

## Evidence policy

F1 is the only tier admitted to the formal flux conclusion: taxonomy is eligible; both compounds are matched by exact annotation or RDKit structural equivalence; native substrate and product exchange reactions exist; and anaerobic maximum product export increases when substrate uptake is enabled. F2–F4 are sensitivity, direct-reaction, or topology-only evidence. F5–F8 are mapping/model failures or unavailable models.

No gap filling was performed. Artificial sink results are retained as stoichiometric-only sensitivity tests and excluded from F1. The medium is a deliberately artificial substrate-only-carbon environment. Therefore a positive result means the reconstruction permits the conversion under those constraints, while a negative result does not prove that the organism cannot perform it biologically.

## Validation

All mandatory checks passed: **{all(checks.values())}**. Machine-readable checks are in `13_validation_checks.json`.

## Manuscript-safe interpretation

The output supports a qualified statement that gutMGene-reported FMH-associated transformations can be stratified by taxon-to-strain mapping, chemical identity, reaction topology, and conditional strain-model flux feasibility. It does not validate all 8,063 composed gutMGene–host-target rows as independent microbial reactions, and it does not establish patient-level activity or treatment efficacy.
"""
    (out.parent / "reports").mkdir(parents=True, exist_ok=True)
    (out.parent / "reports" / "AGORA2_COBRAPY_FEASIBILITY_REPORT.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    global PUBCHEM, TAXON_MAP, MODEL_QC, RUN_MANIFEST
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--pubchem", type=Path, required=True, help="Compound structure table")
    parser.add_argument("--taxon-map", type=Path, required=True, help="Taxon/model mapping table")
    parser.add_argument("--model-qc", type=Path, required=True, help="Model load QC table")
    parser.add_argument("--manifest", type=Path, required=True, help="Flux run manifest")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    PUBCHEM, TAXON_MAP, MODEL_QC, RUN_MANIFEST = args.pubchem, args.taxon_map, args.model_qc, args.manifest
    args.output.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.input, low_memory=False)
    compounds_raw = pd.read_csv(PUBCHEM)
    compounds = audit_pubchem(compounds_raw)
    compounds.to_csv(args.output / "06_compound_structure_master_90.csv", index=False)
    by_cid = {
        int(row["cid"]): row for row in compounds.to_dict(orient="records")
    }

    substrate = mapping_columns(frame, "substrate", by_cid)
    product = mapping_columns(frame, "product", by_cid)
    frame = pd.concat([frame, substrate, product], axis=1)
    frame["taxonomy_analysis_eligible"] = bool_series(frame, "taxonomy_analysis_eligible")
    frame["formal_flux_eligible"] = (
        frame["taxonomy_analysis_eligible"]
        & bool_series(frame, "substrate_structure_verified")
        & bool_series(frame, "product_structure_verified")
        & bool_series(frame, "anaerobic_exchange_ready")
        & bool_series(frame, "anaerobic_substrate_dependent_flux_feasible")
    )
    frame["cobrapy_support_tier"] = frame.apply(classify_support, axis=1)
    frame.to_csv(args.output / "07_route_strain_flux_matrix_final.csv", index=False)

    crosswalk = build_crosswalk(frame)
    crosswalk.to_csv(args.output / "08_compound_agora2_crosswalk.csv", index=False)
    coverage = build_compound_coverage(frame, compounds)
    coverage.to_csv(args.output / "08b_compound_coverage_90.csv", index=False)

    routes = summarize_routes(frame)
    routes.to_csv(args.output / "09_taxon_route_summary_300.csv", index=False)

    taxon_map = pd.read_csv(TAXON_MAP)
    taxon_stats = (
        routes.groupby("microbe_name", sort=True)
        .agg(
            route_count=("taxon_route_id", "nunique"),
            t2d_core_route_count=("in_t2d_core_75", "sum"),
            routes_with_formal_flux=("any_formal_anaerobic_flux", "sum"),
            strain_route_rows=("strain_models_considered", "sum"),
            formal_flux_strain_rows=("strain_models_formal_anaerobic_flux", "sum"),
        )
        .reset_index()
    )
    taxon_summary = taxon_map.merge(taxon_stats, on="microbe_name", how="left", validate="one_to_one")
    for col in [
        "route_count",
        "t2d_core_route_count",
        "routes_with_formal_flux",
        "strain_route_rows",
        "formal_flux_strain_rows",
    ]:
        taxon_summary[col] = taxon_summary[col].fillna(0).astype(int)
    taxon_summary.to_csv(args.output / "10_taxon_flux_summary_147.csv", index=False)

    core = frame.loc[bool_series(frame, "in_t2d_core_75")].copy()
    core.to_csv(args.output / "11_t2d_core_route_strain_matrix.csv", index=False)

    qc = pd.read_csv(MODEL_QC)
    with open(RUN_MANIFEST, encoding="utf-8") as handle:
        run_manifest = json.load(handle)
    metrics = {
        "source_chain_rows": int(run_manifest["chain_rows"]),
        "taxa": int(frame["microbe_name"].nunique()),
        "mapped_taxa": int(taxon_map["analysis_eligible"].astype(bool).sum()),
        "taxon_routes": int(routes["taxon_route_id"].nunique()),
        "strain_models": int(frame["agora2_model_file"].replace("", pd.NA).nunique()),
        "models_loaded": int(qc["load_status"].eq("loaded").sum()),
        "models_encoding_sanitized": int(qc.get("load_mode", pd.Series(dtype=str)).eq("encoding_sanitized_copy").sum()),
        "compounds": int(compounds["cid"].nunique()),
        "matrix_rows": int(len(frame)),
        "both_structures_verified": int(
            (bool_series(frame, "substrate_structure_verified") & bool_series(frame, "product_structure_verified")).sum()
        ),
        "formal_flux_rows": int(bool_series(frame, "formal_flux_eligible").sum()),
        "formal_flux_routes": int(routes["any_formal_anaerobic_flux"].sum()),
        "t2d_formal_flux_rows": int(bool_series(core, "formal_flux_eligible").sum()),
        "t2d_formal_flux_routes": int(
            routes.loc[bool_series(routes, "in_t2d_core_75"), "any_formal_anaerobic_flux"].sum()
        ),
    }
    pd.DataFrame(
        [{"metric": key, "value": value} for key, value in metrics.items()]
    ).to_csv(args.output / "12_analysis_summary.csv", index=False)

    checks = {
        "source_chain_is_8063": metrics["source_chain_rows"] == 8063,
        "taxa_are_147": metrics["taxa"] == 147,
        "compound_union_is_90": metrics["compounds"] == 90,
        "all_pubchem_structures_parse": bool(compounds["rdkit_parse_ok"].all()),
        "taxon_routes_are_300": metrics["taxon_routes"] == 300,
        "taxon_summary_has_147_rows": len(taxon_summary) == 147,
        "all_relevant_models_parsed": metrics["models_loaded"] == metrics["strain_models"],
        "formal_flux_requires_verified_structures": bool(
            (~bool_series(frame, "formal_flux_eligible")
             | (bool_series(frame, "substrate_structure_verified")
                & bool_series(frame, "product_structure_verified"))).all()
        ),
        "formal_flux_requires_native_exchanges": bool(
            (~bool_series(frame, "formal_flux_eligible") | bool_series(frame, "anaerobic_exchange_ready")).all()
        ),
        "formal_flux_requires_taxonomy_mapping": bool(
            (~bool_series(frame, "formal_flux_eligible") | bool_series(frame, "taxonomy_analysis_eligible")).all()
        ),
        "final_rows_equal_raw_rows": len(frame) == len(pd.read_csv(args.input, usecols=["matrix_row_id"])),
    }
    (args.output / "13_validation_checks.json").write_text(
        json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_report(args.output, metrics, checks)

    manifest = {
        "created_at_utc": pd.Timestamp.utcnow().isoformat(),
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "rdkit": rdBase.rdkitVersion,
        "input_sha256": {
            str(args.input): sha256(args.input),
            str(PUBCHEM): sha256(PUBCHEM),
            str(TAXON_MAP): sha256(TAXON_MAP),
            str(MODEL_QC): sha256(MODEL_QC),
        },
        "metrics": metrics,
        "checks": checks,
    }
    (args.output / "14_finalization_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
