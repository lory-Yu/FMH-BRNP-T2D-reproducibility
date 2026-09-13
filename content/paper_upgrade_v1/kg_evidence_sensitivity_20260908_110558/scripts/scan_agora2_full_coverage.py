#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT = SCRIPT_DIR.parent


def inventory_hash(files: list[Path]) -> str:
    h = hashlib.sha256()
    for path in files:
        stat = path.stat()
        h.update(f"{path.name}\t{stat.st_size}\t{stat.st_mtime_ns}\n".encode())
    return h.hexdigest()


def patterns(chem_path: Path) -> dict[str, list[bytes]]:
    chem = pd.read_csv(chem_path, sep="\t", dtype=str).fillna("")
    cids = [str(x).encode() for x in chem["pubchem_cid"]]
    keys = [str(x).lower().encode() for x in chem["inchikey"] if x]
    inchis = [str(x).lower().encode() for x in chem["inchi"] if x]
    return {
        "route_name": [x.encode() for x in [
            "ginsenoside", "ginsenosidase", "protopanaxadiol", "compound k",
            "compound_k", "compound-k", "ginsenoside rb1", "ginsenoside_rd",
            "ginsenoside f2", "ginsenoside compound k",
        ]],
        "pubchem_cid": cids,
        "inchikey": keys,
        "inchi": inchis,
        "ec_191": [b"3.2.1.191", b"3_2_1_191", b"3-2-1-191"],
        "ec_192": [b"3.2.1.192", b"3_2_1_192", b"3-2-1-192"],
        "ec_193": [b"3.2.1.193", b"3_2_1_193", b"3-2-1-193"],
        "ec_194": [b"3.2.1.194", b"3_2_1_194", b"3-2-1-194"],
        "ec_195": [b"3.2.1.195", b"3_2_1_195", b"3-2-1-195"],
        "ec_generic_21": [b"3.2.1.21", b"3_2_1_21", b"3-2-1-21"],
    }


PATTERNS: dict[str, list[bytes]] = {}
MAX_PATTERN = 0


def scan_file(path: Path) -> dict:
    flags = {key: False for key in PATTERNS}
    row = {"model_file": path.name, "bytes": path.stat().st_size, "read_ok": False, "read_error": ""}
    carry = b""
    try:
        with path.open("rb") as handle:
            while True:
                block = handle.read(4 * 1024 * 1024)
                if not block:
                    break
                data = (carry + block).lower()
                for key, terms in PATTERNS.items():
                    if not flags[key] and any(term in data for term in terms):
                        flags[key] = True
                carry = data[-(MAX_PATTERN - 1) :] if MAX_PATTERN > 1 else b""
        row["read_ok"] = True
    except Exception as exc:
        row["read_error"] = f"{type(exc).__name__}: {exc}"
    row.update(flags)
    row["any_route_specific_ec"] = any(flags[f"ec_{x}"] for x in range(191, 196))
    row["any_target_entity_identifier_or_name"] = any(flags[x] for x in ["route_name", "pubchem_cid", "inchikey", "inchi"])
    row["target_positive"] = row["any_route_specific_ec"] or row["any_target_entity_identifier_or_name"]
    return row


def targeted_parse(path: Path) -> dict:
    result = {"model_file": path.name, "xml_parse_ok": False, "xml_parse_error": "", "target_species_count": 0, "target_reaction_count": 0, "target_reactions_with_gpr": 0, "target_species_ids": "", "target_reaction_ids": ""}
    try:
        from lxml import etree
        tree = etree.parse(str(path))
        species_ids, reaction_ids, reactions_gpr = [], [], 0
        needles = [p for key in ["route_name", "pubchem_cid", "inchikey", "inchi"] for p in PATTERNS[key]]
        for elem in tree.iter():
            local = etree.QName(elem).localname
            if local not in {"species", "reaction"}:
                continue
            text = etree.tostring(elem, encoding="utf-8").lower()
            if not any(p in text for p in needles):
                continue
            identifier = elem.get("id", "") or elem.get("name", "")
            if local == "species":
                species_ids.append(identifier)
            else:
                reaction_ids.append(identifier)
                if any(etree.QName(x).localname in {"geneProductAssociation", "geneProductRef"} for x in elem.iter()):
                    reactions_gpr += 1
        result.update({"xml_parse_ok": True, "target_species_count": len(species_ids), "target_reaction_count": len(reaction_ids), "target_reactions_with_gpr": reactions_gpr, "target_species_ids": ";".join(species_ids[:100]), "target_reaction_ids": ";".join(reaction_ids[:100])})
    except Exception as exc:
        result["xml_parse_error"] = f"{type(exc).__name__}: {exc}"
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan a user-supplied AGORA2 SBML directory for target entities and EC tokens.")
    parser.add_argument("--models", type=Path, required=True, help="Directory containing the 7,302 AGORA2 XML models.")
    parser.add_argument("--chem", type=Path, required=True, help="TSV containing pubchem_cid, inchikey and inchi query fields.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output root; 03_agora2 and logs are created below it.")
    return parser.parse_args()


def main() -> None:
    global PATTERNS, MAX_PATTERN
    args = parse_args()
    models = args.models.resolve()
    chem = args.chem.resolve()
    out = args.output.resolve()
    agora_out = out / "03_agora2"
    logs_out = out / "logs"
    agora_out.mkdir(parents=True, exist_ok=True)
    logs_out.mkdir(parents=True, exist_ok=True)
    PATTERNS = patterns(chem)
    MAX_PATTERN = max(len(p) for group in PATTERNS.values() for p in group)

    files = sorted(models.glob("*.xml"))
    if len(files) != 7302:
        raise RuntimeError(f"Expected 7302 XML models, found {len(files)}")
    workers = min(8, max(1, os.cpu_count() or 1))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(scan_file, files))
    scan = pd.DataFrame(rows)
    scan.to_csv(agora_out / "01_full_model_coverage_scan.tsv", sep="\t", index=False)

    positive_files = [models / name for name in scan.loc[scan["target_positive"], "model_file"]]
    parsed_rows = [targeted_parse(path) for path in positive_files]
    parsed = pd.DataFrame(parsed_rows, columns=["model_file", "xml_parse_ok", "xml_parse_error", "target_species_count", "target_reaction_count", "target_reactions_with_gpr", "target_species_ids", "target_reaction_ids"])
    parsed.to_csv(agora_out / "02_target_positive_xml_parse.tsv", sep="\t", index=False)

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "model_directory": "<user-supplied>",
        "model_directory_name": models.name,
        "compound_query_sha256": hashlib.sha256(chem.read_bytes()).hexdigest(),
        "scan_method": "streaming case-insensitive raw SBML token scan; targeted XML parse only for positive files",
        "xml_files_expected": 7302,
        "xml_files_scanned": int(len(scan)),
        "files_read_ok": int(scan["read_ok"].sum()),
        "files_read_failed": int((~scan["read_ok"]).sum()),
        "total_bytes_scanned": int(scan["bytes"].sum()),
        "inventory_sha256_names_sizes_mtimes": inventory_hash(files),
        "models_with_route_name_terms": int(scan["route_name"].sum()),
        "models_with_any_target_pubchem_cid_token": int(scan["pubchem_cid"].sum()),
        "models_with_any_target_inchikey": int(scan["inchikey"].sum()),
        "models_with_any_target_inchi": int(scan["inchi"].sum()),
        "models_with_EC_3_2_1_191": int(scan["ec_191"].sum()),
        "models_with_EC_3_2_1_192": int(scan["ec_192"].sum()),
        "models_with_EC_3_2_1_193": int(scan["ec_193"].sum()),
        "models_with_EC_3_2_1_194": int(scan["ec_194"].sum()),
        "models_with_EC_3_2_1_195": int(scan["ec_195"].sum()),
        "models_with_generic_EC_3_2_1_21_token": int(scan["ec_generic_21"].sum()),
        "target_positive_models": int(scan["target_positive"].sum()),
        "target_positive_xml_parse_attempted": int(len(parsed)),
        "target_positive_xml_parse_failed": int((~parsed["xml_parse_ok"]).sum()) if len(parsed) else 0,
        "target_species_found_after_xml_parse": int(parsed["target_species_count"].sum()) if len(parsed) else 0,
        "target_reactions_found_after_xml_parse": int(parsed["target_reaction_count"].sum()) if len(parsed) else 0,
        "target_reactions_with_gpr": int(parsed["target_reactions_with_gpr"].sum()) if len(parsed) else 0,
        "interpretation_boundary": "A rule-bounded non-match in this local AGORA2 release is a model-coverage result, not evidence that the represented strain lacks real ginsenoside metabolism.",
    }
    (logs_out / "agora2_full_scan_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
