# Evidence-constrained FMH ranking audit: frozen specification

Date: 2026-09-08

## Inputs

- Frozen 104-material direct ranking: `model_output/kg_tier1_top20/KG_all_FMH_ranked.csv`.
- Metabolism-aware v2 paths/ranks: `paper_upgrade_v1/metabolism_aware_v2_20260907_013732`.
- Official herb–compound table: `99_整合输入库/03_herb_compound_layer/2_3_herb_compound_clean.csv`.
- AGORA2 SBML directory: the locally deposited 7,302-model directory recorded in the earlier run manifest.
- BbBgl/BglX evidence audit: `paper_upgrade_v1/evidence_completion_human_gut_20260908_005251`.

All inputs are read only. All new files are written under this directory.

## Chemical identity levels

1. **Exact structure**: sanitized RDKit isomeric canonical SMILES/InChIKey.
2. **Connectivity identity**: canonical structure with stereochemistry removed.
3. **Parent identity**: largest fragment, normalized and uncharged before canonicalization.
4. **Tautomer-parent identity**: canonical tautomer of the standardized parent.

Records are not silently merged across these levels. Multi-fragment structures, charges,
unassigned stereocentres, within-CID structure conflicts, cross-CID exact duplicates and
tautomer/connectivity collapses are reported separately.

## Ranking sensitivity scenarios

- `original_primary`: deposited v2 primary path weights, unchanged.
- `entity_resolved`: retain a path only when parent and product CIDs have a parseable,
  internally connectivity-consistent structure.
- `entity_resolved_structure_deduplicated`: additionally collapse chemically identical
  parent–product–target routes using standardized parent structure keys, retaining the
  maximum evidence weight rather than summing database duplicates.
- `bbBgl_dispute_removed`: remove any Rb1/Rd/F2/Compound-K route attributed to the
  disputed BbBgl molecular identity. If no such route exists in v2, report zero impact.

The existing residualization formula is reused: log1p route score is regressed on
log1p herb compound count and log1p route-substrate opportunities. RCT outcomes are not
used to select filters, weights or scenarios.

## AGORA2 scan

All 7,302 local SBML files are scanned in a streaming, case-insensitive manner for:

- route compound names and conservative synonyms;
- PubChem CIDs and InChIKey/InChI identifiers for Rb1, Rd, F2 and Compound K;
- EC 3.2.1.191–195 and the generic EC 3.2.1.21;
- reaction/metabolite annotations containing ginsenoside/protopanaxadiol terms.

Only target-positive files are subjected to targeted XML parsing/GPR extraction. A raw
string non-match is reported as a version- and rule-bounded coverage result, not as proof
that a strain lacks the biological capability.

## Prohibited interpretations

- Structural identity does not demonstrate herb content, bioavailability or enzyme activity.
- Model coverage does not demonstrate or refute strain activity.
- Ranking change does not demonstrate clinical efficacy.
- BbBgl sequence conflict does not invalidate the reported biochemical observation; it
  prevents unambiguous attribution of that observation to a public protein accession.
