# External-input schemas

This document describes inputs that are required by optional upstream workflows but are not redistributed. It is not a licence grant for the source databases.

## Metabolism-aware prioritisation input directory

`run_metabolism_aware_v2.py --input-dir DIR` expects these seven files:

1. `KG_all_FMH_ranked.csv`
2. `2_3_herb_compound_clean.csv`
3. `Supplementary_Table_Xb_gut_microbial_links.csv`
4. `Supplementary_Table_Xb_raw_microbial_329.csv`
5. `Supplementary_Table_Xc_host_enzymatic_links.csv`
6. `KG_edge_evidence.tsv`
7. `KG_entity_and_mapping_audit.tsv`

The integrated herb–compound table is not redistributed because its 85,110 rows combine records from databases with different source terms. The program validates required columns and writes the complete route, ranking and validation outputs. A locally reconstructed input must preserve original identifiers, source-database attribution and the material-name harmonisation used in the study.

## Ranking evidence-sensitivity inputs

`run_ranking_evidence_sensitivity.py` reads:

- `01_routes/herb_transformation_paths_v2.tsv` and `02_rankings/KG_metabolism_aware_v2_two_axis.tsv` below the supplied v2 directory;
- `07_route_entity_crosswalk.tsv` and `11_route_entity_crosswalk_final.tsv` below the supplied chemistry-QC directory.

Row-level route and occurrence tables are excluded by release policy and must be reconstructed locally under the source providers' terms. Frozen safe derived ranking outputs are provided separately.

## RDKit entity audit inputs

The optional chemistry audit requires the six explicit paths shown by `run_rdkit_entity_audit.py --help`. `--herb-chem` must point to the locally reconstructed integrated table. `--pubchem-raw` is a local API cache; it is neither required to be published nor treated as independent evidence. `--fetch-missing` retrieves only missing public PubChem records.

## AGORA2 input

`scan_agora2_full_coverage.py --models DIR` expects a user-obtained directory containing the 7,302 AGORA2 SBML/XML models used in the study. The repository does not redistribute these models. The query table must contain `pubchem_cid`, `inchikey` and `inchi` columns; a four-compound PubChem-derived query table is included.

## Wei 2025 inputs

The project directory passed to the R scripts must contain:

- `metadata/Wei2025_D1_D50_phenotype_SRA_crosswalk.tsv`;
- `metadata/Wei2025_ENA_fastq_manifest.tsv`;
- paired FASTQ files under `raw/`, named as recorded in the manifest;
- the DADA2-formatted SILVA NR99 138.2 taxonomy files expected by `06_run_dada2.R`.

FASTQ and SILVA files are obtained from their original repositories and are not mirrored here.
