# FMH BRNP T2D data and code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23097072.svg)](https://doi.org/10.5281/zenodo.23097072)

This repository supports the manuscript *Evidence gating exposes representation gaps in gut microbial biotransformation inference: four computational case studies*. The BMC Microbiology release `v1.1.0-bmc-microbiology` is archived at Zenodo (https://doi.org/10.5281/zenodo.23097072). It contains portable analysis code, configurations, author-generated derived results and source data for Figures 1–4. It is assembled from an explicit public-release allowlist; manuscripts, author forms, reviewer-specific assessments, signatures and restricted third-party inputs are excluded.

The current BMC release provides direct reproduction paths for the redistributed analyses and Figures 1–4. Legacy RCT synthesis materials and Figures 2–6 from the earlier release remain in the repository. Metabolism-aware v2 is **CONDITIONALLY REPRODUCIBLE**: its code, schemas, provenance information, safe derived outputs, checksums and reconstruction instructions are provided, but a complete from-scratch run requires local reconstruction of a non-redistributed third-party database-derived transformation input. Optional upstream reconstruction, AGORA2 model scanning and raw-data reprocessing are classified **REQUIRES EXTERNAL INPUT**. See `RUN_REPRODUCTION.md`.

Certain row-level integrated inputs derived from third-party databases are not redistributed because accessibility does not necessarily imply redistribution permission. Their source databases, input schema, provenance and reconstruction procedure are documented. In particular, neither `herb_transformation_paths_v2.tsv` nor a complete or reduced `2_3_herb_compound_clean.csv` is included.

## Integrity and privacy

This is the public clean repository described in `FINAL_PUBLIC_RELEASE_AUDIT.md`. Signed reviewer forms, trial PDFs, raw GBD/IDF extracts, commercial database exports, FASTQ, SILVA and AGORA2 model copies are excluded from the clean package.

Review `DATA_SOURCES.md`, `INPUT_SCHEMAS.md`, `PUBLIC_RELEASE_ALLOWLIST.tsv` and `RUN_REPRODUCTION.md` before rerunning a conditional upstream workflow.

For a clean-candidate integrity check, run `sha256sum -c CLEAN_SHA256SUMS.txt`. `CLEAN_CANDIDATE_MANIFEST.tsv` records each released file's size and SHA256.

## Licence

Author-written code is licensed under the MIT licence (`LICENSE`). Explicitly identified, author-generated and redistributable derived data are licensed under CC BY 4.0 within the scope defined in `DATA_LICENSE.md`. All three authors confirmed these licence choices on 13 September 2026. Third-party materials and records with unclear redistribution rights are excluded or remain under their original providers' terms.

## Release status

Current release: `v1.1.0-bmc-microbiology`.

- GitHub: https://github.com/lory-Yu/FMH-BRNP-T2D-reproducibility/releases/tag/v1.1.0-bmc-microbiology
- Zenodo version DOI: https://doi.org/10.5281/zenodo.23097072
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.22731840

The earlier `v1.0.0-submission` archive remains available at https://doi.org/10.5281/zenodo.22731841.
