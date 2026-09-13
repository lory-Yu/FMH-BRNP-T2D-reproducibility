# FMH BRNP T2D data and code

This release-candidate package supports the manuscript *Biotransformation-resolved network pharmacology for food-medicine homology materials in type 2 diabetes*. It contains portable analysis code, configurations, author-generated derived results and Figures 2–6. It is assembled from `PUBLIC_RELEASE_ALLOWLIST.tsv`; manuscripts, author forms, reviewer-specific assessments, signatures and restricted third-party inputs are excluded.

The RCT synthesis and Figures 2–6 are **DIRECTLY REPRODUCIBLE** from redistributed structured inputs. Metabolism-aware v2 is **CONDITIONALLY REPRODUCIBLE**: its code, schemas, provenance information, safe derived outputs, checksums and reconstruction instructions are provided, but a complete from-scratch run requires local reconstruction of a non-redistributed third-party database-derived transformation input. Optional upstream reconstruction, AGORA2 model scanning and raw-data reprocessing are classified **REQUIRES EXTERNAL INPUT**. See `RUN_REPRODUCTION.md`.

Certain row-level integrated inputs derived from third-party databases are not redistributed because accessibility does not necessarily imply redistribution permission. Their source databases, input schema, provenance and reconstruction procedure are documented. In particular, neither `herb_transformation_paths_v2.tsv` nor a complete or reduced `2_3_herb_compound_clean.csv` is included.

## Integrity and privacy

Do not make this private engineering repository public. Use the clean candidate described in `FINAL_PUBLIC_RELEASE_AUDIT.md`. Signed reviewer forms, trial PDFs, raw GBD/IDF extracts, commercial database exports, FASTQ, SILVA and AGORA2 model copies are excluded from the clean package.

Review `DATA_SOURCES.md`, `INPUT_SCHEMAS.md`, `PUBLIC_RELEASE_ALLOWLIST.tsv` and `RUN_REPRODUCTION.md` before rerunning a conditional upstream workflow.

For a clean-candidate integrity check, run `sha256sum -c CLEAN_SHA256SUMS.txt`. `CLEAN_CANDIDATE_MANIFEST.tsv` records each released file's size and SHA256.

## Licence

Author-written code is licensed under the MIT licence (`LICENSE`). Explicitly identified, author-generated and redistributable derived data are licensed under CC BY 4.0 within the scope defined in `DATA_LICENSE.md`. All three authors confirmed these licence choices on 13 September 2026. Third-party materials and records with unclear redistribution rights are excluded or remain under their original providers' terms.

## Release status

Target tag: v1.0.0-submission. No public repository, GitHub release or Zenodo deposit has been created. The final audit deliberately stops before those actions.

Zenodo automatic GitHub archiving requires a public repository. If the GitHub repository remains private, use a manual Zenodo deposit with Restricted or Embargoed file access instead.
