# FMH BRNP T2D data and code

This is the private engineering repository supporting the manuscript *Biotransformation-resolved network pharmacology for food-medicine homology materials in type 2 diabetes*. It contains portable analysis code, configurations, author-generated derived results and Figures 2–6. A separate clean public candidate is built from `PUBLIC_RELEASE_ALLOWLIST.tsv`; manuscripts, author forms, reviewer-specific assessments, signatures and restricted third-party inputs are not copied into that candidate.

The RCT synthesis and Figures 2–6 are reproducible from redistributed structured inputs. Metabolism-aware v2, ranking sensitivity, chemical-identity audit, AGORA2, Wei 2025 and three-cohort workflows have explicitly documented external-input boundaries. See `RUN_REPRODUCTION.md`.

## Integrity and privacy

Do not make this private engineering repository public. Use the clean candidate described in `FINAL_PUBLIC_RELEASE_AUDIT.md`. Signed reviewer forms, trial PDFs, raw GBD/IDF extracts, commercial database exports, FASTQ, SILVA and AGORA2 model copies are excluded from the clean package.

Review `DATA_SOURCES.md`, `INPUT_SCHEMAS.md`, `PUBLIC_RELEASE_ALLOWLIST.tsv` and `RUN_REPRODUCTION.md` before rerunning a conditional upstream workflow.

For a clean-candidate integrity check, run `sha256sum -c CLEAN_SHA256SUMS.txt`. `CLEAN_CANDIDATE_MANIFEST.tsv` records each released file's size and SHA256.

## Licence

Author-written code is MIT licensed (LICENSE). Author-generated derived data that the authors are entitled to distribute are CC BY 4.0 (LICENSE-DATA). Third-party materials remain under their original terms.

## Release status

Target tag: v1.0.0-submission. The private GitHub repository exists, but no public repository, GitHub release or Zenodo deposit is created by this build. The final audit deliberately stops before those actions.

Zenodo automatic GitHub archiving requires a public repository. If the GitHub repository remains private, use a manual Zenodo deposit with Restricted or Embargoed file access instead.
