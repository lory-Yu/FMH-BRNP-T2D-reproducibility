# Final public-release audit

Audit date: 13 September 2026

Internal archive: retained separately as PRIVATE; it is not a public-release source

Clean candidate: `FMH-BRNP-T2D-public-candidate` (local only; no remote configured)

## Decision

**READY_FOR_PUBLIC_RELEASE = NO**

All technical portability and clean-candidate checks described below passed. Public release is held for one governance action: all authors must confirm the MIT licence for author-written code and CC BY 4.0 for explicitly identified distributable author-generated data, together with the final CRediT and AI disclosure.

`AUTHOR_APPROVAL_PENDING` is the only release blocker. `THIRD_PARTY_ROW_LEVEL_DATA_EXCLUDED_BY_POLICY` is a settled data-governance choice, not a technical or scientific failure. The adopted policy is not to redistribute `herb_transformation_paths_v2.tsv` or complete/minimal versions of `2_3_herb_compound_clean.csv`.

- Scientific reproduction failure: **NONE**
- Technical blockers: **NONE**
- Pending governance action: **co-author approval of licensing, CRediT and AI disclosure**

The missing GitHub release and Zenodo DOI are intentional next-stage actions, not failed QA. No repository visibility, release or Zenodo action was performed.

## Completed fixes

| Workflow | Status | Verification |
|---|---|---|
| Metabolism-aware v2 | portable CLI complete | full and minimal local-input reruns reproduced frozen core outputs; invariant test passed |
| RCT synthesis | DIRECTLY REPRODUCIBLE from redistributed structured input | core result byte-identical; SHA256 `24a84bded9beaa811186b9657dfdfc683aad0cc8573395698359dc8b39f3d463` |
| Ranking evidence sensitivity | portable CLI complete | five frozen outputs byte-identical |
| AGORA2 scan | portable CLI with user-supplied models | 7,302 models, zero read failures; two frozen outputs byte-identical |
| RDKit entity audit | portable optional upstream CLI | eight frozen outputs byte-identical with local non-redistributed inputs |
| Figures 2–6 | DIRECTLY REPRODUCIBLE from deposited source data | PDF text identical for all five; Figure 4 only showed negligible raster rounding documented in `FIGURE_REPRODUCTION_QA.md` |

No scientific scoring rule, model, random seed, route definition, pooled effect, confidence interval, rank or figure label was changed by these portability patches.

## Reproducibility status by workflow

- **DIRECTLY REPRODUCIBLE:** RCT synthesis; Figures 2–6.
- **CONDITIONALLY REPRODUCIBLE:** metabolism-aware v2 and ranking sensitivity, because the required integrated row-level transformation input is not redistributed.
- **REQUIRES EXTERNAL INPUT:** AGORA2 scan; Wei 2025 DADA2 and downstream analyses; optional RDKit entity audit; end-to-end three-cohort modelling.
- **Deposited results but incomplete end-to-end public rerun:** final adjusted three-cohort modelling. The clean package contains adjusted cohort effects, pooled results, leave-one-out results and session information. The legacy R script was excluded because it is not the final adjusted abundance-ratio workflow and points to absent RDS objects.

## Restricted-data status

The clean candidate contains none of the following:

- RCT full-text PDFs or signed/reviewer-specific RoB 2 and GRADE files;
- manuscript, author submission forms, internal extraction logs or QC decisions;
- GBD/IHME or IDF raw extracts;
- commercial/restricted raw database exports or the 85,110-row integrated herb–compound table;
- FASTQ/FQ, filtered read mirrors or SILVA database files;
- AGORA2 XML/SBML/MAT/ZIP models;
- API keys, tokens, passwords, private keys or `.env` files.

`herb_transformation_paths_v2.tsv` is marked `EXCLUDE` in the allowlist and was not copied, because it joins route records to multi-database herb occurrences whose accessibility does not establish redistribution permission. Its exclusion is the adopted release policy.

## Absolute-path audit

No personal absolute path or `~/.codex` dependency occurs in the clean candidate's scientific content. Twelve residual hits in the private engineering repository are classified in `LOCAL_PATH_AUDIT.tsv` as:

- two superseded figure wrappers — `INTERNAL_ONLY`, excluded;
- one internal RCT ledger builder — `INTERNAL_ONLY`, excluded;
- seven paths in a frozen historical v2 run manifest — `HISTORICAL_RECORD`, excluded without altering history;
- one local packaging audit literal — `AUDIT_LITERAL`, excluded;
- one legacy three-cohort path — `INTERNAL_ONLY`, excluded.

## Secret, path and large-file checks

- secret-pattern scan: no credential or private-key match in the clean candidate; documented audit patterns and the local variable name `token` were the only broad-pattern false positives;
- personal email scan: no match;
- personal path scan under `content/`: no match;
- prohibited raw-file scan: no FASTQ/FQ, RDS, XML, SBML, MAT, ZIP or DOCX;
- PDFs: exactly the five author-generated final Figures 2–6;
- largest release file: `01_full_model_coverage_scan.tsv`, 936,191 bytes; no large raw-data object was copied.

## Environment status

- tested Python: 3.13.12;
- recorded Python packages: numpy 2.4.3, pandas 2.3.3, scipy 1.17.1, matplotlib 3.10.8, requests 2.32.5, RDKit 2025.9.6, Pillow 12.1.1, python-docx 1.2.0 and PyMuPDF 1.24.14;
- portable environment: `environment.yml`;
- Python package list: `requirements.txt`;
- exact Wei Linux environment: `content/wei_2025_calibration/config/environment_locked.yml`, with DADA2 1.38.0 and no personal prefix.

All public Python scripts passed `python -m compileall -q content` in the clean candidate.

## Clean-candidate validation

The candidate was assembled only from rows marked `KEEP` in `PUBLIC_RELEASE_ALLOWLIST.tsv`; `EXCLUDE` and `OPTIONAL` rows were not copied. It was initialised as a new Git repository without old history or a remote.

Commands executed from the clean candidate included:

```bash
python -m compileall -q content
python content/paper_upgrade_v1/scripts/run_rct_rebuild.py \
  --input content/tier1_screening_output/b47_deficiency/meta_outputs/meta_input_clean.csv \
  --output /tmp/fmh-public-qa-*/rct
python content/paper_upgrade_v1/13_biotransformation_resolved_manuscript_v1/figures_submission_detail_optimized_v1/scripts/make_figure2_detail_optimized_standalone.py --output /tmp/fmh-public-qa-*/Figure2
python content/paper_upgrade_v1/13_biotransformation_resolved_manuscript_v1/figures_submission_detail_optimized_v1/scripts/make_figure3_source.py --output /tmp/fmh-public-qa-*/Figure3
FMH_FIGURE_OUT=/tmp/fmh-public-qa-*/Figure4 python content/paper_upgrade_v1/13_biotransformation_resolved_manuscript_v1/figures_submission_detail_optimized_v1/scripts/make_figure4_detail_optimized.py
FMH_FIGURE_OUT=/tmp/fmh-public-qa-*/Figure5 python content/paper_upgrade_v1/13_biotransformation_resolved_manuscript_v1/figures_submission_detail_optimized_v1/scripts/make_figure5_source.py
FMH_FIGURE_OUT=/tmp/fmh-public-qa-*/Figure6 python content/paper_upgrade_v1/15_chinese_medicine_submission_v1/figures/scripts/make_figure6_final_signed.py
```

The RCT output matched the frozen SHA256 exactly. `pdftotext` comparison reported `PDF_TEXT_IDENTICAL` for Figures 2, 3, 4, 5 and 6.

### Post-policy final QA rerun

After adopting the do-not-redistribute policy for the third-party-derived row-level tables, final QA was rerun on 13 September 2026:

- allowlist match: 120 `KEEP` paths and 120 candidate files; no missing or unlisted file;
- integrity manifest: 118 file records plus header; checksum list: 119 entries including the manifest but excluding the checksum file itself;
- `sha256sum -c CLEAN_SHA256SUMS.txt`: all entries passed;
- Python compile audit: passed with the byte-code cache redirected outside the candidate;
- Git remote: none configured;
- prohibited exact-file scan: no `herb_transformation_paths_v2.tsv`, `2_3_herb_compound_clean.csv`, internal extraction/QC log or submission DOCX;
- restricted-format scan: no FASTQ/FQ, SILVA database file, AGORA2 XML/SBML/MAT/ZIP model or `.env` file;
- PDF inventory: only the five author-generated final Figures 2–6;
- high-confidence credential/private-key scan: no match;
- personal email scan: no match;
- large-file audit: no file exceeded 5 MiB; the largest remained the 936,191-byte derived AGORA2 scan table.

The filename `make_figure6_final_signed.py` contains the word `signed` but is a portable figure builder; inspection found no signed form or reviewer-specific record in the candidate.

## Relevant private-repository commits

- `dd08927` — Add missing metabolism-aware v2 transformation paths
- `b997895` — Make metabolism-aware v2 workflow portable
- `6bd1b42` — Remove local paths from Wei2025 environment
- `49de02e` — Make RCT meta-analysis workflow portable
- `cc2ed08` — Make ranking sensitivity workflow portable
- `00d90be` — Make AGORA2 scan workflow portable
- `bdb4da6` — Make chemical entity audit portable
- `6b9b2e3` — Complete portable Figure 2–6 workflows
- `f6d5987` — Document reproducible workflows and release boundaries

## Files excluded and why

The authoritative file-by-file decision is `PUBLIC_RELEASE_ALLOWLIST.tsv`. Broad categories are manuscript/author internal records, reviewer-specific or signed assessments, internal extraction/QC provenance, restricted raw data, source-term-pending joins, historical personal-path manifests, obsolete wrappers and local packaging tools.

## Required author action before changing the decision to YES

1. All three authors must confirm the MIT code licence, scoped CC BY 4.0 data licence, final CRediT statement and final AI-assisted-use disclosure.
2. After that confirmation, rerun clean-candidate integrity, allowlist, secret, restricted-data and checksum checks and change the decision to `READY_FOR_PUBLIC_RELEASE = YES` only if all still pass.
3. Stop after reporting `PUBLIC RELEASE CANDIDATE PASSED FINAL QA`; do not create a repository, Release or Zenodo deposit without a separate explicit authorisation.
