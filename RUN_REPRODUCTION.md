# Reproduction guide

## Scope

This package reproduces the RCT meta-analysis from its structured endpoint table and Figures 2–6 from frozen study-generated source data. The metabolism-aware v2, ranking-sensitivity, chemical-identity and AGORA2 workflows are also supplied, but upstream runs that require third-party integrated tables or AGORA2 models are conditional on user-provided inputs. The three-cohort and Wei workflows require original public/third-party data that are not mirrored here.

Scientific boundaries are unchanged: the outputs prioritise hypotheses and calibrate evidence; they do not validate enzyme activity, Rb1-to-Compound K flux, treatment efficacy or patient-specific recommendations.

## Environment

Tested Python: 3.13.12. Install the recorded Python environment with either:

```bash
conda env create -f environment.yml
conda activate fmh-brnp-t2d
```

or, in an isolated Python 3.13 environment:

```bash
python -m pip install -r requirements.txt
```

For the exact Linux Wei DADA2 environment, use:

```bash
conda env create -f content/wei_2025_calibration/config/environment_locked.yml
conda activate wei2025-dada2
```

The portable `environment.yml` includes R 4.5, DADA2 1.38.0, BiocParallel, vegan, ggplot2, jsonlite and metafor. The exact locked Wei file is the authoritative 16S environment record.

## RCT meta-analysis — fully reproducible

```bash
python content/paper_upgrade_v1/scripts/run_rct_rebuild.py \
  --input content/tier1_screening_output/b47_deficiency/meta_outputs/meta_input_clean.csv \
  --output reproduced/rct
sha256sum reproduced/rct/meta_results_revalidated.tsv
```

Expected core SHA256:

```text
24a84bded9beaa811186b9657dfdfc683aad0cc8573395698359dc8b39f3d463
```

The optional `--study-master` and `--old-results` inputs only regenerate ancillary trial characteristics and implementation comparisons; they are not required for the pooled estimates.

## Metabolism-aware prioritisation v2 — conditional upstream reproduction

Reconstruct the seven inputs listed in `INPUT_SCHEMAS.md`, then run:

```bash
python content/paper_upgrade_v1/metabolism_aware_v2_20260907_013732/scripts/run_metabolism_aware_v2.py \
  --input-dir /path/to/reconstructed_v2_inputs \
  --output reproduced/metabolism_aware_v2
python content/paper_upgrade_v1/metabolism_aware_v2_20260907_013732/tests/test_v2_outputs.py \
  reproduced/metabolism_aware_v2
```

The frozen direct-edge comparator is not overwritten. The random seed is 20260907, the topology null has 1,000 iterations and clinical outcomes are not used to tune weights or alpha.

## Ranking evidence sensitivity — conditional upstream reproduction

```bash
python content/paper_upgrade_v1/kg_evidence_sensitivity_20260908_110558/scripts/run_ranking_evidence_sensitivity.py \
  --v2-dir /path/to/reproduced_or_frozen_v2 \
  --qc-dir /path/to/chemistry_qc \
  --output reproduced/ranking_sensitivity
```

The expected outputs are `01_ranking_scenarios_long.tsv`, `02_ranking_scenario_summary.tsv`, `03_herb_rank_shifts.tsv`, `04_top20_by_scenario.tsv` and `05_regression_diagnostics.tsv`. The private verification regenerated all five byte-identically.

## AGORA2 coverage scan — reproducible with user-supplied models

```bash
python content/paper_upgrade_v1/kg_evidence_sensitivity_20260908_110558/scripts/scan_agora2_full_coverage.py \
  --models /path/to/AGORA2/sbml_files \
  --chem content/paper_upgrade_v1/kg_evidence_sensitivity_20260908_110558/03_agora2/agora2_query_compounds_pubchem.tsv \
  --output reproduced/agora2
```

The frozen run read 7,302 models with zero read failures. Both output TSVs were reproduced byte-identically. Absence of a matched model under these rules is not evidence that a biological strain cannot catalyse the reaction.

## Optional RDKit entity audit

Run `run_rdkit_entity_audit.py --help` and supply all six external/local inputs. The integrated herb–compound input and raw PubChem cache are not redistributed. The private verification reproduced eight frozen core outputs byte-identically. This step is optional upstream provenance, not a prerequisite for rerunning the deposited RCT synthesis or figures.

## Wei 2025 workflow — accession based

Download the 100 paired FASTQ files listed in `content/wei_2025_calibration/metadata/Wei2025_ENA_fastq_manifest.tsv`, verify each recorded MD5, place them under `content/wei_2025_calibration/raw/`, and obtain the DADA2-formatted SILVA NR99 138.2 files. Then run:

```bash
Rscript content/wei_2025_calibration/scripts/06_run_dada2.R content/wei_2025_calibration
Rscript content/wei_2025_calibration/scripts/07_analyze_rb1_groups.R content/wei_2025_calibration
Rscript content/wei_2025_calibration/scripts/12_rarefaction_sensitivity.R content/wei_2025_calibration
```

The repository contains accession metadata and final derived result tables, not FASTQ, filtered reads, SILVA files or intermediate RDS objects.

## Three-cohort analysis

The public candidate contains the adjusted cohort-effect and meta-analysis tables and the recorded session information. The legacy `01_butyrate_core_meta.R` is not the source of the final adjusted abundance-ratio result and is excluded from the clean candidate. Reproduction of cohort-level modelling requires the original curatedMetagenomicData objects, cohort metadata and covariate-processing code; this remains a documented reproducibility limitation rather than being disguised as a runnable workflow.

## Figure regeneration — fully reproducible from deposited derived source data

```bash
FIGROOT=content/paper_upgrade_v1/13_biotransformation_resolved_manuscript_v1/figures_submission_detail_optimized_v1
python "$FIGROOT/scripts/make_figure2_detail_optimized_standalone.py" --output reproduced/Figure2
python "$FIGROOT/scripts/make_figure3_source.py" --output reproduced/Figure3
FMH_FIGURE_OUT=reproduced/Figure4 python "$FIGROOT/scripts/make_figure4_detail_optimized.py"
FMH_FIGURE_OUT=reproduced/Figure5 python "$FIGROOT/scripts/make_figure5_source.py"
FMH_FIGURE_OUT=reproduced/Figure6 python content/paper_upgrade_v1/15_chinese_medicine_submission_v1/figures/scripts/make_figure6_final_signed.py
```

Each builder exports editable PDF/SVG and raster previews/TIFF as implemented. Private QA found identical PDF text for Figures 2–6; 150-dpi rendering was identical for Figures 2, 3, 5 and 6, while Figure 4 differed by at most 1/255 in 0.000859% of pixels because of rendering-rounding only. See `FIGURE_REPRODUCTION_QA.md`.

## Validation

```bash
python -m compileall -q content
sha256sum -c CLEAN_SHA256SUMS.txt
rg -n '/home/|~/.codex|BEGIN (RSA|OPENSSH|EC) PRIVATE KEY|ghp_|github_pat_' .
```

For core result comparisons, prefer byte-level SHA256. If a manifest differs only by a path or timestamp, compare row counts, keys, ordering and numeric columns and document the reason. Stop if a core scientific value changes.

## Restricted and non-redistributed inputs

See `DATA_SOURCES.md`, `INPUT_SCHEMAS.md` and `PUBLIC_RELEASE_ALLOWLIST.tsv`. In particular, the package does not include RCT PDFs, full-text manuscripts, raw GBD/IDF extracts, the 85,110-row integrated herb–compound table, FASTQ, SILVA, AGORA2 model files, commercial database exports, reviewer-specific forms, signed records or secrets.
