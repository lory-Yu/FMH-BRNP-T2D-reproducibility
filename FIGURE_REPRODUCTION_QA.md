# Figure 2–6 reproduction QA

## Scope and evidence contract

The public figure workflow uses Python/matplotlib only. Each figure is generated from frozen, study-derived source tables; no values are parsed from the submitted PDF files.

| Figure | Results-level question | Evidence role | Interpretation boundary |
|---|---|---|---|
| 2 | Does evidence-qualified transformation information alter direct-edge prioritisation, and how robust is that alteration? | Framework implementation, conditional AGORA2 coverage, two-axis comparison and stress tests | Hypothesis prioritisation, not clinical efficacy or measured ginsenoside conversion |
| 3 | How much of the chemical and enzyme route evidence survives identity, provenance and model-coverage audit? | Entity adjudication, ranking sensitivity, sequence anchoring and AGORA2 coverage | Evidence audit, not enzyme-activity validation |
| 4 | Are butyrate-pathway associations consistent across T2D cohorts and metformin-related sensitivity analyses? | Cohort composition, cohort effects, random-effects summaries and bounded sensitivity analysis | Functional potential, not butyrate concentration, flux or causal metformin interaction |
| 5 | Which Wei 2025 converter-group microbial patterns are recovered after public-data reprocessing? | QC, richness, beta diversity, candidate genera and global multiplicity | Same-source reanalysis, not independent external validation or enzyme assay |
| 6 | How mature and robust is the clinical evidence after outcome-specific synthesis and TRIAL028 sensitivity analysis? | Screening, evidence attribution, outcome-specific forests, sensitivity and RoB 2 | Clinical calibration, not validation of the KG ranking |

## Reproduction comparison

The five plotting workflows were executed in a clean test-output directory. The submitted PDFs and rebuilt PDFs had identical page dimensions and identical text extracted with layout preservation.

| Figure | PDF text | 150-dpi rendered comparison | Interpretation |
|---|---|---|---|
| 2 | identical | pixel-identical | exact visual reproduction |
| 3 | identical | pixel-identical | exact visual reproduction |
| 4 | identical | maximum difference 1/255 in 0.000859% of pixels | renderer rounding only; scientific geometry unchanged |
| 5 | identical | pixel-identical | exact visual reproduction |
| 6 | identical | pixel-identical | exact visual reproduction |

PDF byte hashes differ because regenerated files contain new document metadata. This does not affect text, geometry or rendered content.

## Automated QA

- Python source preflight: Figures 2, 3, 4 and 6 passed 21/21 checks. Figure 5 passed all scientific and export checks with one non-blocking warning because its physical width is established by the rendered PDF rather than a separately declared scalar.
- Rendered font floor: all five PDFs are auditable and contain no text below 5 pt. Minimum detected sizes were 5.4, 5.1, 5.5, 5.1 and 5.1 pt for Figures 2–6, respectively.
- Collision audit: all five PDFs passed with zero failures and zero warnings.
- Panel alignment: every plotting script generated a passing alignment JSON before export.

## Source-data exclusions

No RCT full text, signed RoB 2/GRADE records, FASTQ data, SILVA database, AGORA2 XML/SBML models or restricted herb–compound database exports are required or included for figure regeneration. Figure 3 reads the frozen QA totals rather than redistributing the full chemical-entity audit table because only the aggregate strict-usable count is displayed.
