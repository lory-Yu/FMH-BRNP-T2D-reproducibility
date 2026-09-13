# Data sources and redistribution status

The labels below distinguish access from redistribution. `PUBLIC` means that a source can be accessed publicly; it does not automatically mean that every bulk export can be republished. `DERIVED` denotes study-generated summaries or analysis outputs. `ACCESSION-ONLY` means that this repository publishes identifiers and acquisition metadata, not the source files. `NOT REDISTRIBUTED` and `RESTRICTED/LOCAL` require users to obtain the source independently under its own terms.

Under the adopted release policy, row-level `herb_transformation_paths_v2.tsv` and complete or reduced versions of `2_3_herb_compound_clean.csv` are **NOT REDISTRIBUTED**. This is a data-governance decision, not a failed scientific reproduction test. The package instead supplies source names, provenance, schemas, safe derived outputs, checksums and reconstruction instructions.

| Source | Project use | Release class | Material provided here | Version or identifier recorded in the project |
|---|---|---|---|---|
| CMAUP, SymMap, HIT, TCMSP, TCMID, TCM-BANK, HERB, ETCM, BATMAN, TCMIP and ccTCM | herb–compound occurrence integration | RESTRICTED/LOCAL | database names, methods, schemas and aggregate derived outputs only | per-source versions were not consistently recorded in the public candidate |
| PubChem | compound identifiers and identity adjudication | PUBLIC / DERIVED | selected identifiers, audit summaries and a four-compound AGORA2 query table | retrieval provenance is recorded in the chemistry audit; raw API cache is not redistributed |
| LOTUS | independent herb–compound occurrence flag | PUBLIC / NOT REDISTRIBUTED | derived verification flags only | version not independently established in this candidate |
| GutMGene | microbial transformation/target source records | PUBLIC / DERIVED | evidence-qualified route summaries subject to source-term review | source identifiers and PMIDs retained in private provenance |
| BRENDA and HMDB | enzyme/metabolite bridges | NOT REDISTRIBUTED | aggregate or evidence-audit outputs only | source versions not independently established in this candidate |
| Open Targets, UniProt and Reactome | direct target/pathway evidence and bridges | PUBLIC / DERIVED | study-generated scores and summaries | versions are those recorded in the frozen analysis provenance |
| DrugBank | DDI penalty source | RESTRICTED/LOCAL | no raw export | version not redistributed |
| ENA/SRA | Wei 2025 amplicon reads | ACCESSION-ONLY | PRJNA1268742 crosswalk, FTP URLs, MD5 and byte counts | project accession PRJNA1268742; 50 donors and 100 paired-end files |
| SILVA | DADA2 taxonomy assignment | NOT REDISTRIBUTED | database citation and configuration only | SILVA NR99 138.2, DADA2-formatted; DOI 10.5281/zenodo.14169026 |
| AGORA2 / Virtual Metabolic Human | genome-scale metabolic models | NOT REDISTRIBUTED | study-generated model-level scan table; no XML/SBML/MAT/ZIP | 7,302-model collection used in the frozen scan |
| Qin 2012, Karlsson 2013 and MetaCardis | curated metagenomic pathway analyses | DERIVED / SOURCE NOT REDISTRIBUTED | adjusted cohort effects, pooled results and session information | cohort labels and sample accounting are documented with the outputs |
| RCT reports and registrations | HbA1c/FPG evidence synthesis | DERIVED / FULL TEXT NOT REDISTRIBUTED | structured endpoint input and pooled results; no PDFs | PROSPERO CRD420251250149 |
| GBD/IHME | Figure 1 burden estimates | NOT REDISTRIBUTED | no raw extract and no Figure 1 source table in the public candidate | source-specific attribution/permission applies |
| IDF Diabetes Atlas | Figure 1 present/projection estimates | NOT REDISTRIBUTED | no raw extract and no Figure 1 source table in the public candidate | source-specific attribution/permission applies |
| Natural Earth | Figure 1 geographic polygons | PUBLIC | not bundled in this candidate | 1:50m administrative boundaries used in the private figure workflow |

## Shared-source rule

Multiple database records that trace to the same original publication are not treated as independent validation. Database-derived associations, homology annotations, molecular modelling and experimental evidence remain separate evidence classes.

## Prohibited substitutions

The accession manifest does not substitute for donor-level Rb1 or Compound K concentrations. Taxonomic profiles do not constitute enzyme-activity or metabolic-flux measurements. AGORA2 coverage does not establish that a strain performs a reaction in vivo.
