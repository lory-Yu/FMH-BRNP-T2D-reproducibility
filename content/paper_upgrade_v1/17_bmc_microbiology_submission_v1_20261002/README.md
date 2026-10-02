# BMC Microbiology release materials

This directory contains the author-generated, redistributable source tables,
figure outputs and portable scripts used for the BMC Microbiology manuscript
“Evidence gating exposes representation gaps in gut microbial biotransformation
inference: four computational case studies”.

The release is intentionally a reconstruction package rather than a mirror of
the authors’ local workspace. It does not contain the manuscript DOCX, cover
letter, raw sequencing reads, SILVA reference files, AGORA2 SBML models,
third-party row-level database exports or internal review records.

## Reproduction boundary

- Figure source tables and scripts are provided under `source_data/` and
  `scripts/`.
- Figure 1 and Figure 2 require locally obtained, license-compliant inputs
  supplied through command-line arguments or local package-relative paths.
- AGORA2 model files and non-redistributed gutMGene-derived inputs must be
  obtained or reconstructed by the user under their source terms.
- Conditional AGORA2 feasibility is a model-compatibility result, not culture
  or in-vivo conversion.
- The Wei tables are a same-source reanalysis of public reads and are not an
  independent validation cohort.

The exact input schema, checksums and interpretation boundaries are documented
in the parent repository. No absolute local filesystem paths are required by
the portable command-line scripts.
