# Metabolism-aware v2: frozen analysis specification

## Purpose

This analysis connects the existing microbial and putative host transformation
tables to the 104-material FMH universe without altering the frozen direct-edge
ranking. It produces an evidence axis and prespecified rank sensitivities. It
does not estimate clinical efficacy, individual response or biochemical flux.

## Frozen comparator

- `model_output/kg_tier1_top20/KG_all_FMH_ranked.csv` is the immutable
  direct-edge baseline (v1).
- A degree-adjusted comparator is derived by regressing its frozen composite
  score on `log1p(n_compounds)` and `log1p(n_t2d_targets)`. This is a
  sensitivity comparator, not a replacement for v1.

## Route construction

Routes are joined by normalized PubChem CID:

`herb --parent CID--> substrate --transformation--> product --association--> target`

The deposited `herbs` strings in terminal transformation tables are not used to
construct herb-route paths. Herb-compound occurrences are reconstructed from
the official herb-compound table.

## Evidence classes

- `human_database_subset`: microbial upstream model is human only; the target
  model is human; transformation and target PMIDs are present; the deposited
  target association is causal; the route is not a self-loop; and the
  herb-compound occurrence has `lotus_verified == True`.
- `mixed_human_mouse`: microbial upstream evidence contains both human and
  mouse records.
- `mouse_only`: microbial upstream evidence contains mouse records only.
- `host_cascade_unverified`: putative host cascade without a mapping conflict;
  primary reaction and terminal target provenance remain incomplete.
- `excluded_mapping_conflict`: audited EC/reaction conflict or transporter
  self-loop; weight is zero.

`human_database_subset` is database-integrated evidence with traceable source
identifiers. It is not relabelled as a fully adjudicated experimental human
mechanism.

## Scores

The primary output is two-dimensional:

- `D_h`: frozen direct-edge score and its degree-adjusted sensitivity.
- `T_h`: route-evidence score derived from unique CID-connected paths.

Route weights are prespecified for sensitivity analysis only:

- human-only microbial model: 1.00
- mixed human/mouse microbial model: 0.70
- mouse-only microbial model: 0.35
- non-conflicting but provenance-incomplete host cascade: 0.20
- mapping conflict: 0
- causal target association multiplier: 1.00
- correlational target association multiplier: 0.50
- LOTUS-verified occurrence multiplier: 1.00
- other database occurrence multiplier: 0.50

Every path is additionally downweighted by
`1 / sqrt(log1p(substrate herb degree) * log1p(product target degree))`.
Rows and repeated publications do not add weight. Source counts are reported
descriptively.

The transformation score is residualized against `log1p(n_compounds)` and
`log1p(n_route_substrates)` to expose remaining topology dependence.

## Ranking outputs

1. Frozen v1 rank (unchanged).
2. Two-axis evidence table (primary scientific output).
3. Evidence-tier-first research-priority ordering.
4. Alpha sensitivity with alpha in {0, 0.25, 0.50, 0.75, 1.00}; alpha is not
   selected using RCT results.

The alpha ranking is exploratory and must not be labelled as validated or
clinically predictive.

## Validation

- Confirm the frozen file is internally ordered by its frozen score.
- Confirm CID join cardinality and retain unmatched/excluded rows.
- Confirm mapping-conflict routes contribute zero.
- Confirm alpha=0 reproduces the degree-adjusted direct comparator.
- Report rank correlations, top-20 overlaps and degree correlations.
- Run a degree-preserving herb-substrate topology null with 1,000 iterations.
- Run prespecified evidence-scenario and alpha sensitivity analyses.
- Do not use RCT outcomes to tune weights, filters or alpha.

## Interpretation boundary

The outputs prioritize hypotheses under the deposited evidence structure. They
do not demonstrate transformation rates, enzyme activity, treatment effects or
patient-specific recommendations.
