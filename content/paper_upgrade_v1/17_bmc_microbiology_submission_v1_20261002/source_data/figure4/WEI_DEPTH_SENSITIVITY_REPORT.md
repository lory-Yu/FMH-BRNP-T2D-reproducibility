# Wei 2025 library-depth sensitivity audit

This low-cost audit uses the frozen clean DADA2 ASV table; FASTQ files were not reprocessed.

- High and Low groups had median retained library sizes of 37,699 and 36,960 reads, respectively (two-sided Mann–Whitney P=0.85).
- In small-sample t-based HC3 models adjusted for age, sex and log library size, the High–Low difference was 78.83 ASVs (95% CI 6.03 to 151.63; P=0.0357; two-metric BH q=0.0387) and 78.41 for Chao1 (95% CI 4.63 to 152.20; P=0.0387; q=0.0387).
- Repeated rarefaction used 500 independently sampled tables at 30,000 reads per donor and seed 20260903. All replicates retained a positive adjusted Observed-ASV difference; 100.0% had P<0.05 (median P=0.0353; maximum P=0.0445).
- Leave-one-donor-out coefficients remained positive (Observed ASV 66.02 to 94.97; Chao1 65.16 to 94.58), although the maximum leave-one-out P values were 0.0735 and 0.0785. This supports directional stability but not immunity to small-sample influence.

The result should be interpreted only as a richness difference between source-defined extreme converter groups. It is not an independent replication of Rb1 conversion, an enzyme measurement, or evidence about T2D or metformin.
