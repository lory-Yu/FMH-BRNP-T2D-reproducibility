#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
  library(BiocParallel)
})

# The managed desktop environment blocks the local socket ports used by the
# default BiocParallel backend. A frozen SerialParam avoids hidden socket
# failures and prioritizes reproducibility over speed for this 50-sample run.
register(SerialParam(progressbar = FALSE), default = TRUE)

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath("wei_2025_calibration")
raw_dir <- file.path(project, "raw")
filtered_dir <- file.path(project, "filtered")
results_dir <- file.path(project, "results")
figures_dir <- file.path(project, "figures")
provenance_dir <- file.path(project, "provenance")
logs_dir <- file.path(project, "logs")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
log_connection <- file(file.path(logs_dir, "06_run_dada2.log"), open = "wt")
sink(log_connection, type = "output")
sink(log_connection, type = "message")

crosswalk <- read.delim(
  file.path(project, "metadata", "Wei2025_D1_D50_phenotype_SRA_crosswalk.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- read.delim(
  file.path(project, "metadata", "Wei2025_ENA_fastq_manifest.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- merge(manifest, crosswalk[, c("donor_id", "donor_number")], by = "donor_id", all.x = TRUE)
manifest <- manifest[order(manifest$donor_number, manifest$mate), ]

samples <- paste0("D", seq_len(50))
fn_f <- file.path(raw_dir, manifest$filename[manifest$mate == 1])
fn_r <- file.path(raw_dir, manifest$filename[manifest$mate == 2])
names(fn_f) <- manifest$donor_id[manifest$mate == 1]
names(fn_r) <- manifest$donor_id[manifest$mate == 2]
fn_f <- fn_f[samples]
fn_r <- fn_r[samples]
stopifnot(all(file.exists(fn_f)), all(file.exists(fn_r)))

filt_f <- file.path(filtered_dir, paste0(samples, "_F_filt.fastq.gz"))
filt_r <- file.path(filtered_dir, paste0(samples, "_R_filt.fastq.gz"))
names(filt_f) <- samples
names(filt_r) <- samples

started <- Sys.time()
set.seed(20250903)

pdf(file.path(figures_dir, "Figure_S1_quality_profiles_raw.pdf"), width = 12, height = 8)
print(plotQualityProfile(fn_f[c(1, 10, 20, 30, 40, 50)]))
print(plotQualityProfile(fn_r[c(1, 10, 20, 30, 40, 50)]))
dev.off()

filter_out <- filterAndTrim(
  fn_f,
  filt_f,
  fn_r,
  filt_r,
  trimLeft = c(20, 20),
  truncLen = c(280, 230),
  maxN = 0,
  maxEE = c(2, 5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE,
  verbose = TRUE
)
rownames(filter_out) <- samples
write.table(
  data.frame(donor_id = rownames(filter_out), filter_out, check.names = FALSE),
  file.path(results_dir, "dada2_filter_tracking.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (any(filter_out[, "reads.out"] == 0)) stop("At least one sample has zero reads after filtering")

err_f <- learnErrors(filt_f, multithread = FALSE, nbases = 1e8, randomize = TRUE)
err_r <- learnErrors(filt_r, multithread = FALSE, nbases = 1e8, randomize = TRUE)
saveRDS(err_f, file.path(results_dir, "dada2_error_model_forward.rds"))
saveRDS(err_r, file.path(results_dir, "dada2_error_model_reverse.rds"))

pdf(file.path(figures_dir, "Figure_S2_error_models.pdf"), width = 12, height = 6)
print(plotErrors(err_f, nominalQ = TRUE))
print(plotErrors(err_r, nominalQ = TRUE))
dev.off()

derep_f <- derepFastq(filt_f, verbose = TRUE)
derep_r <- derepFastq(filt_r, verbose = TRUE)
names(derep_f) <- samples
names(derep_r) <- samples
dada_f <- dada(derep_f, err = err_f, multithread = FALSE, pool = "pseudo")
dada_r <- dada(derep_r, err = err_r, multithread = FALSE, pool = "pseudo")
mergers <- mergePairs(
  dada_f, derep_f, dada_r, derep_r,
  minOverlap = 20,
  maxMismatch = 0,
  verbose = TRUE
)
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)

get_n <- function(x) sum(getUniques(x))
tracking <- cbind(
  filter_out,
  denoisedF = sapply(dada_f, get_n),
  denoisedR = sapply(dada_r, get_n),
  merged = sapply(mergers, get_n),
  nonchim = rowSums(seqtab_nochim)
)
tracking <- data.frame(donor_id = rownames(tracking), tracking, check.names = FALSE)
write.table(tracking, file.path(results_dir, "dada2_full_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

asv_ids <- paste0("ASV", seq_len(ncol(seqtab_nochim)))
sequences <- colnames(seqtab_nochim)
colnames(seqtab_nochim) <- asv_ids
asv_counts <- data.frame(donor_id = rownames(seqtab_nochim), seqtab_nochim, check.names = FALSE)
write.table(asv_counts, file.path(results_dir, "asv_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.vector(rbind(paste0(">", asv_ids), sequences)), file.path(results_dir, "asv_sequences.fasta"))

silva_genus <- file.path(project, "reference", "silva_nr99_v138.2_toGenus_trainset.fa.gz")
silva_species <- file.path(project, "reference", "silva_v138.2_assignSpecies.fa.gz")
stopifnot(file.exists(silva_genus), file.exists(silva_species))
taxa <- assignTaxonomy(sequences, silva_genus, multithread = FALSE, tryRC = TRUE)
taxa <- addSpecies(taxa, silva_species, allowMultiple = FALSE, tryRC = TRUE)
rownames(taxa) <- asv_ids

taxa_df <- data.frame(asv_id = rownames(taxa), sequence = sequences, taxa, check.names = FALSE)
write.table(taxa_df, file.path(results_dir, "asv_taxonomy_silva138_2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

keep <- !is.na(taxa[, "Kingdom"]) & taxa[, "Kingdom"] %in% c("Bacteria", "Archaea")
if ("Order" %in% colnames(taxa)) keep <- keep & (is.na(taxa[, "Order"]) | taxa[, "Order"] != "Chloroplast")
if ("Family" %in% colnames(taxa)) keep <- keep & (is.na(taxa[, "Family"]) | taxa[, "Family"] != "Mitochondria")
seqtab_clean <- seqtab_nochim[, keep, drop = FALSE]
taxa_clean <- taxa[keep, , drop = FALSE]

genus_labels <- taxa_clean[, "Genus"]
missing_genus <- is.na(genus_labels) | genus_labels == ""
family_labels <- taxa_clean[, "Family"]
usable_family <- missing_genus & !is.na(family_labels) & family_labels != ""
genus_labels[usable_family] <- paste0("norank_f_", family_labels[usable_family])
genus_labels[missing_genus & !usable_family] <- "Unclassified_genus"
genus_counts <- rowsum(t(seqtab_clean), group = genus_labels, reorder = FALSE)
genus_counts <- t(genus_counts)
genus_out <- data.frame(donor_id = rownames(genus_counts), genus_counts, check.names = FALSE)
write.table(genus_out, file.path(results_dir, "genus_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

genus_rel <- genus_counts / rowSums(genus_counts)
genus_rel_out <- data.frame(donor_id = rownames(genus_rel), genus_rel, check.names = FALSE)
write.table(genus_rel_out, file.path(results_dir, "genus_relative_abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(
  list(
    seqtab_raw = seqtab,
    seqtab_nochim = seqtab_nochim,
    seqtab_clean = seqtab_clean,
    taxonomy = taxa,
    taxonomy_clean = taxa_clean,
    genus_counts = genus_counts,
    metadata = crosswalk,
    tracking = tracking
  ),
  file.path(results_dir, "wei2025_dada2_objects.rds")
)

summary <- list(
  started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
  ended = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  samples = nrow(seqtab_nochim),
  input_read_pairs = sum(tracking$reads.in),
  filtered_read_pairs = sum(tracking$reads.out),
  merged_read_pairs = sum(tracking$merged),
  nonchimeric_read_pairs = sum(tracking$nonchim),
  asvs_before_taxonomic_filter = ncol(seqtab_nochim),
  asvs_after_taxonomic_filter = ncol(seqtab_clean),
  genera = ncol(genus_counts),
  parameters = list(
    trim_left = c(20, 20), trunc_len = c(280, 230), max_ee = c(2, 5),
    trunc_q = 2, min_overlap = 20, pool = "pseudo", chimera = "consensus",
    parallel_backend = "BiocParallel::SerialParam"
  )
)
write_json(summary, file.path(results_dir, "dada2_run_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(provenance_dir, "dada2_sessionInfo.txt"))
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")
sink(type = "message")
sink(type = "output")
close(log_connection)
