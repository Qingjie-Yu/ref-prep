#!/usr/bin/env Rscript
# tf_merge_to_grl.R — step 4: merge ENCODE peak files per TF into a GRangesList
#
# Usage:
#   Rscript tf_merge_to_grl.R --manifest ref_mm10/manifest.tsv --outprefix ref_mm10/mm10_TF \
#       [--min-support 1] [--std-chrom TRUE] [--write-bed FALSE] [--genome mm10]
#
# For each target: read every peak file (BED, first 3 cols), concatenate, reduce()
# into non-overlapping regions, and record how many distinct experiments support
# each merged region. --min-support k keeps only regions seen in >= k experiments
# (k = 1 reproduces the old "union" behaviour).
#
# Outputs:
#   <outprefix>.grl.rds           GRangesList, one element per TF; mcols: n_exp, n_peaks
#   <outprefix>.merge_stats.tsv   per-TF peak counts / widths before & after merging
#   <outprefix>_bed/<TF>.bed      (optional) one BED per TF

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

`%||%` <- function(x, y) if (is.null(x)) y else x   # base R < 4.4 lacks it

# ---- args (--key value or --key=value; no extra packages needed) -------------
parse_args <- function(x) {
  out <- list(); i <- 1
  while (i <= length(x)) {
    k <- x[[i]]
    if (!startsWith(k, "--")) stop("unexpected argument: ", k)
    k <- sub("^--", "", k)
    if (grepl("=", k)) { kv <- strsplit(k, "=", fixed = TRUE)[[1]]; out[[kv[1]]] <- kv[2]; i <- i + 1 }
    else { out[[k]] <- x[[i + 1]]; i <- i + 2 }
  }
  out
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(a$manifest) || is.null(a$outprefix)) {
  cat("Usage: Rscript tf_merge_to_grl.R --manifest manifest.tsv --outprefix out/prefix",
      "[--min-support 1] [--std-chrom TRUE] [--write-bed FALSE] [--genome mm10]\n")
  quit(status = 1)
}
min_support <- as.integer(a[["min-support"]] %||% 1L)
std_chrom   <- as.logical(a[["std-chrom"]] %||% TRUE)
write_bed   <- as.logical(a[["write-bed"]] %||% FALSE)
genome      <- a$genome %||% NA_character_
dir.create(dirname(a$outprefix), showWarnings = FALSE, recursive = TRUE)

# ---- manifest -----------------------------------------------------------------
dt <- fread(a$manifest, sep = "\t", colClasses = "character")
dt <- dt[file_format == "bed" & nzchar(target) & nzchar(local_path)]
missing <- dt[!file.exists(local_path)]
if (nrow(missing)) warning(nrow(missing), " files in manifest are not on disk; skipped")
dt <- dt[file.exists(local_path)]
if (!nrow(dt)) stop("No local BED files found. Did you run extract with --download?")
message(sprintf("%d files, %d experiments, %d targets",
                nrow(dt), uniqueN(dt$experiment_accession), uniqueN(dt$target)))

# ---- read BED3 -> GRanges (fread reads .gz directly if R.utils is installed) ---
read_bed3 <- function(path, exp_acc) {
  x <- fread(path, header = FALSE, sep = "\t", select = 1:3, fill = TRUE, quote = "",
             colClasses = "character")
  s0 <- suppressWarnings(as.integer(x$V2)); e1 <- suppressWarnings(as.integer(x$V3))
  ok <- !is.na(s0) & !is.na(e1)                 # drops track/browser/header lines
  # BED [start, end) 0-based  ->  GRanges [start+1, end] 1-based
  GRanges(x$V1[ok], IRanges(s0[ok] + 1L, e1[ok]), exp = rep(exp_acc, sum(ok)))
}

std_re <- "^chr([0-9]+|X|Y)$"   # drop chrM, random, Un, alt contigs

merge_one <- function(sub) {
  gr <- do.call(c, Map(read_bed3, sub$local_path, sub$experiment_accession, USE.NAMES = FALSE))
  if (std_chrom) gr <- gr[grepl(std_re, as.character(seqnames(gr)))]
  n_raw <- length(gr); w_raw <- width(gr)

  red <- reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
  exps <- extractList(gr$exp, red$revmap)
  red$n_exp   <- lengths(unique(exps))       # distinct experiments supporting region
  red$n_peaks <- lengths(red$revmap)         # raw peaks collapsed into region
  red$revmap  <- NULL
  n_merged_all <- length(red)
  red <- red[red$n_exp >= min_support]

  list(gr = red, stats = data.table(
    target = sub$target[1],
    n_files = nrow(sub),
    n_experiments = uniqueN(sub$experiment_accession),
    n_peaks_raw = n_raw,
    median_width_raw = if (n_raw) median(w_raw) else NA_real_,
    n_regions_merged = n_merged_all,
    n_regions_kept = length(red),
    median_width_kept = if (length(red)) median(width(red)) else NA_real_,
    bp_covered = sum(as.numeric(width(red))),
    frac_regions_multi_exp = if (n_merged_all) mean(red$n_exp > 1) else NA_real_
  ))
}

tfs <- sort(unique(dt$target))
res <- vector("list", length(tfs)); names(res) <- tfs
for (i in seq_along(tfs)) {
  res[[i]] <- merge_one(dt[target == tfs[i]])
  if (i %% 25 == 0 || i == length(tfs)) message(sprintf("  merged %d/%d", i, length(tfs)))
}

grl <- GRangesList(lapply(res, `[[`, "gr"))
if (!is.na(genome)) genome(grl) <- genome
stats <- rbindlist(lapply(res, `[[`, "stats"))

# drop TFs left empty by --min-support
empty <- names(grl)[lengths(grl) == 0]
if (length(empty)) {
  message(length(empty), " targets have no regions at min-support=", min_support, "; dropped")
  grl <- grl[lengths(grl) > 0]
}

saveRDS(grl, paste0(a$outprefix, ".grl.rds"), compress = "xz")
fwrite(stats, paste0(a$outprefix, ".merge_stats.tsv"), sep = "\t")

if (write_bed) {
  bed_dir <- paste0(a$outprefix, "_bed"); dir.create(bed_dir, showWarnings = FALSE)
  for (tf in names(grl)) {
    g <- grl[[tf]]
    fwrite(data.table(as.character(seqnames(g)), start(g) - 1L, end(g), tf, g$n_exp),
           file.path(bed_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", tf), ".bed")),
           sep = "\t", col.names = FALSE)
  }
}

message(sprintf("Saved GRangesList: %d TFs, %s regions -> %s.grl.rds",
                length(grl), format(sum(lengths(grl)), big.mark = ","), a$outprefix))
