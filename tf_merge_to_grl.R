#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage:\n  Rscript tf_merge_to_grl.R <csv_path> <download_root_dir> <out_rds>\n")
  quit(status = 1)
}

csv_path <- args[[1]]
root_dir <- args[[2]]
out_rds  <- args[[3]]

dt <- fread(csv_path)

dt <- dt[file_format == "bed" & !is.na(target) & !is.na(experiment_accession) & !is.na(href)]
dt[, local_path := file.path(root_dir, experiment_accession, basename(href))]
dt <- dt[file.exists(local_path)]
if (nrow(dt) == 0) stop("No local bed files found after filtering. Check csv/root_dir.")

# 只读 bed 的前三列，忽略后面所有列
read_bed3_as_gr <- function(path) {
  x <- fread(path, header = FALSE, sep = "\t", data.table = FALSE, fill = TRUE, quote = "")
  if (ncol(x) < 3) stop("BED has <3 columns: ", path)

  chr <- x[[1]]
  start0 <- suppressWarnings(as.integer(x[[2]]))
  end1   <- suppressWarnings(as.integer(x[[3]]))

  ok <- !is.na(chr) & !is.na(start0) & !is.na(end1)
  chr <- chr[ok]; start0 <- start0[ok]; end1 <- end1[ok]

  # BED: 0-based half-open [start, end); GRanges: 1-based closed [start+1, end]
  GRanges(seqnames = chr, ranges = IRanges(start = start0 + 1L, end = end1))
}

spl <- split(dt$local_path, dt$target)

grl <- GRangesList(lapply(spl, function(paths) {
  gr_all <- do.call(c, lapply(unique(paths), read_bed3_as_gr))
  reduce(gr_all, ignore.strand = TRUE)  # 一行 merge
}))
names(grl) <- names(spl)

saveRDS(grl, out_rds, compress = "xz")
cat(sprintf("Saved GRangesList: %d TFs -> %s\n", length(grl), out_rds))