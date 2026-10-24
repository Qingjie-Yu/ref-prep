suppressPackageStartupMessages({
  library(glue)
  library(S4Vectors)
  library(data.table)
  library(IRanges)
  library(GenomicRanges)
})


BED_STD_COLS <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")

#' BEDn+m -> GRanges
bed_to_gr <- function(bed_path, extra_cols = NULL) {
  df <- fread(bed_path, header=FALSE, sep="\t", data.table = FALSE)
  nc <- ncol(df)
  n <- if (!is.null(extra_cols)) nc - length(extra_cols) else nc
  stopifnot(n >= 3L, n <= 12L)
  ne <- nc - n

  std <- BED_STD_COLS[seq_len(n)]
  names(df) <- c(std, extra_cols)

  strand <- if (n >= 6) {
    s <- as.character(df$strand); s[!s %in% c("+", "-")] <- "*"
    s
  } else "*"

  gr <- GRanges(
    seqnames = df$chrom,
    ranges = IRanges(start = df$chromStart + 1L, end = df$chromEnd),
    strand = strand
  )
  mc <- df[, setdiff(names(df), c("chrom", "chromStart", "chromEnd", "strand")), drop=FALSE]
  if ("thickStart" %in% names(mc)) mc$thickStart <- mc$thickStart + 1L
  mcols(gr) <- mc
  gr
}

#' GRanges -> BEDn+m
gr_to_bed <- function(gr, bed_path, extra_cols = NULL){
  mc <- mcols(gr)
  cn <- colnames(mc)

  opt <- BED_STD_COLS[4:12]
  present <- opt %in% cn | opt == "strand"
  n <- 3L + if (all(present)) 9L else which(!present)[1] - 1L
  std <- BED_STD_COLS[seq_len(n)]
  message(sprintf("BED%d+%d", n, length(extra_cols)))

  bed <- data.frame(
    chrom = as.character(seqnames(gr)),
    chromStart = start(gr) - 1L,
    chromEnd = end(gr)
  )
  for (col in std[-(1:3)]) {
    bed[[col]] <- switch(col,
    strand = { s <- as.character(strand(gr)); s[s == "*"] <- "."; s },
    thickStart = as.integer(mc$thickStart) - 1L,
    as.vector(mc[[col]])
    )
  }
  for (col in extra_cols) bed[[col]] <- as.vector(mc[[col]])
  fwrite(bed, bed_path, sep = "\t", col.names = FALSE, quote = FALSE, na = ".")
  invisible(bed)
}

#' Quick check .rds file
view_rds <- function(file_path) {
  obj <- readRDS(file_path)
  try(
    print(head(obj))
  )
  try(
    print(length(obj))
  )
  return(obj)
}

#' Quick check overleap
check_overlap <- function(gr1, gr2) {
  gene_hits <- findOverlaps(gr1, gr2)
  if (length(gene_hits) > 0) {
    cat(sprintf("Found %d overlapping pairs\n", length(gene_hits)))
  } else {
    cat("No overlapping found\n")
  }
}

check_overlap_self <- function(gr) {
  hits <- GenomicRanges::findOverlaps(gr, gr)
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  if (length(hits) > 0) {
    cat(sprintf(
      "Found %d internal overlapping interval pairs\n",
      length(hits)
    ))
  } else {
    cat("No internal overlaps found\n")
  }
}
