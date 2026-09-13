#' Create a host gene-density track
#'
#' @param bin_width Width of host genomic bins in base pairs.
#' @param height Track height.
#' @param fill Fill color for density bars.
#' @param gene_types Optional gene-type filter.
#' @param label Optional track label.
#' @param label_cex Track-label text size.
#' @param label_col Track-label color.
#' @return A `vi_track` object.
#' @export
track_host_gene_density <- function(bin_width = 1e6, height = 0.12,
                                    fill = "#B2182B", gene_types = NULL,
                                    label = NULL, label_cex = 0.8,
                                    label_col = "grey30") {
  if (!is.numeric(bin_width) || length(bin_width) != 1L || is.na(bin_width) || bin_width <= 0) {
    stop("bin_width must be a single positive number.", call. = FALSE)
  }
  structure(list(type = "host_gene_density", bin_width = as.numeric(bin_width),
                 height = height, fill = fill, gene_types = gene_types,
                 label = label, label_cex = label_cex, label_col = label_col),
            class = "vi_track")
}

parse_host_attributes <- function(x, gff = FALSE) {
  fields <- strsplit(x, ";", fixed = TRUE)[[1L]]
  out <- list()
  for (field in trimws(fields)) {
    if (!nzchar(field)) next
    if (gff) {
      pair <- strsplit(field, "=", fixed = TRUE)[[1L]]
      if (length(pair) >= 2L) out[[pair[1L]]] <- utils::URLdecode(paste(pair[-1L], collapse = "="))
    } else {
      pair <- strsplit(field, "[[:space:]]+", perl = TRUE)[[1L]]
      if (length(pair) >= 2L) out[[pair[1L]]] <- gsub('^"|"$', "", paste(pair[-1L], collapse = " "))
    }
  }
  out
}

read_host_genes <- function(path, sequences) {
  ext <- tolower(tools::file_ext(path))
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else path
  raw <- utils::read.table(con, sep = "\t", header = FALSE, quote = "", comment.char = "",
                           stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(raw) < 9L) stop("Host annotation must be a GTF or GFF3 file with nine columns.", call. = FALSE)
  raw <- raw[raw[[3L]] == "gene", , drop = FALSE]
  if (!nrow(raw)) stop("Host annotation contains no gene records.", call. = FALSE)
  gff <- ext %in% c("gff", "gff3") || grepl("\\.gff3?\\.gz$", path, ignore.case = TRUE)
  attrs <- lapply(raw[[9L]], parse_host_attributes, gff = gff)
  field <- function(keys) vapply(attrs, function(x) {
    for (key in keys) if (!is.null(x[[key]]) && nzchar(x[[key]])) return(x[[key]])
    NA_character_
  }, character(1L))
  out <- data.frame(chr = trimws(gsub("(?i)^chr", "", raw[[1L]], perl = TRUE)),
                    start = as.numeric(raw[[4L]]), end = as.numeric(raw[[5L]]),
                    strand = as.character(raw[[7L]]),
                    gene_id = field(c("gene_id", "ID", "locus_tag")),
                    gene_name = field(c("gene_name", "Name", "gene", "gene_id", "ID")),
                    gene_type = field(c("gene_type", "gene_biotype", "biotype")),
                    stringsAsFactors = FALSE)
  out <- out[out$chr %in% sequences$chr & !is.na(out$start) & !is.na(out$end) & out$end >= out$start, , drop = FALSE]
  if (!nrow(out)) stop("No annotation genes matched host sequences.", call. = FALSE)
  out
}

fetch_ucsc_refseq_genes <- function(assembly, sequences) {
  url <- sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/database/ncbiRefSeq.txt.gz", assembly)
  target <- tempfile(fileext = ".txt.gz")
  on.exit(unlink(target), add = TRUE)
  utils::download.file(url, target, quiet = TRUE, mode = "wb")
  con <- gzfile(target, open = "rt")
  on.exit(close(con), add = TRUE)
  raw <- utils::read.table(con, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  if (ncol(raw) < 13L) stop("Unable to parse UCSC ncbiRefSeq annotation.", call. = FALSE)
  out <- data.frame(chr = trimws(gsub("(?i)^chr", "", raw[[3L]], perl = TRUE)),
                    start = as.numeric(raw[[5L]]) + 1, end = as.numeric(raw[[6L]]),
                    strand = as.character(raw[[4L]]), gene_id = as.character(raw[[2L]]),
                    gene_name = as.character(raw[[13L]]), gene_type = NA_character_, stringsAsFactors = FALSE)
  out <- out[out$chr %in% sequences$chr & !is.na(out$start) & !is.na(out$end) & out$end >= out$start, , drop = FALSE]
  if (!nrow(out)) stop("No UCSC ncbiRefSeq genes matched host sequences.", call. = FALSE)
  out
}

ensure_host_genes <- function(host) {
  if (!inherits(host, "vi_host_genome")) stop("host must be a vi_host_genome object.", call. = FALSE)
  if (!is.null(host$genes)) return(host)
  if (!is.null(host$sources$annotation_file)) {
    host$genes <- read_host_genes(host$sources$annotation_file, host$sequences)
  } else if (identical(host$metadata$source, "UCSC") && !is.null(host$metadata$assembly)) {
    host$genes <- fetch_ucsc_refseq_genes(host$metadata$assembly, host$sequences)
  } else {
    stop("host gene density requires a local GTF/GFF3 annotation or a UCSC host assembly.", call. = FALSE)
  }
  host
}

make_host_gene_density <- function(host, bin_width, gene_types = NULL) {
  host <- ensure_host_genes(host)
  genes <- host$genes
  if (!is.null(gene_types)) genes <- genes[genes$gene_type %in% gene_types, , drop = FALSE]
  out <- lapply(seq_len(nrow(host$sequences)), function(i) {
    seq_row <- host$sequences[i, , drop = FALSE]
    starts <- seq(seq_row$start, seq_row$end - 1, by = bin_width)
    ends <- pmin(starts + bin_width, seq_row$end)
    count <- integer(length(starts))
    g <- genes[genes$chr == seq_row$chr, , drop = FALSE]
    if (nrow(g)) {
      idx <- findInterval(g$start - seq_row$start, starts - seq_row$start, rightmost.closed = TRUE)
      idx <- pmax(1L, pmin(idx, length(count)))
      count <- tabulate(idx, nbins = length(count))
    }
    data.frame(chr = seq_row$chr, start = starts, end = ends, count = count, stringsAsFactors = FALSE)
  })
  list(host = host, data = do.call(rbind, out))
}

prepare_host_reference_tracks <- function(tracks, host) {
  out <- tracks
  for (i in seq_along(out)) {
    if (identical(out[[i]]$type, "host_gene_density")) {
      density <- make_host_gene_density(host, out[[i]]$bin_width, out[[i]]$gene_types)
      host <- density$host
      out[[i]]$data <- density$data
    }
  }
  list(host = host, tracks = out)
}

draw_host_gene_density <- function(track, cfg) {
  data <- track$data
  max_count <- max(data$count, na.rm = TRUE)
  if (!is.finite(max_count) || max_count < 1) max_count <- 1
  circlize::circos.genomicTrackPlotRegion(
    data, ylim = c(0, max_count), track.height = track$height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      if (identical(circlize::CELL_META$sector.index, cfg$virus_name)) return(invisible(NULL))
      circlize::circos.genomicRect(region, value, ybottom = 0, ytop = value$count,
                                   col = track$fill, border = NA)
      if (!is.null(track$label)) draw_track_label(track$label, y = max_count * 0.8,
                                                   cex = track$label_cex, col = track$label_col)
    }
  )
  invisible(NULL)
}



