# Host and virus reference constructors --------------------------------------

normalize_host_sequences <- function(df) {
  if (!is.data.frame(df) || !all(c("chr", "start", "end") %in% names(df))) {
    stop("Host sequences must be a data frame with chr, start, and end columns.", call. = FALSE)
  }
  out <- df[, c("chr", "start", "end"), drop = FALSE]
  out$chr <- trimws(gsub("(?i)^chr", "", as.character(out$chr), perl = TRUE))
  out$start <- suppressWarnings(as.numeric(out$start))
  out$end <- suppressWarnings(as.numeric(out$end))
  out <- out[!is.na(out$chr) & nzchar(out$chr) & !is.na(out$start) & !is.na(out$end), , drop = FALSE]
  if (!nrow(out) || anyDuplicated(out$chr) || any(out$end <= out$start)) {
    stop("Host sequences must have unique names and positive lengths.", call. = FALSE)
  }
  out
}

read_host_fai <- function(path) {
  x <- utils::read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         comment.char = "", quote = "")
  if (ncol(x) < 2L) stop("FAI file must contain sequence name and length columns.", call. = FALSE)
  normalize_host_sequences(data.frame(chr = x[[1]], start = 0, end = x[[2]], stringsAsFactors = FALSE))
}

read_host_sequences <- function(sequences) {
  if (is.data.frame(sequences)) return(normalize_host_sequences(sequences))
  if (!is.character(sequences) || length(sequences) != 1L || is.na(sequences)) {
    stop("sequences must be a chromosome table, a local file path, or a UCSC assembly.", call. = FALSE)
  }
  if (file.exists(sequences)) {
    if (tolower(tools::file_ext(sequences)) == "fai") return(read_host_fai(sequences))
    x <- utils::read.table(sequences, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           check.names = FALSE)
    return(normalize_host_sequences(x))
  }
  fetch_ucsc_chrom_sizes(sequences)
}

resolve_ucsc_assembly <- function(host) {
  key <- tolower(trimws(as.character(host)))
  aliases <- c(human = "hg38", homo_sapiens = "hg38", mouse = "mm39",
               mus_musculus = "mm39", rat = "rn7", rattus_norvegicus = "rn7",
               zebrafish = "danRer11", danio_rerio = "danRer11", fruitfly = "dm6",
               drosophila = "dm6", yeast = "sacCer3")
  if (length(key) != 1L || is.na(key) || !nzchar(key)) stop("assembly must be one non-empty string.", call. = FALSE)
  unname(if (key %in% names(aliases)) aliases[[key]] else host)
}

filter_primary_chromosomes <- function(df, host, assembly = NULL) {
  # UCSC download contains many alternate/haplotype sequences. Keep its primary assembly.
  if (tolower(resolve_ucsc_assembly(host)) == "hg38") {
    order_chr <- c(as.character(1:22), "X", "Y")
    out <- df[df$chr %in% order_chr, , drop = FALSE]
    idx <- match(order_chr, out$chr)
    return(out[idx[!is.na(idx)], , drop = FALSE])
  }
  primary <- df[!grepl("_|random|Un|alt|fix|hap", df$chr, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(primary)) stop("No primary UCSC sequences remained after filtering.", call. = FALSE)
  primary
}

fetch_ucsc_chrom_sizes <- function(host) {
  assembly <- resolve_ucsc_assembly(host)
  url <- sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes", assembly, assembly)
  fetched <- tryCatch(utils::read.table(url, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                                        comment.char = "", quote = ""), error = function(e) NULL)
  if (is.null(fetched) || ncol(fetched) < 2L) {
    stop("Unable to fetch UCSC chromosome sizes for assembly: ", assembly, call. = FALSE)
  }
  df <- normalize_host_sequences(data.frame(chr = fetched[[1]], start = 0, end = fetched[[2]], stringsAsFactors = FALSE))
  filter_primary_chromosomes(df, host = host, assembly = assembly)
}

host_sequence_colors <- function(sequences, colors = NULL) {
  ids <- as.character(sequences$chr)
  out <- stats::setNames(grDevices::hcl.colors(length(ids), palette = "Dynamic"), ids)
  if (is.null(colors)) return(out)
  if (is.null(names(colors)) || anyDuplicated(names(colors)) || any(!names(colors) %in% ids)) {
    stop("colors must be a named vector using host sequence names.", call. = FALSE)
  }
  out[names(colors)] <- unname(colors)
  out
}

#' Construct a host genome reference
#'
#' @param sequences A UCSC assembly name, local `.fai`/TSV path, or data frame
#'   with `chr`, `start`, and `end` columns.
#' @param annotation Optional local host annotation path, retained as source
#'   metadata for annotation-dependent tracks.
#' @param colors Optional named overrides for automatically generated sequence colors.
#' @return A `vi_host_genome` object.
#' @export
host_genome <- function(sequences = "hg38", annotation = NULL, colors = NULL) {
  if (!is.null(annotation) && (!is.character(annotation) || length(annotation) != 1L || !file.exists(annotation))) {
    stop("annotation must be an existing local file path.", call. = FALSE)
  }
  seqs <- read_host_sequences(sequences)
  local_file <- if (is.character(sequences) && length(sequences) == 1L && file.exists(sequences)) sequences else NULL
  remote <- is.null(local_file) && is.character(sequences) && length(sequences) == 1L && !is.data.frame(sequences)
  structure(list(
    sequences = seqs,
    genes = NULL,
    sources = list(sequence_file = local_file, annotation_file = annotation),
    metadata = list(source = if (remote) "UCSC" else "local",
                    assembly = if (remote) resolve_ucsc_assembly(sequences) else NULL),
    cache = list(gene_density = list()),
    colors = host_sequence_colors(seqs, colors)
  ), class = "vi_host_genome")
}

#' Construct a virus genome reference
#'
#' @param name Optional plotting sequence name. It must match virus_chr in integration records; defaults to the accession when omitted.
#' @param length Virus genome length in base pairs.
#' @param features Optional virus feature table or `vi_virus_features` object.
#' @param accession Optional NCBI nucleotide accession used to retrieve features and recorded as reference provenance.
#' @return A `vi_virus_genome` object.
#' @export
virus_genome <- function(name = NULL, length = NULL, features = NULL, accession = NULL) {
  if (!is.null(accession)) {
    if (is.null(features) || is.null(length)) {
      info <- virus_info(accession)
      if (is.null(features)) features <- virus_features(info)
      if (is.null(length)) length <- annotation_virus_length(info)
    }
    if (is.null(name)) name <- accession
  }
  if (is.null(length) && is.numeric(name) && length(name) == 1L && !is.null(names(name))) {
    length <- unname(name); name <- names(name)
  }
  name <- trimws(gsub("(?i)^chr", "", as.character(name), perl = TRUE))
  if (length(name) != 1L || is.na(name) || !nzchar(name)) stop("Virus name must be a single non-empty string.", call. = FALSE)
  if (!is.numeric(length) || length(length) != 1L || is.na(length) || length <= 0) stop("Virus length must be a single positive number.", call. = FALSE)
  feature_obj <- if (is.null(features)) NULL else virus_features(features)
  if (!is.null(feature_obj)) feature_obj$chr <- rep(name, nrow(feature_obj))
  structure(list(
    sequences = data.frame(chr = name, start = 0, end = as.numeric(length), stringsAsFactors = FALSE),
    name = name, length = as.numeric(length), features = feature_obj,
    metadata = list(accession = accession, source = if (is.null(accession)) "local" else "NCBI")
  ), class = "vi_virus_genome")
}

resolve_host_chrom_sizes <- function(host = NULL, chrom_file = NULL) {
  if (!is.null(host) && !is.null(chrom_file)) {
    stop("Please supply either host or chrom_file, not both.", call. = FALSE)
  }
  if (inherits(host, "vi_host_genome")) return(host$sequences)
  if (!is.null(chrom_file)) return(host_genome(sequences = chrom_file)$sequences)
  if (is.data.frame(host) || is.character(host)) return(host_genome(sequences = host)$sequences)
  stop("host must be a vi_host_genome object, a sequence table, or a reference path.", call. = FALSE)
}

validate_virus_info <- function(virus_name, virus_length, host_chr) {
  if (virus_name %in% host_chr) stop("virus_name conflicts with a host chromosome name: ", virus_name, call. = FALSE)
  list(virus_name = virus_name, virus_length = virus_length)
}




