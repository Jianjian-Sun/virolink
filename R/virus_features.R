as_vi_virus_features <- function(x = NULL, feature = NULL, start = NULL,
                                 end = NULL, type = NULL,
                                 allow_circular = FALSE) {
  if (is.null(x)) {
    x <- data.frame(
      feature = feature,
      start = start,
      end = end,
      stringsAsFactors = FALSE
    )
    if (!is.null(type)) {
      x$type <- type
    }
  }

  if (!is.data.frame(x)) {
    stop("x must be a data frame or NULL when feature/start/end are supplied.", call. = FALSE)
  }

  df <- as.data.frame(x, stringsAsFactors = FALSE)
  rename_one <- function(target, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0) {
      return(NULL)
    }
    df[[target]] <<- df[[hit[1]]]
    invisible(NULL)
  }

  rename_one("feature", c("feature", "gene", "name", "label"))
  rename_one("start", c("start", "start_pos", "from"))
  rename_one("end", c("end", "stop", "end_pos", "to"))
  rename_one("type", c("type", "feature_type"))

  required_cols <- c("feature", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Virus features must contain columns: feature, start, end", call. = FALSE)
  }

  if (!"type" %in% colnames(df)) {
    df$type <- "feature"
  }

  out <- df[, c("feature", "start", "end", "type"), drop = FALSE]
  out$feature <- as.character(out$feature)
  out$start <- suppressWarnings(as.numeric(out$start))
  out$end <- suppressWarnings(as.numeric(out$end))
  out$type <- as.character(out$type)

  keep <- !is.na(out$feature) & nzchar(out$feature) &
    !is.na(out$start) & !is.na(out$end)
  if (isTRUE(allow_circular)) {
    keep <- keep & out$start > 0 & out$end > 0 & out$end != out$start
  } else {
    keep <- keep & out$end > out$start
  }
  if (!all(keep)) {
    warning(sum(!keep), " virus features were dropped because they were incomplete or invalid.")
  }
  out <- out[keep, , drop = FALSE]

  if (nrow(out) == 0) {
    stop("No valid virus features were found.", call. = FALSE)
  }

  out <- out[order(out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("vi_virus_features", "data.frame")
  out
}

#' Built-in HBV feature annotations
#'
#' @param version HBV annotation version. Currently only \code{"ayw"} is supported.
#' @return A \code{vi_virus_features} object.
#' @export
hbv_features <- function(version = "ayw") {
  version <- match.arg(version, "ayw")
  as_vi_virus_features(data.frame(
    feature = c("preS1", "preS2", "S", "X", "preC", "C", "P", "P"),
    start = c(2850, 1, 155, 1374, 1814, 1901, 2307, 1),
    end = c(3204, 154, 835, 1838, 1900, 2452, 3215, 1623),
    type = c("surface", "surface", "surface", "regulatory", "core", "core", "polymerase", "polymerase"),
    stringsAsFactors = FALSE
  ))
}

#' Download virus annotation information from NCBI
#'
#' Fetches the complete GFF3 annotation for a nucleotide accession. Use
#' [virus_features()] to convert selected records into the standard object used
#' by virolink virus tracks. This uses the public NCBI E-utilities URL directly,
#' so no additional R package is required for downloading the annotation.
#'
#' @param accession NCBI nucleotide accession, for example `NC_001526.4`.
#' @param out_file Optional path where the downloaded GFF3 text is saved.
#' @param timeout Timeout in seconds for the URL backend.
#' @return A list containing `accession`, `gff3_url`, raw `gff3` text,
#'   `gff3_file`, and `metadata`.
#' @export
virus_info <- function(accession,
                       out_file = NULL,
                       timeout = 60) {
  if (!is.null(out_file) && (!is.character(out_file) || length(out_file) != 1L ||
                             is.na(out_file) || !nzchar(out_file))) {
    stop("out_file must be NULL or one non-empty file path.", call. = FALSE)
  }
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    stop("timeout must be one positive number.", call. = FALSE)
  }

  url <- ncbi_gff3_url(accession)
  gff3 <- fetch_ncbi_gff3_url(url, timeout = timeout)

  if (!is.null(out_file)) {
    writeLines(gff3, out_file, useBytes = TRUE)
  }

  list(
    accession = accession,
    gff3_url = url,
    gff3 = gff3,
    gff3_file = out_file,
    metadata = parse_gff3_metadata(gff3)
  )
}

#' Create or extract virus feature annotations
#'
#' When `annotation` is a result from [virus_info()], this parses the downloaded
#' GFF3 and returns selected feature records for plotting. When `annotation` is a
#' data frame, or when `feature`, `start`, and `end` are supplied directly, this
#' standardizes those records into the same `vi_virus_features` object.
#'
#' @param annotation A list returned by [virus_info()], a data frame with feature
#'   intervals, or `NULL` when columns are supplied separately.
#' @param type Optional GFF3 feature type filter for [virus_info()] results, or
#'   optional feature type labels when columns are supplied directly.
#' @param feature Optional feature name filter for [virus_info()] results, or
#'   feature names when columns are supplied directly.
#' @param start Feature start positions when columns are supplied directly.
#' @param end Feature end positions when columns are supplied directly.
#' @return A `vi_virus_features` object.
#' @export
virus_features <- function(annotation = NULL, type = NULL, feature = NULL,
                           start = NULL, end = NULL) {
  if (!is.list(annotation) || is.null(annotation$gff3) ||
      !is.character(annotation$gff3) || length(annotation$gff3) != 1L) {
    return(as_vi_virus_features(annotation, feature = feature, start = start,
                                end = end, type = type))
  }

  virus_length <- annotation_virus_length(annotation)
  parse_ncbi_gff3(
    annotation$gff3,
    annotation$accession,
    virus_length = virus_length,
    type = type,
    feature = feature
  )
}

#' Construct the NCBI E-utilities URL for a nucleotide GFF3 annotation
#'
#' @param accession NCBI nucleotide accession.
#' @return A length-one character URL.
#' @keywords internal
ncbi_gff3_url <- function(accession) {
  paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "?db=nuccore&id=", utils::URLencode(accession, reserved = TRUE),
    "&rettype=gff3&retmode=text"
  )
}

fetch_ncbi_gff3_url <- function(url, timeout) {
  target <- tempfile(fileext = ".gff3")
  on.exit(unlink(target), add = TRUE)
  old_timeout <- getOption("timeout")
  options(timeout = timeout)
  on.exit(options(timeout = old_timeout), add = TRUE)

  utils::download.file(url, target, quiet = TRUE, mode = "wb")
  paste(readLines(target, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

parse_ncbi_gff3 <- function(gff3, accession = NULL, virus_length = NULL,
                            type = NULL, feature = NULL) {
  if (!is.character(gff3) || length(gff3) != 1L || is.na(gff3)) {
    stop("gff3 must be one character string.", call. = FALSE)
  }
  lines <- strsplit(gff3, "\n", fixed = TRUE)[[1L]]
  lines <- sub("\r$", "", lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0L) {
    stop("The GFF3 annotation contains no feature records.", call. = FALSE)
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  valid_fields <- vapply(fields, length, integer(1L)) >= 9L
  if (!all(valid_fields)) {
    warning(sum(!valid_fields), " malformed GFF3 records were omitted.")
    fields <- fields[valid_fields]
  }
  if (length(fields) == 0L) {
    stop("The GFF3 annotation contains no valid feature records.", call. = FALSE)
  }

  table <- do.call(rbind, lapply(fields, function(x) {
    c(seqid = x[1L], source = x[2L], type = x[3L], start = x[4L],
      end = x[5L], score = x[6L], strand = x[7L], phase = x[8L],
      attributes = paste(x[9:length(x)], collapse = "\t"))
  }))
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  if (!is.null(accession)) {
    table <- table[table$seqid %in% c(accession, sub("\\..*$", "", accession)), , drop = FALSE]
  }
  if (nrow(table) == 0L) {
    stop("No feature records were found for the requested accession.", call. = FALSE)
  }

  attrs <- lapply(table$attributes, parse_gff3_attributes)
  name_from_attrs <- function(x) {
    for (key in c("gene", "locus_tag", "Name", "product", "ID")) {
      if (!is.null(x[[key]]) && nzchar(x[[key]])) return(x[[key]])
    }
    NA_character_
  }
  labels <- vapply(attrs, name_from_attrs, character(1L))
  starts <- suppressWarnings(as.numeric(table$start))
  ends <- suppressWarnings(as.numeric(table$end))
  out <- data.frame(
    feature = labels,
    start = starts,
    end = ends,
    type = table$type,
    stringsAsFactors = FALSE
  )

  if (!is.null(type)) {
    out <- out[out$type %in% type, , drop = FALSE]
  }
  if (!is.null(feature)) {
    out <- out[out$feature %in% feature, , drop = FALSE]
  }
  if (nrow(out) == 0L) {
    stop("No virus features matched the requested filters.", call. = FALSE)
  }

  keep <- !is.na(out$feature) & nzchar(out$feature) &
    !is.na(out$start) & !is.na(out$end)
  if (!all(keep)) {
    warning(sum(!keep), " GFF3 features were dropped because their name or interval was invalid.")
  }
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0L) stop("No valid named virus features were found.", call. = FALSE)

  if (!is.null(virus_length)) {
    out <- expand_circular_virus_features(out, virus_length)
  }

  as_vi_virus_features(out)
}

annotation_virus_length <- function(annotation) {
  if (!is.list(annotation)) {
    return(NULL)
  }

  if (is.data.frame(annotation$metadata) && nrow(annotation$metadata) > 0L &&
      "end" %in% colnames(annotation$metadata)) {
    virus_length <- suppressWarnings(max(as.numeric(annotation$metadata$end), na.rm = TRUE))
    if (is.finite(virus_length) && virus_length > 0) {
      return(virus_length)
    }
  }

  if (!is.null(annotation$gff3) && is.character(annotation$gff3) && length(annotation$gff3) == 1L) {
    metadata <- parse_gff3_metadata(annotation$gff3)
    if (nrow(metadata) > 0L && "end" %in% colnames(metadata)) {
      virus_length <- suppressWarnings(max(as.numeric(metadata$end), na.rm = TRUE))
      if (is.finite(virus_length) && virus_length > 0) {
        return(virus_length)
      }
    }
  }

  NULL
}

expand_circular_virus_features <- function(df, virus_length) {
  virus_length <- suppressWarnings(as.numeric(virus_length))
  if (length(virus_length) != 1L || is.na(virus_length) || virus_length <= 0) {
    return(df)
  }

  out <- vector("list", nrow(df))
  idx <- 0L
  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    if (!is.na(row$start) && !is.na(row$end) && row$end < row$start) {
      idx <- idx + 1L
      out[[idx]] <- transform(row, end = virus_length)
      idx <- idx + 1L
      out[[idx]] <- transform(row, start = 1, end = row$end)
    } else if (!is.na(row$end) && row$end > virus_length) {
      idx <- idx + 1L
      out[[idx]] <- transform(row, end = virus_length)
      idx <- idx + 1L
      out[[idx]] <- transform(row, start = 1, end = row$end - virus_length)
    } else {
      idx <- idx + 1L
      out[[idx]] <- row
    }
  }

  do.call(rbind, out[seq_len(idx)])
}

parse_gff3_attributes <- function(value) {
  parts <- strsplit(value, ";", fixed = TRUE)[[1L]]
  result <- list()
  for (part in parts) {
    pair <- strsplit(part, "=", fixed = TRUE)[[1L]]
    if (length(pair) < 2L) next
    key <- trimws(pair[1L])
    val <- utils::URLdecode(paste(pair[-1L], collapse = "="))
    result[[key]] <- val
  }
  result
}

parse_gff3_metadata <- function(gff3) {
  lines <- strsplit(gff3, "\n", fixed = TRUE)[[1L]]
  lines <- sub("\r$", "", lines)
  regions <- lines[startsWith(lines, "##sequence-region")]
  if (length(regions) == 0L) return(data.frame())
  parts <- strsplit(regions, "[[:space:]]+")
  data.frame(
    seqid = vapply(parts, `[`, character(1L), 2L),
    start = as.numeric(vapply(parts, `[`, character(1L), 3L)),
    end = as.numeric(vapply(parts, `[`, character(1L), 4L)),
    stringsAsFactors = FALSE
  )
}

normalize_virus_features <- function(features, virus_length = NULL) {
  allow_circular <- !is.null(virus_length)
  if (inherits(features, "vi_virus_features")) {
    df <- as_vi_virus_features(features, allow_circular = allow_circular)
  } else {
    df <- as_vi_virus_features(features, allow_circular = allow_circular)
  }

  if (!is.null(virus_length)) {
    df <- expand_circular_virus_features(df, virus_length)
    return(as_vi_virus_features(df))
  }

  df
}

layout_virus_features <- function(features, virus_length = NULL) {
  df <- as.data.frame(normalize_virus_features(features), stringsAsFactors = FALSE)
  if (!is.null(virus_length)) {
    wrap <- df$end <= df$start
    df$end[wrap] <- df$end[wrap] + virus_length
  }

  df$level <- 1L
  level_ends <- numeric()
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    for (level in seq_along(level_ends)) {
      if (df$start[i] >= level_ends[level]) {
        df$level[i] <- level
        level_ends[level] <- df$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      level_ends <- c(level_ends, df$end[i])
      df$level[i] <- length(level_ends)
    }
  }

  df
}
