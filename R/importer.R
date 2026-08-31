#' Read integration caller output
#'
#' Import output from supported upstream viral-integration callers and convert
#' it to the package-standard `vi_integrations` table.
#'
#' @param path Path to an upstream caller output file.
#' @param format Caller output format. Use `"auto"` to detect supported bundled
#'   formats from column names.
#' @param sample Optional sample label. When `NULL`, no sample column is added
#'   and no label is inferred from the file name or path.
#' @param sep Field separator for text files.
#' @param ... Additional arguments passed to `utils::read.table()` for text
#'   files or `readxl::read_excel()` for Excel workbooks.
#' @return A `vi_integrations` object.
#' @export
read_caller_integrations <- function(path,
                                     format = c("auto", "ctat_vif", "virusfinder2", "virus_clip", "bsvf"),
                                     sample = NULL,
                                     sep = "\t",
                                     ...) {
  format <- match.arg(format)
  df <- read_import_table(path, sep = sep, ...)

  if (format == "auto") {
    format <- detect_integration_caller_format(df)
  }

  switch(
    format,
    ctat_vif = ctat_vif_to_integrations(df, path = path, sample = sample),
    virusfinder2 = virusfinder2_to_integrations(df, path = path, sample = sample),
    virus_clip = virus_clip_to_integrations(df, path = path, sample = sample),
    bsvf = bsvf_to_integrations(df, path = path, sample = sample)
  )
}

#' Read CTAT-VIF integration candidates
#'
#' CTAT-VIF columns are mapped as follows: `chrA`/`coordA` and `chrB`/`coordB`
#' are used to identify host and virus breakpoints, `total` is used as the
#' standard `support_reads` value, `contig` is retained as `raw_record_id`,
#' `prelim.primary_brkpt_type` as `breakpoint_type`, and `prelim.total`,
#' `split`, `span`, and `total` are retained as caller-specific support fields.
#'
#' @param path Path to a CTAT-VIF `insertion_candidates.tsv` file.
#' @param sample Optional sample label. When `NULL`, no sample column is added
#'   and no label is inferred from the file name or path.
#' @param sep Field separator.
#' @param ... Additional arguments passed to `utils::read.table()`.
#' @return A `vi_integrations` object.
#' @export
read_ctat_vif_integrations <- function(path, sample = NULL, sep = "\t", ...) {
  df <- read_import_table(path, sep = sep, ...)
  ctat_vif_to_integrations(df, path = path, sample = sample)
}

#' Read VirusFinder2 integration sites
#'
#' @param path Path to a VirusFinder2 integration-sites file.
#' @param sample Optional sample label. When `NULL`, no sample column is added
#'   and no label is inferred from the file name or path.
#' @param sep Field separator for text files.
#' @param ... Additional arguments passed to `utils::read.table()` for text
#'   files or `readxl::read_excel()` for Excel workbooks.
#' @return A `vi_integrations` object.
#' @export
read_virusfinder2_integrations <- function(path, sample = NULL, sep = "\t", ...) {
  df <- read_import_table(path, sep = sep, ...)
  virusfinder2_to_integrations(df, path = path, sample = sample)
}

#' Read Virus-Clip integration sites
#'
#' Virus-Clip columns are mapped using the paired chromosome, position, strand,
#' and supporting-read fields. When present, `confidence` is preserved as an
#' extra column on the standardized `vi_integrations` table.
#'
#' @param path Path to a Virus-Clip output file.
#' @param sample Optional sample label. When `NULL`, no sample column is added
#'   and no label is inferred from the file name or path.
#' @param sep Field separator for text files.
#' @param ... Additional arguments passed to `utils::read.table()` for text
#'   files or `readxl::read_excel()` for Excel workbooks.
#' @return A `vi_integrations` object.
#' @export
read_virus_clip_integrations <- function(path, sample = NULL, sep = "\t", ...) {
  df <- read_import_table(path, sep = sep, ...)
  virus_clip_to_integrations(df, path = path, sample = sample)
}

#' Read BSVF integration sites
#'
#' BSVF columns are mapped using the host chromosome and breakpoint together
#' with the virus interval and support count. Because BSVF reports a virus
#' range rather than a single virus coordinate, the standardized `virus_pos`
#' is set to the midpoint of `virus-start` and `virus-end`. The original virus
#' interval and cluster name are preserved as extra columns.
#'
#' @param path Path to a BSVF output file.
#' @param sample Optional sample label. When `NULL`, no sample column is added
#'   and no label is inferred from the file name or path.
#' @param sep Field separator for text files.
#' @param ... Additional arguments passed to `utils::read.table()` for text
#'   files or `readxl::read_excel()` for Excel workbooks.
#' @return A `vi_integrations` object.
#' @export
read_bsvf_integrations <- function(path, sample = NULL, sep = "\t", ...) {
  df <- read_import_table(path, sep = sep, ...)
  bsvf_to_integrations(df, path = path, sample = sample)
}

read_import_table <- function(path, sep = "\t", ...) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  if (is_excel_like_file(path)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Reading Excel caller outputs requires the readxl package.",
        call. = FALSE
      )
    }
    return(as.data.frame(readxl::read_excel(path, ...), stringsAsFactors = FALSE))
  }

  utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    ...
  )
}

is_excel_like_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xls", "xlsx")) {
    return(TRUE)
  }

  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 4)
  length(magic) >= 2L && identical(toupper(paste(format(magic[1:2]), collapse = "")), "504B")
}

detect_integration_caller_format <- function(df) {
  nms <- colnames(df)
  if (all(c("chrA", "coordA", "orientA", "chrB", "coordB", "orientB") %in% nms)) {
    return("ctat_vif")
  }
  if (all(c(
    "Chromosome 1", "Position 1", "Strand 1",
    "Chromosome 2", "Position 2", "Strand 2"
  ) %in% nms)) {
    if ("confidence" %in% nms) {
      return("virus_clip")
    }
    return("virusfinder2")
  }
  if (all(c(
    "Chr", "breakpoint", "virus-start", "virus-end",
    "virusstrand", "how-many-reads-support", "cluster-name"
  ) %in% nms)) {
    return("bsvf")
  }

  stop(
    "Could not detect caller format from column names. ",
    "Please set format to one of: ctat_vif, virusfinder2, virus_clip, bsvf.",
    call. = FALSE
  )
}

ctat_vif_to_integrations <- function(df, path, sample = NULL) {
  ensure_import_columns(
    df,
    c(
      "contig", "chrA", "coordA", "orientA", "chrB", "coordB", "orientB",
      "prelim.primary_brkpt_type", "prelim.total", "split", "span", "total"
    ),
    format = "CTAT-VIF"
  )

  support_col <- if ("total" %in% colnames(df)) {
    "total"
  } else if ("prelim.total" %in% colnames(df)) {
    "prelim.total"
  } else {
    NULL
  }

  standardize_paired_breakpoints(
    df = df,
    chr1 = "chrA",
    pos1 = "coordA",
    strand1 = "orientA",
    chr2 = "chrB",
    pos2 = "coordB",
    strand2 = "orientB",
    support_col = support_col,
    sample = sample,
    raw_id_col = "contig",
    extra_cols = c(
      breakpoint_type = "prelim.primary_brkpt_type",
      prelim_support_reads = "prelim.total",
      split_reads = "split",
      span_reads = "span",
      total_reads = "total"
    ),
    extra_numeric_cols = c(
      "prelim_support_reads", "split_reads", "span_reads", "total_reads"
    )
  )
}

virusfinder2_to_integrations <- function(df, path, sample = NULL) {
  ensure_import_columns(
    df,
    c(
      "Chromosome 1", "Position 1", "Strand 1",
      "Chromosome 2", "Position 2", "Strand 2"
    ),
    format = "VirusFinder2"
  )

  support_col <- "#Supporting reads (pair+softclipping)"
  if (!support_col %in% colnames(df)) {
    support_col <- NULL
  }

  standardize_paired_breakpoints(
    df = df,
    chr1 = "Chromosome 1",
    pos1 = "Position 1",
    strand1 = "Strand 1",
    chr2 = "Chromosome 2",
    pos2 = "Position 2",
    strand2 = "Strand 2",
    support_col = support_col,
    sample = sample,
    raw_id_col = NULL
  )
}

virus_clip_to_integrations <- function(df, path, sample = NULL) {
  ensure_import_columns(
    df,
    c(
      "Chromosome 1", "Position 1", "Strand 1",
      "Chromosome 2", "Position 2", "Strand 2"
    ),
    format = "Virus-Clip"
  )

  support_col <- "#Supporting reads (pair+softclipping)"
  if (!support_col %in% colnames(df)) {
    support_col <- NULL
  }

  standardize_paired_breakpoints(
    df = df,
    chr1 = "Chromosome 1",
    pos1 = "Position 1",
    strand1 = "Strand 1",
    chr2 = "Chromosome 2",
    pos2 = "Position 2",
    strand2 = "Strand 2",
    support_col = support_col,
    sample = sample,
    raw_id_col = NULL,
    preserve_cols = "confidence"
  )
}

bsvf_to_integrations <- function(df, path, sample = NULL) {
  ensure_import_columns(
    df,
    c(
      "Chr", "breakpoint", "virus-start", "virus-end",
      "virusstrand", "how-many-reads-support", "cluster-name"
    ),
    format = "BSVF"
  )

  support_col <- "how-many-reads-support"
  virus_start <- suppressWarnings(as.numeric(df[["virus-start"]]))
  virus_end <- suppressWarnings(as.numeric(df[["virus-end"]]))
  virus_pos <- round((virus_start + virus_end) / 2)

  out <- data.frame(
    host_chr = df[["Chr"]],
    host_pos = suppressWarnings(as.numeric(df[["breakpoint"]])),
    host_strand = "*",
    virus_chr = "virus",
    virus_pos = virus_pos,
    virus_strand = normalize_import_strand(df[["virusstrand"]]),
    support_reads = parse_support_value(df[[support_col]]),
    raw_record_id = as.character(df[["cluster-name"]]),
    virus_start = virus_start,
    virus_end = virus_end,
    cluster_name = as.character(df[["cluster-name"]]),
    stringsAsFactors = FALSE
  )

  if (!is.null(sample)) {
    out$sample <- validate_import_sample(sample)
  }

  as_integrations(out)
}

standardize_paired_breakpoints <- function(df,
                                           chr1,
                                           pos1,
                                           strand1,
                                           chr2,
                                           pos2,
                                           strand2,
                                           support_col,
                                           sample,
                                           raw_id_col = NULL,
                                           extra_cols = NULL,
                                           extra_numeric_cols = character(),
                                           preserve_cols = character()) {
  side1_host <- is_host_sequence(df[[chr1]])
  side2_host <- is_host_sequence(df[[chr2]])
  keep <- xor(side1_host, side2_host)

  if (!all(keep)) {
    warning(
      sum(!keep),
      " caller records were dropped because host and virus sides could not be distinguished.",
      call. = FALSE
    )
  }

  df <- df[keep, , drop = FALSE]
  side1_host <- side1_host[keep]
  if (nrow(df) == 0L) {
    stop("No valid host-virus integration records were found.", call. = FALSE)
  }

  out <- data.frame(
    host_chr = ifelse(side1_host, df[[chr1]], df[[chr2]]),
    host_pos = ifelse(side1_host, df[[pos1]], df[[pos2]]),
    host_strand = normalize_import_strand(ifelse(side1_host, df[[strand1]], df[[strand2]])),
    virus_chr = ifelse(side1_host, df[[chr2]], df[[chr1]]),
    virus_pos = ifelse(side1_host, df[[pos2]], df[[pos1]]),
    virus_strand = normalize_import_strand(ifelse(side1_host, df[[strand2]], df[[strand1]])),
    support_reads = if (is.null(support_col)) rep(NA_real_, nrow(df)) else parse_support_value(df[[support_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(sample)) {
    out$sample <- validate_import_sample(sample)
  }

  if (!is.null(raw_id_col)) {
    out$raw_record_id <- as.character(df[[raw_id_col]])
  }
  if (!is.null(extra_cols)) {
    for (target in names(extra_cols)) {
      source <- unname(extra_cols[[target]])
      out[[target]] <- if (target %in% extra_numeric_cols) {
        parse_support_value(df[[source]])
      } else {
        as.character(df[[source]])
      }
    }
  }
  if (length(preserve_cols) > 0L) {
    preserve_cols <- intersect(preserve_cols, colnames(df))
    for (col in preserve_cols) {
      out[[col]] <- df[[col]]
    }
  }

  as_integrations(out)
}

ensure_import_columns <- function(df, columns, format) {
  missing_cols <- setdiff(columns, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(
      format,
      " output must contain columns: ",
      paste(columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

is_host_sequence <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("(?i)^chr", "", x, perl = TRUE)
  grepl("^[0-9]+$|^X$|^Y$|^M$|^MT$", x, ignore.case = TRUE)
}

parse_support_value <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x <- trimws(as.character(x))
  plus_pair <- grepl("^\\d+(\\.\\d+)?\\s*\\+\\s*\\d+(\\.\\d+)?$", x)
  out <- suppressWarnings(as.numeric(x))

  if (any(plus_pair, na.rm = TRUE)) {
    out[plus_pair] <- vapply(
      strsplit(x[plus_pair], "\\+"),
      function(parts) sum(suppressWarnings(as.numeric(trimws(parts))), na.rm = TRUE),
      numeric(1)
    )
  }

  missing <- is.na(out)
  if (any(missing)) {
    first_number <- sub("^.*?([0-9]+(?:\\.[0-9]+)?).*$", "\\1", x[missing], perl = TRUE)
    out[missing] <- suppressWarnings(as.numeric(first_number))
  }

  out
}

normalize_import_strand <- function(x) {
  x <- trimws(as.character(x))
  x[!x %in% c("+", "-", "*")] <- "*"
  x
}

validate_import_sample <- function(sample) {
  sample <- as.character(sample)
  if (length(sample) != 1L || is.na(sample) || !nzchar(sample)) {
    stop("sample must be NULL or a single non-empty string.", call. = FALSE)
  }
  sample
}
