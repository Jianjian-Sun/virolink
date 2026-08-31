#' Create host feature annotations
#'
#' @param x A data frame with feature intervals, or `NULL` when columns are
#'   supplied separately.
#' @param chr Host chromosome names.
#' @param start Feature start positions.
#' @param end Feature end positions.
#' @param feature Feature labels.
#' @param type Optional feature type labels.
#' @param strand Optional strand labels.
#' @return A `vi_host_features` object.
#' @export
host_features <- function(x = NULL, chr = NULL, start = NULL, end = NULL,
                          feature = NULL, type = NULL, strand = NULL) {
  if (is.null(x)) {
    x <- data.frame(
      chr = chr,
      start = start,
      end = end,
      feature = feature,
      stringsAsFactors = FALSE
    )
    if (!is.null(type)) {
      x$type <- type
    }
    if (!is.null(strand)) {
      x$strand <- strand
    }
  }

  if (!is.data.frame(x)) {
    stop("x must be a data frame or NULL when chr/start/end/feature are supplied.",
         call. = FALSE)
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

  rename_one("chr", c("chr", "chrom", "chromosome", "seqnames"))
  rename_one("start", c("start", "start_pos", "from"))
  rename_one("end", c("end", "stop", "end_pos", "to"))
  rename_one("feature", c("feature", "gene", "name", "label"))
  rename_one("type", c("type", "feature_type"))
  rename_one("strand", c("strand"))

  required_cols <- c("chr", "start", "end", "feature")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Host features must contain columns: chr, start, end, feature",
         call. = FALSE)
  }
  if (!"type" %in% colnames(df)) {
    df$type <- "feature"
  }
  if (!"strand" %in% colnames(df)) {
    df$strand <- "."
  }

  out <- df[, c("chr", "start", "end", "feature", "type", "strand"), drop = FALSE]
  out$chr <- trimws(gsub("(?i)^chr", "", as.character(out$chr)))
  out$start <- suppressWarnings(as.numeric(out$start))
  out$end <- suppressWarnings(as.numeric(out$end))
  out$feature <- as.character(out$feature)
  out$type <- as.character(out$type)
  out$strand <- as.character(out$strand)

  keep <- !is.na(out$chr) & nzchar(out$chr) &
    !is.na(out$start) & !is.na(out$end) & out$end > out$start &
    !is.na(out$feature) & nzchar(out$feature)
  if (!all(keep)) {
    warning(sum(!keep), " host features were dropped because they were incomplete or invalid.",
            call. = FALSE)
  }
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0L) {
    stop("No valid host features were found.", call. = FALSE)
  }

  out <- out[order(out$chr, out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("vi_host_features", "data.frame")
  out
}

#' Plot local context around a host breakpoint
#'
#' Display integration records near a selected host breakpoint and optionally
#' overlay local host feature intervals such as genes, transcripts, regulatory
#' elements, or user-defined evidence tracks.
#'
#' @param integrations Integration records accepted by [as_integrations()].
#' @param chr Host chromosome containing the focal breakpoint.
#' @param pos Host position of the focal breakpoint.
#' @param window Number of base pairs shown on each side of `pos`.
#' @param annotations Optional host feature intervals from [host_features()].
#' @param group_by Optional column in `integrations` mapped to point color.
#' @param size_by Optional numeric column in `integrations` mapped to point size.
#' @param colors Optional vector of colors. A named vector can be used to set
#'   colors for specific groups.
#' @param annotation_fill Optional feature colors by annotation type.
#' @param point_size Point size used when `size_by` is `NULL`.
#' @param size_range Point-size range used when `size_by` is supplied.
#' @param label_annotations Logical. Whether to label annotation intervals.
#' @param show_feature_legend Logical. Whether to show the annotation type
#'   legend when host features are supplied.
#' @return A `ggplot` object.
#' @export
plot_locus_context <- function(integrations, chr, pos, window = 5000,
                               annotations = NULL, group_by = NULL,
                               size_by = "support_reads", colors = NULL,
                               annotation_fill = NULL, point_size = 2.5,
                               size_range = c(1.5, 5),
                               label_annotations = TRUE,
                               show_feature_legend = TRUE) {
  integrations <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)
  chr <- normalize_locus_chr(chr)
  pos <- normalize_locus_pos(pos)
  if (!is.numeric(window) || length(window) != 1L || is.na(window) || window <= 0) {
    stop("window must be a single positive number.", call. = FALSE)
  }
  validate_single_column_arg(group_by, integrations, "group_by", allow_null = TRUE)
  validate_single_column_arg(size_by, integrations, "size_by", allow_null = TRUE)
  if (!is.null(size_by) && !is.numeric(integrations[[size_by]])) {
    stop("size_by must refer to a numeric column.", call. = FALSE)
  }
  if (!is.numeric(point_size) || length(point_size) != 1L || is.na(point_size) ||
      point_size <= 0) {
    stop("point_size must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(size_range) || length(size_range) != 2L ||
      any(is.na(size_range)) || any(size_range <= 0) || size_range[1] > size_range[2]) {
    stop("size_range must contain two positive numbers in increasing order.", call. = FALSE)
  }
  if (!is.logical(label_annotations) || length(label_annotations) != 1L ||
      is.na(label_annotations)) {
    stop("label_annotations must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_feature_legend) || length(show_feature_legend) != 1L ||
      is.na(show_feature_legend)) {
    stop("show_feature_legend must be TRUE or FALSE.", call. = FALSE)
  }

  xmin <- max(0, pos - window)
  xmax <- pos + window
  locus_df <- integrations[
    integrations$host_chr == chr &
      integrations$host_pos >= xmin &
      integrations$host_pos <= xmax,
    ,
    drop = FALSE
  ]
  if (nrow(locus_df) == 0L) {
    stop("No integration records were found in the requested locus window.",
         call. = FALSE)
  }
  locus_df$track_y <- 1

  feature_df <- NULL
  if (!is.null(annotations)) {
    feature_df <- prepare_locus_features(
      annotations = annotations,
      chr = chr,
      xmin = xmin,
      xmax = xmax,
      annotation_fill = annotation_fill
    )
  }

  point_mapping <- ggplot2::aes(x = .data$host_pos, y = .data$track_y)
  if (!is.null(group_by) && !is.null(size_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$host_pos,
      y = .data$track_y,
      color = .data[[group_by]],
      size = .data[[size_by]]
    )
  } else if (!is.null(group_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$host_pos,
      y = .data$track_y,
      color = .data[[group_by]]
    )
  } else if (!is.null(size_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$host_pos,
      y = .data$track_y,
      size = .data[[size_by]]
    )
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = pos, linetype = "dashed", color = "grey45") +
    ggplot2::geom_hline(yintercept = 1, color = "grey85", linewidth = 0.3)

  if (!is.null(feature_df) && nrow(feature_df) > 0L) {
    feature_fill <- attr(feature_df, "fill_map")
    p <- p +
      ggplot2::geom_rect(
        data = feature_df,
        ggplot2::aes(xmin = .data$plot_start, xmax = .data$plot_end,
                     ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$type),
        color = "white",
        linewidth = 0.25,
        alpha = 0.45
      ) +
      ggplot2::scale_fill_manual(
        values = feature_fill,
        name = "Host feature type",
        guide = if (show_feature_legend) "legend" else "none"
      )
    if (isTRUE(label_annotations)) {
      p <- p +
        ggplot2::geom_text(
          data = feature_df,
          ggplot2::aes(x = (.data$plot_start + .data$plot_end) / 2,
                       y = .data$ymax + 0.025,
                       label = .data$feature),
          size = 2.8,
          color = "grey25",
          vjust = 0,
          check_overlap = TRUE
        )
    }
  }

  if (is.null(group_by) && is.null(size_by)) {
    p <- p + ggplot2::geom_point(
      data = locus_df,
      mapping = point_mapping,
      color = "#0072B2",
      size = point_size,
      alpha = 0.8,
      na.rm = TRUE
    )
  } else if (is.null(group_by)) {
    p <- p + ggplot2::geom_point(
      data = locus_df,
      mapping = point_mapping,
      color = "#0072B2",
      alpha = 0.8,
      na.rm = TRUE
    )
  } else if (is.null(size_by)) {
    p <- p + ggplot2::geom_point(
      data = locus_df,
      mapping = point_mapping,
      size = point_size,
      alpha = 0.8,
      na.rm = TRUE
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = locus_df,
      mapping = point_mapping,
      alpha = 0.8,
      na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      limits = c(xmin, xmax),
      labels = pretty_bp_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(1, 0.45),
      labels = c("Integrations", "Host features"),
      limits = c(0, 1.2),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = paste0("chr", chr, " position"), y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = if (is.null(group_by) && is.null(size_by) &&
                            is.null(feature_df)) "none" else "right"
    )

  if (!is.null(group_by)) {
    p <- p + ggplot2::labs(color = group_by)
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_color_manual(values = colors)
    }
  }
  if (!is.null(size_by)) {
    p <- p +
      ggplot2::scale_size_continuous(range = size_range) +
      ggplot2::labs(size = size_by)
  }

  p
}

normalize_locus_chr <- function(chr) {
  chr <- trimws(gsub("(?i)^chr", "", as.character(chr)))
  if (length(chr) != 1L || is.na(chr) || !nzchar(chr)) {
    stop("chr must be a single non-empty chromosome name.", call. = FALSE)
  }
  chr
}

normalize_locus_pos <- function(pos) {
  pos <- suppressWarnings(as.numeric(pos))
  if (length(pos) != 1L || is.na(pos) || pos < 0) {
    stop("pos must be a single non-negative number.", call. = FALSE)
  }
  pos
}

prepare_locus_features <- function(annotations, chr, xmin, xmax,
                                   annotation_fill = NULL) {
  feature_df <- as.data.frame(host_features(annotations), stringsAsFactors = FALSE)
  feature_df <- feature_df[
    feature_df$chr == chr &
      feature_df$end >= xmin &
      feature_df$start <= xmax,
    ,
    drop = FALSE
  ]
  if (nrow(feature_df) == 0L) {
    return(feature_df)
  }

  feature_df$plot_start <- pmax(feature_df$start, xmin)
  feature_df$plot_end <- pmin(feature_df$end, xmax)
  feature_df <- layout_locus_features(feature_df)
  fill_map <- if (is.null(annotation_fill)) {
    default_locus_feature_colors(unique(feature_df$type))
  } else {
    normalize_named_colors(annotation_fill, unique(feature_df$type))
  }
  attr(feature_df, "fill_map") <- fill_map
  feature_df
}

default_locus_feature_colors <- function(types) {
  known_colors <- c(
    "gene" = "#DDEAD7",
    "transcript" = "#D8E8F2",
    "regulatory" = "#F2E2D8",
    "enhancer" = "#F2E2D8",
    "promoter" = "#E8DDF2",
    "repeat" = "#F2EFD8",
    "feature" = "#E8E8E8"
  )
  fill_map <- known_colors[types]
  missing <- is.na(fill_map)
  if (any(missing)) {
    fill_map[missing] <- grDevices::hcl.colors(sum(missing), "Pastel 1")
  }
  stats::setNames(unname(fill_map), types)
}

layout_locus_features <- function(feature_df) {
  feature_df <- feature_df[order(feature_df$plot_start, feature_df$plot_end), , drop = FALSE]
  feature_df$level <- 1L
  level_ends <- numeric()
  for (i in seq_len(nrow(feature_df))) {
    placed <- FALSE
    for (level in seq_along(level_ends)) {
      if (feature_df$plot_start[i] >= level_ends[level]) {
        feature_df$level[i] <- level
        level_ends[level] <- feature_df$plot_end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      level_ends <- c(level_ends, feature_df$plot_end[i])
      feature_df$level[i] <- length(level_ends)
    }
  }
  feature_df$ymin <- 0.18 + (feature_df$level - 1) * 0.12
  feature_df$ymax <- feature_df$ymin + 0.08
  feature_df
}

pretty_bp_labels <- function(x) {
  ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 3), " Mb"), paste0(round(x / 1e3, 1), " kb"))
}
