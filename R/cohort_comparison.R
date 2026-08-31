#' Plot cohort-level integration summaries
#'
#' Compare integration burden or supporting-read burden across samples, cohorts,
#' or any other user-supplied column in the integration table. When `group_by`
#' is `NULL`, all events are combined into one overall group. A second column
#' can be used to show stacked composition within each group.
#'
#' @param integrations Integration records accepted by [as_integrations()].
#' @param group_by Column used for the main comparison groups.
#' @param fill_by Optional column used for stacked composition within each
#'   group.
#' @param measure Quantity to summarize: number of integration events
#'   (`"events"`) or total supporting reads (`"support_reads"`).
#' @param normalize Logical. If `TRUE`, show within-group proportions instead
#'   of raw totals.
#' @param top_n Optional number of `fill_by` categories to keep. Remaining
#'   categories are collapsed into `"Other"`.
#' @param colors Optional fill colors. A named vector can be used to set colors
#'   for specific `fill_by` values.
#' @return A `ggplot` object.
#' @export
plot_cohort_comparison <- function(integrations, group_by = NULL,
                                   fill_by = NULL,
                                   measure = c("events", "support_reads"),
                                   normalize = FALSE,
                                   top_n = NULL,
                                   colors = NULL) {
  measure <- match.arg(measure)
  integrations <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)

  validate_single_column_arg(group_by, integrations, "group_by", allow_null = TRUE)
  validate_single_column_arg(fill_by, integrations, "fill_by", allow_null = TRUE)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("normalize must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(top_n)) {
    if (is.null(fill_by)) {
      stop("top_n can only be used when fill_by is supplied.", call. = FALSE)
    }
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
        top_n < 1 || top_n != floor(top_n)) {
      stop("top_n must be a single positive integer.", call. = FALSE)
    }
  }

  plot_data <- prepare_cohort_summary(
    integrations = integrations,
    group_by = group_by,
    fill_by = fill_by,
    measure = measure,
    normalize = normalize,
    top_n = top_n
  )

  y_label <- if (normalize) {
    "Proportion of integration burden"
  } else if (measure == "events") {
    "Integration events"
  } else {
    "Supporting reads"
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$group, y = .data$value))

  if (is.null(fill_by)) {
    p <- p +
      ggplot2::geom_col(fill = "#0072B2", width = 0.75) +
      ggplot2::theme(legend.position = "none")
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(fill = .data$fill), width = 0.75) +
      ggplot2::labs(fill = fill_by)
    fill_levels <- levels(plot_data$fill)
    fill_colors <- if (!is.null(colors)) {
      normalize_named_colors(colors, fill_levels)
    } else {
      normalize_named_colors(grDevices::hcl.colors(length(fill_levels), "Dark 3"), fill_levels)
    }
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  }

  if (normalize) {
    p <- p + ggplot2::scale_y_continuous(labels = percent_labels)
  }

  p +
    ggplot2::labs(
      x = if (is.null(group_by)) "All integrations" else group_by,
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}

prepare_cohort_summary <- function(integrations, group_by, fill_by, measure,
                                   normalize, top_n) {
  df <- integrations
  if (is.null(group_by)) {
    df$.group <- "All integrations"
    group_levels <- "All integrations"
  } else {
    df$.group <- as.character(df[[group_by]])
    df$.group[is.na(df$.group) | !nzchar(df$.group)] <- "Missing"
    group_levels <- unique(df$.group)
  }

  if (is.null(fill_by)) {
    df$.fill <- "integration"
    fill_levels <- "integration"
  } else {
    df$.fill <- as.character(df[[fill_by]])
    df$.fill[is.na(df$.fill) | !nzchar(df$.fill)] <- "Missing"
    fill_levels <- unique(df$.fill)
  }

  df$.value <- if (measure == "events") {
    1
  } else {
    suppressWarnings(as.numeric(df$support_reads))
  }
  df <- df[!is.na(df$.value), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No valid values were available for cohort comparison.", call. = FALSE)
  }

  if (!is.null(top_n) && length(fill_levels) > top_n) {
    fill_totals <- stats::aggregate(
      df$.value,
      by = list(fill = df$.fill),
      FUN = sum
    )
    names(fill_totals)[2] <- "total"
    fill_totals <- fill_totals[order(fill_totals$total, decreasing = TRUE), , drop = FALSE]
    keep_fill <- utils::head(fill_totals$fill, top_n)
    df$.fill[!df$.fill %in% keep_fill] <- "Other"
    fill_levels <- c(keep_fill, "Other")
  }

  out <- stats::aggregate(
    df$.value,
    by = list(group = df$.group, fill = df$.fill),
    FUN = sum
  )
  names(out)[3] <- "value"

  if (normalize) {
    group_totals <- stats::aggregate(out$value, by = list(group = out$group), FUN = sum)
    names(group_totals)[2] <- "total"
    out <- merge(out, group_totals, by = "group", all.x = TRUE, sort = FALSE)
    out$value <- out$value / out$total
    out$total <- NULL
  }

  out$group <- factor(out$group, levels = group_levels, ordered = TRUE)
  out$fill <- factor(out$fill, levels = fill_levels, ordered = TRUE)
  out <- out[order(out$group, out$fill), , drop = FALSE]
  rownames(out) <- NULL
  out
}

percent_labels <- function(x) {
  paste0(round(x * 100), "%")
}
