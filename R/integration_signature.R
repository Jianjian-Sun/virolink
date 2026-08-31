#' Plot integration signatures for clonality or breakpoint orientation
#'
#' Summarize integration records as evidence views for candidate clone-like
#' breakpoint clusters or viral-strand orientation. The clonality view groups
#' nearby breakpoints using a simple distance window; it is intended as a
#' visualization of recurrent evidence, not as an automatic biological
#' clonality call.
#'
#' @param integrations Integration records accepted by [as_integrations()].
#' @param mode Signature view to draw. `"clonality"` summarizes candidate
#'   recurrent breakpoint clusters; `"strand"` summarizes `virus_strand`
#'   orientation.
#' @param group_by Optional column used for stacked composition.
#' @param cluster_by Clustering rule for `mode = "clonality"`. `"host"` groups
#'   nearby host breakpoints; `"host_virus"` additionally requires nearby viral
#'   breakpoints.
#' @param cluster_window Distance in base pairs used to group nearby
#'   breakpoints.
#' @param measure Quantity to summarize: number of integration events
#'   (`"events"`) or total supporting reads (`"support_reads"`).
#' @param top_n Number of candidate clusters shown in `mode = "clonality"`.
#' @param colors Optional fill colors. A named vector can be used to set colors
#'   for specific groups.
#' @return A `ggplot` object.
#' @export
plot_integration_signature <- function(integrations,
                                       mode = c("clonality", "strand"),
                                       group_by = NULL,
                                       cluster_by = c("host", "host_virus"),
                                       cluster_window = 1000,
                                       measure = c("support_reads", "events"),
                                       top_n = 10,
                                       colors = NULL) {
  mode <- match.arg(mode)
  cluster_by <- match.arg(cluster_by)
  measure <- match.arg(measure)
  integrations <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)

  validate_single_column_arg(group_by, integrations, "group_by", allow_null = TRUE)
  if (!is.numeric(cluster_window) || length(cluster_window) != 1L ||
      is.na(cluster_window) || cluster_window <= 0) {
    stop("cluster_window must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
      top_n < 1 || top_n != floor(top_n)) {
    stop("top_n must be a single positive integer.", call. = FALSE)
  }

  if (mode == "clonality") {
    plot_data <- prepare_clonality_signature(
      integrations = integrations,
      group_by = group_by,
      cluster_by = cluster_by,
      cluster_window = cluster_window,
      measure = measure,
      top_n = top_n
    )
    y_label <- if (measure == "support_reads") "Supporting reads" else "Integration events"
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$cluster, y = .data$value))
    if (is.null(group_by)) {
      p <- p +
        ggplot2::geom_col(fill = "#0072B2", width = 0.75) +
        ggplot2::theme(legend.position = "none")
    } else {
      p <- p +
        ggplot2::geom_col(ggplot2::aes(fill = .data$group), width = 0.75) +
        ggplot2::labs(fill = group_by)
      p <- p + ggplot2::scale_fill_manual(
        values = signature_fill_colors(plot_data$group, colors)
      )
    }

    return(p +
      ggplot2::labs(x = "Candidate breakpoint cluster", y = y_label) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
      ))
  }

  plot_data <- prepare_strand_signature(
    integrations = integrations,
    group_by = group_by,
    measure = measure
  )
  y_label <- if (measure == "support_reads") "Supporting reads" else "Integration events"
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$strand, y = .data$value))
  if (is.null(group_by)) {
    p <- p +
      ggplot2::geom_col(fill = "#0072B2", width = 0.75) +
      ggplot2::theme(legend.position = "none")
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(fill = .data$group), position = "dodge", width = 0.75) +
      ggplot2::labs(fill = group_by) +
      ggplot2::scale_fill_manual(values = signature_fill_colors(plot_data$group, colors))
  }

  p +
    ggplot2::labs(x = "Virus strand", y = y_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

prepare_clonality_signature <- function(integrations, group_by, cluster_by,
                                        cluster_window, measure, top_n) {
  clustered <- assign_breakpoint_clusters(
    integrations = integrations,
    cluster_by = cluster_by,
    cluster_window = cluster_window
  )
  clustered$.group <- signature_group_values(clustered, group_by)
  clustered$.value <- signature_measure_values(clustered, measure)
  clustered <- clustered[!is.na(clustered$.value), , drop = FALSE]
  if (nrow(clustered) == 0L) {
    stop("No valid values were available for integration signatures.", call. = FALSE)
  }

  cluster_totals <- stats::aggregate(
    clustered$.value,
    by = list(cluster_id = clustered$.cluster_id),
    FUN = sum
  )
  names(cluster_totals)[2] <- "total"
  cluster_totals <- cluster_totals[order(cluster_totals$total, decreasing = TRUE), , drop = FALSE]
  keep_clusters <- utils::head(cluster_totals$cluster_id, top_n)
  clustered <- clustered[clustered$.cluster_id %in% keep_clusters, , drop = FALSE]

  out <- stats::aggregate(
    clustered$.value,
    by = list(cluster_id = clustered$.cluster_id, group = clustered$.group),
    FUN = sum
  )
  names(out)[3] <- "value"

  cluster_labels <- unique(clustered[, c(".cluster_id", ".cluster_label"), drop = FALSE])
  out <- merge(out, cluster_labels, by.x = "cluster_id", by.y = ".cluster_id",
               all.x = TRUE, sort = FALSE)
  cluster_order <- cluster_totals$cluster_id[cluster_totals$cluster_id %in% keep_clusters]
  label_order <- cluster_labels$.cluster_label[match(cluster_order, cluster_labels$.cluster_id)]
  out$cluster <- factor(out$.cluster_label, levels = label_order, ordered = TRUE)
  out$group <- factor(out$group, levels = unique(clustered$.group), ordered = TRUE)
  out <- out[order(out$cluster, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}

prepare_strand_signature <- function(integrations, group_by, measure) {
  df <- integrations
  df$.group <- signature_group_values(df, group_by)
  df$.strand <- as.character(df$virus_strand)
  df$.strand[is.na(df$.strand) | !nzchar(df$.strand)] <- "Missing"
  df$.value <- signature_measure_values(df, measure)
  df <- df[!is.na(df$.value), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No valid values were available for integration signatures.", call. = FALSE)
  }

  out <- stats::aggregate(
    df$.value,
    by = list(strand = df$.strand, group = df$.group),
    FUN = sum
  )
  names(out)[3] <- "value"
  strand_levels <- c("+", "-", setdiff(unique(df$.strand), c("+", "-")))
  out$strand <- factor(out$strand, levels = strand_levels, ordered = TRUE)
  out$group <- factor(out$group, levels = unique(df$.group), ordered = TRUE)
  out <- out[order(out$strand, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}

assign_breakpoint_clusters <- function(integrations, cluster_by, cluster_window) {
  df <- integrations[order(integrations$host_chr, integrations$host_pos, integrations$virus_pos), ,
                     drop = FALSE]
  df$.cluster_id <- integer(nrow(df))
  df$.cluster_label <- character(nrow(df))

  cluster_id <- 0L
  cluster_chr <- NULL
  cluster_host_start <- NA_real_
  cluster_host_end <- NA_real_
  cluster_virus_min <- NA_real_
  cluster_virus_max <- NA_real_

  for (i in seq_len(nrow(df))) {
    same_chr <- !is.null(cluster_chr) && identical(df$host_chr[i], cluster_chr)
    host_close <- same_chr && (df$host_pos[i] - cluster_host_start) <= cluster_window
    virus_close <- TRUE
    if (cluster_by == "host_virus" && same_chr) {
      next_virus_min <- min(cluster_virus_min, df$virus_pos[i], na.rm = TRUE)
      next_virus_max <- max(cluster_virus_max, df$virus_pos[i], na.rm = TRUE)
      virus_close <- (next_virus_max - next_virus_min) <= cluster_window
    }

    if (!same_chr || !host_close || !virus_close) {
      cluster_id <- cluster_id + 1L
      cluster_chr <- df$host_chr[i]
      cluster_host_start <- df$host_pos[i]
      cluster_host_end <- df$host_pos[i]
      cluster_virus_min <- df$virus_pos[i]
      cluster_virus_max <- df$virus_pos[i]
    } else {
      cluster_host_end <- max(cluster_host_end, df$host_pos[i], na.rm = TRUE)
      cluster_virus_min <- min(cluster_virus_min, df$virus_pos[i], na.rm = TRUE)
      cluster_virus_max <- max(cluster_virus_max, df$virus_pos[i], na.rm = TRUE)
    }

    df$.cluster_id[i] <- cluster_id
    df$.cluster_label[i] <- format_cluster_label(
      chr = cluster_chr,
      start = cluster_host_start,
      end = cluster_host_end
    )
  }

  cluster_start <- stats::aggregate(
    df$host_pos,
    by = list(cluster_id = df$.cluster_id, chr = df$host_chr),
    FUN = min
  )
  names(cluster_start)[3] <- "start"
  cluster_end <- stats::aggregate(
    df$host_pos,
    by = list(cluster_id = df$.cluster_id),
    FUN = max
  )
  names(cluster_end)[2] <- "end"
  cluster_meta <- merge(cluster_start, cluster_end, by = "cluster_id", sort = FALSE)
  cluster_meta$label <- mapply(
    format_cluster_label,
    chr = cluster_meta$chr,
    start = cluster_meta$start,
    end = cluster_meta$end,
    USE.NAMES = FALSE
  )
  df$.cluster_label <- cluster_meta$label[match(df$.cluster_id, cluster_meta$cluster_id)]
  df
}

format_cluster_label <- function(chr, start, end) {
  mid <- (start + end) / 2
  if (abs(mid) >= 1e6) {
    paste0("chr", chr, ":", round(mid / 1e6, 3), " Mb")
  } else {
    paste0("chr", chr, ":", round(mid / 1e3, 1), " kb")
  }
}

signature_group_values <- function(df, group_by) {
  if (is.null(group_by)) {
    return(rep("integration", nrow(df)))
  }
  group <- as.character(df[[group_by]])
  group[is.na(group) | !nzchar(group)] <- "Missing"
  group
}

signature_measure_values <- function(df, measure) {
  if (measure == "events") {
    return(rep(1, nrow(df)))
  }
  suppressWarnings(as.numeric(df$support_reads))
}

signature_fill_colors <- function(groups, colors = NULL) {
  groups <- levels(factor(groups))
  if (!is.null(colors)) {
    return(normalize_named_colors(colors, groups))
  }
  normalize_named_colors(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
}
