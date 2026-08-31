#' Plot the host-genome integration landscape
#'
#' Display viral integration events across a host genome as individual sites
#' or binned event density. Chromosomes are placed on a continuous genomic axis
#' in the order supplied by the host genome.
#'
#' @param integrations Integration records accepted by [as_integrations()].
#' @param host Host genome specification. Can be a built-in host name, a host
#'   genome table, or a [host_genome()] object.
#' @param chrom_file Optional path to a chromosome-size table. Supply either
#'   `host` or `chrom_file`, not both.
#' @param mode Display individual integration sites (`"sites"`) or counts in
#'   genomic bins (`"density"`).
#' @param group_by Optional column in `integrations` mapped to color.
#' @param bin_size Bin width in base pairs for `mode = "density"`.
#' @param colors Optional vector of colors. A named vector can be used to set
#'   colors for specific groups.
#' @param point_size Point size for `mode = "sites"`.
#' @return A `ggplot` object.
#' @export
plot_integration_landscape <- function(integrations, host = "hg38",
                                       chrom_file = NULL,
                                       mode = c("sites", "density"),
                                       group_by = NULL,
                                       bin_size = 5e6,
                                       colors = NULL,
                                       point_size = 2.2) {
  mode <- match.arg(mode)
  integrations <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)

  if (!is.null(group_by)) {
    if (!is.character(group_by) || length(group_by) != 1L || is.na(group_by) ||
        !nzchar(group_by)) {
      stop("group_by must be NULL or a single column name.", call. = FALSE)
    }
    if (!group_by %in% colnames(integrations)) {
      stop("Unknown group_by column: ", group_by, call. = FALSE)
    }
  }
  if (!is.numeric(bin_size) || length(bin_size) != 1L || is.na(bin_size) ||
      bin_size <= 0) {
    stop("bin_size must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(point_size) || length(point_size) != 1L || is.na(point_size) ||
      point_size <= 0) {
    stop("point_size must be a single positive number.", call. = FALSE)
  }

  host_df <- resolve_host_chrom_sizes(host = host, chrom_file = chrom_file)
  host_df <- host_df[!duplicated(host_df$chr), , drop = FALSE]
  host_df$length <- host_df$end - host_df$start
  if (any(host_df$length <= 0)) {
    stop("Host chromosome end positions must be greater than start positions.", call. = FALSE)
  }
  host_df$offset <- c(0, utils::head(cumsum(host_df$length), -1L))
  host_df$genome_start <- host_df$offset
  host_df$genome_end <- host_df$offset + host_df$length
  host_df$midpoint <- (host_df$genome_start + host_df$genome_end) / 2
  host_df$band <- rep(c("#F2F2F2", "#FFFFFF"), length.out = nrow(host_df))

  chr_index <- match(integrations$host_chr, host_df$chr)
  in_host <- !is.na(chr_index)
  in_range <- in_host
  in_range[in_host] <- integrations$host_pos[in_host] >= host_df$start[chr_index[in_host]] &
    integrations$host_pos[in_host] <= host_df$end[chr_index[in_host]]
  keep <- in_host & in_range
  if (!all(keep)) {
    warning(sum(!keep), " integration records were outside the host genome and were omitted.",
            call. = FALSE)
  }
  integrations <- integrations[keep, , drop = FALSE]
  chr_index <- chr_index[keep]
  if (nrow(integrations) == 0L) {
    stop("No integration records matched the host genome.", call. = FALSE)
  }

  integrations$genome_pos <- host_df$offset[chr_index] +
    integrations$host_pos - host_df$start[chr_index]
  integrations$chromosome <- factor(
    integrations$host_chr,
    levels = host_df$chr,
    ordered = TRUE
  )

  plot_data <- integrations
  y_label <- "Supporting reads"
  if (mode == "density") {
    plot_data$bin <- floor((plot_data$host_pos - host_df$start[chr_index]) / bin_size)
    split_cols <- c("chromosome", "bin")
    if (!is.null(group_by)) {
      split_cols <- c(split_cols, group_by)
    }
    counts <- stats::aggregate(
      rep(1L, nrow(plot_data)),
      by = plot_data[split_cols],
      FUN = sum
    )
    names(counts)[ncol(counts)] <- "count"
    count_chr_index <- match(as.character(counts$chromosome), host_df$chr)
    counts$genome_start <- host_df$offset[count_chr_index] + counts$bin * bin_size
    counts$genome_end <- host_df$offset[count_chr_index] + pmin(
      (counts$bin + 1) * bin_size,
      host_df$length[count_chr_index]
    )
    if (is.null(group_by)) {
      counts$ymin <- 0
      counts$ymax <- counts$count
    } else {
      bin_id <- interaction(counts$chromosome, counts$bin, drop = TRUE)
      counts$ymax <- stats::ave(counts$count, bin_id, FUN = cumsum)
      counts$ymin <- counts$ymax - counts$count
    }
    plot_data <- counts
    y_label <- "Integration events per bin"
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = host_df,
      ggplot2::aes(xmin = .data$genome_start, xmax = .data$genome_end,
                   ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = host_df$band,
      color = NA
    )

  if (mode == "sites") {
    if (is.null(group_by)) {
      p <- p + ggplot2::geom_point(
        data = plot_data,
        mapping = ggplot2::aes(x = .data$genome_pos, y = .data$support_reads),
        size = point_size, alpha = 0.8, color = "#0072B2", na.rm = TRUE
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = plot_data,
        mapping = ggplot2::aes(
          x = .data$genome_pos,
          y = .data$support_reads,
          color = .data[[group_by]]
        ),
        size = point_size, alpha = 0.8, na.rm = TRUE
      )
    }
  } else {
    if (is.null(group_by)) {
      p <- p + ggplot2::geom_rect(
        data = plot_data,
        mapping = ggplot2::aes(
          xmin = .data$genome_start,
          xmax = .data$genome_end,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        alpha = 0.78, fill = "#0072B2", color = NA
      )
    } else {
      p <- p + ggplot2::geom_rect(
        data = plot_data,
        mapping = ggplot2::aes(
          xmin = .data$genome_start,
          xmax = .data$genome_end,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data[[group_by]]
        ),
        alpha = 0.78, color = NA
      )
    }
  }

  if (!is.null(colors)) {
    if (mode == "sites") {
      p <- p + ggplot2::scale_color_manual(values = colors)
    } else {
      p <- p + ggplot2::scale_fill_manual(values = colors)
    }
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = host_df$midpoint,
      labels = host_df$chr,
      expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) +
    ggplot2::labs(x = "Host chromosome", y = y_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = if (is.null(group_by)) "none" else "right"
    )

  if (!is.null(group_by)) {
    if (mode == "sites") {
      p <- p + ggplot2::labs(color = group_by)
    } else {
      p <- p + ggplot2::labs(fill = group_by)
    }
  }

  p
}
