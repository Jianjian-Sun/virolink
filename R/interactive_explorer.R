#' Create an interactive integration explorer
#'
#' Build a lightweight HTML explorer that combines an interactive host-virus
#' breakpoint map with a searchable integration table. Points use stable
#' record identifiers as `ggiraph` data IDs, and the same identifiers are shown
#' in the table for cross-reference.
#'
#' @param integrations Integration records accepted by [as_integrations()].
#' @param host Host genome specification. Can be a built-in host name, a host
#'   genome table, or a [host_genome()] object.
#' @param virus Virus genome specification. Can be a [virus_genome()] object or
#'   a named numeric vector such as `c(HBV = 3215)`.
#' @param features Optional virus feature intervals from [virus_features()].
#' @param chrom_file Optional path to a chromosome-size table. Supply either
#'   `host` or `chrom_file`, not both.
#' @param group_by Optional column in `integrations` mapped to point color.
#' @param size_by Optional numeric column in `integrations` mapped to point size.
#' @param table_columns Optional columns to show in the data table. Core
#'   columns are used by default.
#' @param colors Optional vector of colors. A named vector can be used to set
#'   colors for specific groups.
#' @param point_size Point size used when `size_by` is `NULL`.
#' @param size_range Point-size range used when `size_by` is supplied.
#' @param width SVG width passed to [ggiraph::girafe()].
#' @param height SVG height passed to [ggiraph::girafe()].
#' @return An HTML tag list containing a `ggiraph` widget and a `DT` table.
#' @export
plot_interactive_explorer <- function(integrations, host = "hg38", virus,
                                      features = NULL, chrom_file = NULL,
                                      group_by = NULL, size_by = "support_reads",
                                      table_columns = NULL, colors = NULL,
                                      point_size = 2.8,
                                      size_range = c(2, 6),
                                      width = 10, height = 5) {
  require_interactive_packages()

  integrations <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)
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
  if (!is.numeric(width) || length(width) != 1L || is.na(width) || width <= 0 ||
      !is.numeric(height) || length(height) != 1L || is.na(height) || height <= 0) {
    stop("width and height must be single positive numbers.", call. = FALSE)
  }

  plot_data <- prepare_interactive_breakpoints(
    integrations = integrations,
    host = host,
    virus = virus,
    chrom_file = chrom_file
  )
  p <- build_interactive_breakpoint_plot(
    plot_data = plot_data,
    features = features,
    group_by = group_by,
    size_by = size_by,
    colors = colors,
    point_size = point_size,
    size_range = size_range
  )

  table_data <- prepare_interactive_table(plot_data$integrations, table_columns)
  widget <- ggiraph::girafe(
    ggobj = p,
    width_svg = width,
    height_svg = height,
    options = list(
      ggiraph::opts_hover(css = "stroke:black;stroke-width:1.5px;"),
      ggiraph::opts_selection(
        type = "single",
        css = "stroke:black;stroke-width:2px;",
        only_shiny = FALSE
      ),
      ggiraph::opts_toolbar(saveaspng = TRUE)
    )
  )
  table <- DT::datatable(
    table_data,
    rownames = FALSE,
    filter = "top",
    selection = "single",
    options = list(pageLength = min(10, nrow(table_data)), scrollX = TRUE)
  )

  htmltools::tagList(
    htmltools::tags$div(class = "virolink-interactive-map", widget),
    htmltools::tags$div(class = "virolink-interactive-table", table)
  )
}

require_interactive_packages <- function() {
  missing <- c(
    if (!requireNamespace("ggiraph", quietly = TRUE)) "ggiraph",
    if (!requireNamespace("DT", quietly = TRUE)) "DT",
    if (!requireNamespace("htmltools", quietly = TRUE)) "htmltools"
  )
  if (length(missing) > 0L) {
    stop(
      "Interactive exploration requires packages: ",
      paste(missing, collapse = ", "),
      ". Please install them first.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

prepare_interactive_breakpoints <- function(integrations, host, virus, chrom_file) {
  virus_obj <- if (inherits(virus, "vi_virus_genome")) {
    virus
  } else {
    virus_genome(virus)
  }
  host_df <- prepare_breakpoint_host_axis(host = host, chrom_file = chrom_file)

  integrations$.record_id <- sprintf("VI%04d", seq_len(nrow(integrations)))
  chr_index <- match(integrations$host_chr, host_df$chr)
  in_host <- !is.na(chr_index)
  in_host_range <- in_host
  in_host_range[in_host] <- integrations$host_pos[in_host] >= host_df$start[chr_index[in_host]] &
    integrations$host_pos[in_host] <= host_df$end[chr_index[in_host]]
  in_virus_range <- integrations$virus_pos >= 0 & integrations$virus_pos <= virus_obj$length
  keep <- in_host & in_host_range & in_virus_range
  if (!all(keep)) {
    warning(sum(!keep), " integration records were outside the host or virus genome and were omitted.",
            call. = FALSE)
  }

  integrations <- integrations[keep, , drop = FALSE]
  chr_index <- chr_index[keep]
  if (nrow(integrations) == 0L) {
    stop("No integration records matched the host and virus genomes.", call. = FALSE)
  }
  integrations$genome_pos <- host_df$offset[chr_index] +
    integrations$host_pos - host_df$start[chr_index]
  integrations$chromosome <- factor(integrations$host_chr, levels = host_df$chr, ordered = TRUE)
  integrations$.tooltip <- make_integration_tooltip(integrations, virus_obj$name)

  list(integrations = integrations, host = host_df, virus = virus_obj)
}

build_interactive_breakpoint_plot <- function(plot_data, features, group_by,
                                              size_by, colors, point_size,
                                              size_range) {
  integrations <- plot_data$integrations
  host_df <- plot_data$host
  virus_obj <- plot_data$virus

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = host_df,
      ggplot2::aes(xmin = .data$genome_start, xmax = .data$genome_end,
                   ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = host_df$band,
      color = NA
    )

  if (!is.null(features)) {
    feature_df <- prepare_breakpoint_features(features, virus_obj$length)
    if (nrow(feature_df) > 0L) {
      p <- p +
        ggplot2::geom_rect(
          data = feature_df,
          ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = .data$start, ymax = .data$end),
          inherit.aes = FALSE,
          fill = feature_df$fill,
          alpha = 0.18,
          color = NA
        )
    }
  }

  point_mapping <- ggplot2::aes(
    x = .data$genome_pos,
    y = .data$virus_pos,
    tooltip = .data$.tooltip,
    data_id = .data$.record_id
  )
  if (!is.null(group_by) && !is.null(size_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$genome_pos,
      y = .data$virus_pos,
      color = .data[[group_by]],
      size = .data[[size_by]],
      tooltip = .data$.tooltip,
      data_id = .data$.record_id
    )
  } else if (!is.null(group_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$genome_pos,
      y = .data$virus_pos,
      color = .data[[group_by]],
      tooltip = .data$.tooltip,
      data_id = .data$.record_id
    )
  } else if (!is.null(size_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$genome_pos,
      y = .data$virus_pos,
      size = .data[[size_by]],
      tooltip = .data$.tooltip,
      data_id = .data$.record_id
    )
  }

  if (is.null(group_by) && is.null(size_by)) {
    p <- p + ggiraph::geom_point_interactive(
      data = integrations,
      mapping = point_mapping,
      color = "#0072B2",
      size = point_size,
      alpha = 0.82,
      na.rm = TRUE
    )
  } else if (is.null(group_by)) {
    p <- p + ggiraph::geom_point_interactive(
      data = integrations,
      mapping = point_mapping,
      color = "#0072B2",
      alpha = 0.82,
      na.rm = TRUE
    )
  } else if (is.null(size_by)) {
    p <- p + ggiraph::geom_point_interactive(
      data = integrations,
      mapping = point_mapping,
      size = point_size,
      alpha = 0.82,
      na.rm = TRUE
    )
  } else {
    p <- p + ggiraph::geom_point_interactive(
      data = integrations,
      mapping = point_mapping,
      alpha = 0.82,
      na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = host_df$midpoint,
      labels = host_df$chr,
      expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, virus_obj$length),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(x = "Host chromosome", y = paste0(virus_obj$name, " position")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = if (is.null(group_by) && is.null(size_by)) "none" else "right"
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

make_integration_tooltip <- function(df, virus_name) {
  tooltip <- paste0(
    "ID: ", df$.record_id,
    "\nHost: chr", df$host_chr, ":", format(df$host_pos, scientific = FALSE, trim = TRUE),
    "\n", virus_name, ": ", format(df$virus_pos, scientific = FALSE, trim = TRUE),
    "\nSupporting reads: ", df$support_reads,
    "\nVirus strand: ", df$virus_strand
  )
  if ("sample" %in% colnames(df)) {
    tooltip <- paste0(tooltip, "\nSample: ", df$sample)
  }
  tooltip
}

prepare_interactive_table <- function(integrations, table_columns = NULL) {
  core_cols <- c(
    ".record_id", "sample", "host_chr", "host_pos",
    "host_strand", "virus_chr", "virus_pos", "virus_strand",
    "support_reads"
  )
  if (is.null(table_columns)) {
    table_columns <- intersect(core_cols, colnames(integrations))
  } else {
    unknown <- setdiff(table_columns, colnames(integrations))
    if (length(unknown) > 0L) {
      stop("Unknown table_columns: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    table_columns <- unique(c(".record_id", table_columns))
  }
  out <- integrations[, table_columns, drop = FALSE]
  names(out)[names(out) == ".record_id"] <- "record_id"
  out
}
