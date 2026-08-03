#' Read integration records
#'
#' @param path Path to a tabular file containing integration records.
#' @param sep Field separator. Defaults to tab.
#' @return A \code{vi_integrations} object.
#' @export
read_integrations <- function(path, sep = "\t") {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  df <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  as_integrations(df)
}

#' Coerce to integration records
#'
#' @param x A data frame, file path, or \code{vi_integrations} object.
#' @return A \code{vi_integrations} object.
#' @export
as_integrations <- function(x) {
  if (inherits(x, "vi_integrations")) {
    return(validate_integrations(x))
  }

  if (is.character(x) && length(x) == 1L) {
    return(read_integrations(x))
  }

  if (!is.data.frame(x)) {
    stop("x must be a data frame, file path, or vi_integrations object.", call. = FALSE)
  }

  validate_integrations(x)
}

#' Validate integration records
#'
#' @param x Integration records.
#' @return A standardized \code{vi_integrations} object.
#' @export
validate_integrations <- function(x) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)

  rename_one <- function(target, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0) {
      return(NULL)
    }
    df[[target]] <<- df[[hit[1]]]
    invisible(NULL)
  }

  rename_one("host_chr", c("host_chr", "chr"))
  rename_one("host_pos", c("host_pos", "host_loc"))
  rename_one("virus_pos", c("virus_pos", "viral_loc"))
  rename_one("sample", c("sample", "Tumor_Sample_Barcode"))
  rename_one("support", c("support", "reads"))
  rename_one("virus_strand", c("virus_strand", "viral_strand"))
  rename_one("method", c("method"))

  required_cols <- c("host_chr", "host_pos", "virus_pos", "sample", "support", "virus_strand", "method")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Integration data must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df$host_chr <- trimws(gsub("(?i)^chr", "", as.character(df$host_chr)))
  df$host_pos <- suppressWarnings(as.numeric(df$host_pos))
  df$virus_pos <- suppressWarnings(as.numeric(df$virus_pos))
  df$sample <- as.character(df$sample)
  df$support <- suppressWarnings(as.numeric(df$support))
  df$virus_strand <- as.character(df$virus_strand)
  df$method <- as.character(df$method)

  keep <- !is.na(df$host_chr) & nzchar(df$host_chr) &
    !is.na(df$host_pos) & !is.na(df$virus_pos) &
    !is.na(df$sample) & nzchar(df$sample)

  if (!all(keep)) {
    warning(sum(!keep), " rows were dropped because they were incomplete or invalid.")
  }

  df <- df[keep, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No valid integration records were found.", call. = FALSE)
  }

  rownames(df) <- NULL
  class(df) <- c("vi_integrations", "data.frame")
  df
}

#' Print integration records
#' @param x A \code{vi_integrations} object.
#' @param ... Ignored.
#' @export
print.vi_integrations <- function(x, ...) {
  cat("<vi_integrations>", nrow(x), "records x", ncol(x), "columns\n")
  utils::str(utils::head(as.data.frame(x), 3))
  invisible(x)
}

integration_plot_df <- function(integrations, virus_name) {
  df <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)

  out <- data.frame(
    seqnames1 = df$host_chr,
    start1 = df$host_pos,
    end1 = df$host_pos,
    seqnames2 = rep(virus_name, nrow(df)),
    start2 = df$virus_pos,
    end2 = df$virus_pos,
    host_chr = df$host_chr,
    host_pos = df$host_pos,
    virus_chr = rep(virus_name, nrow(df)),
    virus_pos = df$virus_pos,
    sample = df$sample,
    support = df$support,
    virus_strand = df$virus_strand,
    method = df$method,
    Label = df$sample,
    Depth = df$support,
    Source = df$method,
    ViralStrand = df$virus_strand,
    stringsAsFactors = FALSE
  )

  out
}

split_named_filters <- function(df, subset = NULL) {
  if (is.null(subset) || length(subset) == 0) {
    return(rep(TRUE, nrow(df)))
  }

  keep <- rep(TRUE, nrow(df))
  for (nm in names(subset)) {
    if (!nm %in% colnames(df)) {
      stop("Unknown subset column: ", nm, call. = FALSE)
    }
    keep <- keep & df[[nm]] %in% subset[[nm]]
  }
  keep
}

#' Create an ideogram track specification
#' @param height Track height.
#' @param grid_col Optional named vector of chromosome colors.
#' @param border_col Border color for the ideogram rectangles.
#' @param label_cex Chromosome label text size.
#' @param axis_label_cex Genomic axis label text size.
#' @param axis Which sectors should show genomic axis labels. One of
#'   \code{"virus"}, \code{"all"}, or \code{"none"}.
#' @return A track object for \code{plot_integrations()}.
#' @export
track_ideogram <- function(height = 0.08, grid_col = NULL, border_col = "white",
                           label_cex = 0.52, axis_label_cex = 0.25,
                           axis = c("virus", "all", "none")) {
  axis <- match.arg(axis)
  structure(
    list(
      type = "ideogram",
      height = height,
      grid_col = grid_col,
      border_col = border_col,
      label_cex = label_cex,
      axis_label_cex = axis_label_cex,
      axis = axis
    ),
    class = "vi_track"
  )
}

#' Create a site/scatter track specification
#' @param height Track height.
#' @param split_by Optional column used to split into multiple tracks.
#' @param color Column name or literal color for points.
#' @param size Column name or numeric constant for point size.
#' @param label Optional track label.
#' @param label_cex Track label text size.
#' @param label_col Track label color.
#' @param point_color Default point color when no mapping is provided.
#' @param baseline_col Baseline color.
#' @param subset Optional named list of filters.
#' @param method_col Optional named vector of colors for explicit mapping.
#' @return A track object for \code{plot_integrations()}.
#' @export
track_sites <- function(height = 0.15, split_by = NULL, color = NULL,
                        size = "support", label = NULL,
                        label_cex = 0.8, label_col = "grey30",
                        point_color = "blue", baseline_col = "grey90",
                        subset = NULL, method_col = NULL) {
  structure(
    list(
      type = "sites",
      height = height,
      split_by = split_by,
      color = color,
      size = size,
      label = label,
      label_cex = label_cex,
      label_col = label_col,
      point_color = point_color,
      baseline_col = baseline_col,
      subset = subset,
      method_col = method_col
    ),
    class = "vi_track"
  )
}

#' Create a density track specification
#' @param height Track height.
#' @param split_by Optional column used to split into multiple tracks.
#' @param bins Number of bins per sector.
#' @param fill Fill color for the histogram bars.
#' @param label Optional track label.
#' @param subset Optional named list of filters.
#' @return A track object for \code{plot_integrations()}.
#' @export
track_density <- function(height = 0.12, split_by = NULL, bins = 50,
                          fill = "#4D4D4D", label = NULL, subset = NULL) {
  structure(
    list(
      type = "density",
      height = height,
      split_by = split_by,
      bins = bins,
      fill = fill,
      label = label,
      subset = subset
    ),
    class = "vi_track"
  )
}

#' Create a virus feature track specification
#' @param features Virus feature intervals from \code{virus_features()}.
#' @param height Track height.
#' @param fill Feature fill color or named vector by feature type.
#' @param label Logical. Whether to label features.
#' @param label_cex Feature label text size.
#' @param label_col Feature label color.
#' @param label_min_width Minimum feature width, in base pairs, required for
#'   drawing an internal feature label.
#' @param border_col Feature border color.
#' @return A track object for \code{plot_integrations()}.
#' @export
track_virus_genes <- function(features, height = 0.1, fill = NULL,
                              label = TRUE, label_cex = 0.32,
                              label_col = "grey20",
                              label_min_width = 250,
                              border_col = "white") {
  structure(
    list(
      type = "virus_genes",
      height = height,
      features = normalize_virus_features(features),
      fill = fill,
      label = label,
      label_cex = label_cex,
      label_col = label_col,
      label_min_width = label_min_width,
      border_col = border_col,
      subset = NULL
    ),
    class = "vi_track"
  )
}

#' Create a virus-position density track specification
#' @param height Track height.
#' @param bins Number of bins across the virus genome.
#' @param fill Histogram fill color.
#' @param label Optional track label.
#' @param label_cex Track label text size.
#' @param label_col Track label color.
#' @param subset Optional named list of filters.
#' @return A track object for \code{plot_integrations()}.
#' @export
track_virus_density <- function(height = 0.12, bins = 50, fill = "#4D4D4D",
                                label = NULL, label_cex = 0.8,
                                label_col = "grey30", subset = NULL) {
  structure(
    list(
      type = "virus_density",
      height = height,
      bins = bins,
      fill = fill,
      label = label,
      label_cex = label_cex,
      label_col = label_col,
      subset = subset
    ),
    class = "vi_track"
  )
}

#' Create a link track specification
#' @param color Column name or literal color for links.
#' @param radius Optional link radius.
#' @param lwd Link width.
#' @param default_col Default link color.
#' @param subset Optional named list of filters.
#' @param method_col Optional named vector of colors for explicit mapping.
#' @return A track object for \code{plot_integrations()}.
#' @export
track_links <- function(color = NULL, radius = NULL, lwd = 0.35,
                        default_col = "grey", subset = NULL, method_col = NULL) {
  structure(
    list(
      type = "links",
      color = color,
      radius = radius,
      lwd = lwd,
      default_col = default_col,
      subset = subset,
      method_col = method_col
    ),
    class = "vi_track"
  )
}

prepare_integration_tracks <- function(tracks, plot_df) {
  out <- list()

  for (track in tracks) {
    if (!inherits(track, "vi_track")) {
      stop("Each track must be created by track_ideogram(), track_sites(), track_density(), track_virus_genes(), track_virus_density(), or track_links().", call. = FALSE)
    }

    if (!is.null(track$split_by)) {
      if (!track$split_by %in% colnames(plot_df)) {
        stop("Unknown split_by column: ", track$split_by, call. = FALSE)
      }

      values <- unique(plot_df[[track$split_by]])
      values <- values[!is.na(values)]
      values <- as.character(values)

      for (value in values) {
        clone <- track
        clone$subset <- c(clone$subset, stats::setNames(list(value), track$split_by))
        if (is.null(clone$label)) {
          clone$label <- value
        }
        clone$split_by <- NULL
        out[[length(out) + 1L]] <- clone
      }
    } else {
      out[[length(out) + 1L]] <- track
    }
  }

  out
}

draw_integration_track <- function(track, plot_df, cfg) {
  data <- plot_df[split_named_filters(plot_df, track$subset), , drop = FALSE]

  if (track$type == "ideogram") {
    grid_col <- track$grid_col
    if (is.null(grid_col) && !is.null(cfg$grid_col)) {
      grid_col <- cfg$grid_col
    }

    return(draw_ideogram(
      height = track$height,
      cfg = cfg,
      grid_col = grid_col,
      border_col = track$border_col,
      label_cex = track$label_cex,
      axis_label_cex = track$axis_label_cex,
      axis = track$axis
    ))
  }

  if (track$type == "sites") {
    if (!is.null(track$color) && length(track$color) == 1L && track$color %in% colnames(data)) {
      data$Source <- as.character(data[[track$color]])
      track_method_col <- NULL
      point_color <- track$point_color
    } else if (!is.null(track$color)) {
      data$Source <- rep("integration", nrow(data))
      track_method_col <- if (!is.null(track$method_col)) {
        track$method_col
      } else {
        stats::setNames(track$color, "integration")
      }
      point_color <- track$point_color
    } else {
      track_method_col <- track$method_col
      point_color <- track$point_color
    }

    if (!is.null(track$size)) {
      if (is.character(track$size) && length(track$size) == 1L && track$size %in% colnames(data)) {
        data$Depth <- suppressWarnings(as.numeric(data[[track$size]]))
      } else if (is.numeric(track$size) && length(track$size) == 1L) {
        data$Depth <- rep(as.numeric(track$size), nrow(data))
      }
    }

    return(draw_scatter(
      data = data,
      height = track$height,
      cfg = cfg,
      track_label = track$label,
      track_label_cex = track$label_cex,
      track_label_col = track$label_col,
      method_col = track_method_col,
      point_color = point_color,
      baseline_col = track$baseline_col
    ))
  }

  if (track$type == "density") {
    return(draw_histogram(
      data = data,
      height = track$height,
      cfg = cfg,
      track_label = track$label,
      bins = track$bins,
      col = track$fill
    ))
  }

  if (track$type == "virus_genes") {
    return(draw_virus_genes(
      features = track$features,
      height = track$height,
      cfg = cfg,
      fill = track$fill,
      label = track$label,
      label_cex = track$label_cex,
      label_col = track$label_col,
      label_min_width = track$label_min_width,
      border_col = track$border_col
    ))
  }

  if (track$type == "virus_density") {
    return(draw_virus_density(
      data = data,
      height = track$height,
      cfg = cfg,
      bins = track$bins,
      col = track$fill,
      track_label = track$label,
      track_label_cex = track$label_cex,
      track_label_col = track$label_col
    ))
  }

  if (track$type == "links") {
    if (!is.null(track$color) && length(track$color) == 1L && track$color %in% colnames(data)) {
      data$Source <- as.character(data[[track$color]])
      link_method_col <- NULL
      default_col <- track$default_col
    } else if (!is.null(track$color)) {
      data$Source <- rep("integration", nrow(data))
      link_method_col <- if (!is.null(track$method_col)) {
        track$method_col
      } else {
        stats::setNames(track$color, "integration")
      }
      default_col <- track$default_col
    } else {
      link_method_col <- track$method_col
      default_col <- track$default_col
    }

    return(draw_link(
      link_data = data,
      cfg = cfg,
      radius = track$radius,
      lwd = track$lwd,
      method_col = link_method_col,
      default_col = default_col
    ))
  }

  stop("Unsupported track type: ", track$type, call. = FALSE)
}

#' Plot virus-host integrations
#'
#' @param integrations Integration records.
#' @param host Host genome specification. Can be a built-in host name, a host
#'   genome table, or a \code{vi_host_genome} object.
#' @param chrom_file Optional path to a host chromosome-size table.
#' @param virus Virus genome specification, typically a named numeric vector of
#'   length 1 such as \code{c(HBV = 3215)}.
#' @param tracks A list of track objects created by \code{track_*()}.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @param clear Logical. Whether to clear the current circlize device first.
#' @param draw Logical. Whether to draw the plot immediately.
#' @return An object of class \code{vi_integration_plot}.
#' @export
plot_integrations <- function(integrations, host = "hg38", virus,
                              tracks = NULL, visual_ratio = 0.1,
                              clear = TRUE, chrom_file = NULL,
                              draw = TRUE) {
  integrations <- as_integrations(integrations)
  virus_obj <- if (inherits(virus, "vi_virus_genome")) {
    virus
  } else {
    virus_genome(virus)
  }
  host_obj <- if (inherits(host, "vi_host_genome")) {
    host
  } else {
    host_genome(host = host, chrom_file = chrom_file)
  }
  cfg <- create_config(
    host = host_obj,
    virus_name = virus_obj$name,
    virus_length = virus_obj$length,
    visual_ratio = visual_ratio
  )
  cfg$grid_col <- host_obj$colors

  plot_df <- integration_plot_df(integrations, virus_name = virus_obj$name)
  if (is.null(tracks)) {
    tracks <- list(
      track_ideogram(),
      track_sites(),
      track_links()
    )
  }

  tracks <- prepare_integration_tracks(tracks, plot_df)

  out <- structure(
    list(
      integrations = integrations,
      host = host_obj,
      virus = virus_obj,
      cfg = cfg,
      plot_df = plot_df,
      tracks = tracks,
      drawn = FALSE
    ),
    class = "vi_integration_plot"
  )

  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) {
    stop("draw must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(draw)) {
    out <- draw_integration_plot(out, clear = clear)
    out$drawn <- TRUE
  }

  out
}

#' @export
print.vi_integration_plot <- function(x, ...) {
  if (isTRUE(x$drawn)) {
    invisible(x)
  } else if (!is.null(x$plot_df) && !is.null(x$cfg)) {
    return(draw_integration_plot(x))
  }
  invisible(x)
}

draw_integration_plot <- function(x, clear = TRUE) {
  if (!inherits(x, "vi_integration_plot")) {
    stop("x must be a vi_integration_plot object.", call. = FALSE)
  }

  if (isTRUE(clear)) {
    circlize::circos.clear()
  }

  n_sectors <- nrow(x$cfg$data)
  gaps <- rep(1, n_sectors)
  virus_idx <- which(x$cfg$data$chr == x$cfg$virus_name)
  if (length(virus_idx) > 0) {
    gaps[virus_idx] <- 10
  }

  circlize::circos.par(
    start.degree = 90,
    gap.degree = gaps,
    cell.padding = c(0, 0, 0, 0),
    points.overflow.warning = FALSE
  )

  circlize::circos.genomicInitialize(
    x$cfg$data[, c("chr", "start", "end"), drop = FALSE],
    plotType = NULL,
    sector.width = x$cfg$widths
  )

  for (track in x$tracks) {
    draw_integration_track(track, x$plot_df, x$cfg)
  }

  legend_spec <- integration_method_legend_spec(x$tracks, x$plot_df)
  if (!is.null(legend_spec)) {
    draw_method_legend(legend_spec)
  }

  x$drawn <- TRUE
  invisible(x)
}

integration_method_legend_spec <- function(tracks, plot_df) {
  for (track in tracks) {
    if (track$type %in% c("sites", "links") &&
        is.character(track$color) &&
        length(track$color) == 1L &&
        track$color %in% colnames(plot_df)) {
      method_col <- if (!is.null(track$method_col)) track$method_col else NULL
      return(resolve_method_colors(plot_df[[track$color]], method_col = method_col))
    }
  }

  NULL
}

legacy_layout_to_tracks <- function(layout_list) {
  out <- list()

  for (task in layout_list) {
    if (task$type == "ideogram") {
      out[[length(out) + 1L]] <- track_ideogram(
        height = if (is.null(task$height)) 0.08 else task$height,
        grid_col = task$grid_col,
        border_col = if (is.null(task$border_col)) "white" else task$border_col,
        label_cex = if (is.null(task$label_cex)) 0.52 else task$label_cex,
        axis_label_cex = if (is.null(task$axis_label_cex)) 0.25 else task$axis_label_cex,
        axis = if (is.null(task$axis)) "virus" else task$axis
      )
    } else if (task$type == "scatter") {
      out[[length(out) + 1L]] <- track_sites(
        height = if (is.null(task$height)) 0.15 else task$height,
        split_by = NULL,
        color = task$color,
        label = task$sample_label,
        point_color = if (is.null(task$point_color)) "blue" else task$point_color,
        baseline_col = if (is.null(task$baseline_col)) "grey90" else task$baseline_col,
        subset = if (is.null(task$sample_label)) NULL else stats::setNames(list(task$sample_label), "sample"),
        method_col = task$method_col
      )
    } else if (task$type == "histogram") {
      out[[length(out) + 1L]] <- track_density(
        height = if (is.null(task$height)) 0.12 else task$height,
        split_by = NULL,
        bins = if (is.null(task$bins)) 50 else task$bins,
        fill = if (is.null(task$col)) "#4D4D4D" else task$col,
        label = task$sample_label,
        subset = if (is.null(task$sample_label)) NULL else stats::setNames(list(task$sample_label), "sample")
      )
    } else if (task$type == "links") {
      out[[length(out) + 1L]] <- track_links(
        color = task$color,
        radius = task$radius,
        lwd = if (is.null(task$lwd)) 0.35 else task$lwd,
        default_col = if (is.null(task$default_col)) "grey" else task$default_col,
        method_col = task$method_col
      )
    }
  }

  out
}

#' Plot oncoprint with annotation tracks
#'
#' @param maf MAF object.
#' @param genes The gene names or the number of genes to display.
#' @param annotations Character vector or list of annotation variables.
#' @param numeric_bins Number of bins for continuous annotation variables.
#' @param annotation_colors Named list of colors for annotations.
#' @param sort_by Annotation variable used to sort samples.
#' @param annotation_levels Named list of explicit level orders.
#' @return A combined plot object.
#' @export
plot_oncoprint <- function(maf, genes = 20, annotations = NULL, numeric_bins = 5,
                           annotation_colors = NULL, sort_by = NULL,
                           annotation_levels = NULL) {
  components <- aplotExtra::oncoplot_components(maf, genes)

  sample_order <- oncoplot_resolve_sample_order(
    maf,
    genes,
    sort_by = sort_by,
    clinical_order = annotation_levels,
    components = components
  )

  p_main <- oncoplot_apply_sample_order(components$main, sample_order)
  p_top <- oncoplot_apply_sample_order(components$sample, sample_order)
  p_right <- components$gene
  p_spacer <- ggplot2::ggplot() + ggfun::theme_transparent()

  pp <- p_main |>
    aplot::insert_top(p_spacer, height = 0.02) |>
    aplot::insert_top(p_top, height = 0.2) |>
    aplot::insert_right(p_right, width = 0.2)

  tracks <- oncoplot_clinical_track(
    maf,
    genes = genes,
    clinical_vars = annotations,
    n_breaks = numeric_bins,
    colors = annotation_colors,
    sample_order = sample_order,
    clinical_order = annotation_levels,
    components = components
  )

  if (length(tracks) > 0) {
    for (var in names(tracks)) {
      pp <- aplot::insert_bottom(pp, tracks[[var]], height = 0.05)
    }
  }

  pp
}

#' Visualize Host-Virus Integration Events
#'
#' Draw a host-virus circos plot on the current graphics device from an
#' integration table and host chromosome-size information. This legacy wrapper
#' delegates to \code{plot_integrations()} while preserving the original
#' \code{visualize_viral_integration()} interface.
#'
#' @param input_file Path to a tab-delimited integration table.
#' @param host Character scalar specifying a standard host genome name or a
#'   UCSC assembly name. Common aliases such as \code{"human"} and
#'   \code{"mouse"} are supported.
#' @param chrom_file Optional path to a host chromosome-size table with columns
#'   \code{chr}, \code{start}, and \code{end}. Use this when \code{host} is not
#'   supplied or when a custom host genome is needed.
#' @param virus_name Name of the virus sequence.
#' @param virus_length Length of the virus sequence in base pairs.
#' @param layout_list A list of legacy track definitions. Supported
#'   \code{type} values are \code{"ideogram"}, \code{"scatter"},
#'   \code{"histogram"}, and \code{"links"}.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @param clear Logical. Whether to clear the existing circlize plot before
#'   drawing.
#' @return Invisibly returns a list with \code{cfg}, \code{gi}, and
#'   \code{data}.
#' @export
visualize_viral_integration <- function(input_file,
                                        host = NULL,
                                        chrom_file = NULL,
                                        virus_name,
                                        virus_length,
                                        layout_list,
                                        visual_ratio = 0.1,
                                        clear = TRUE) {
  integrations <- read_integrations(input_file)
  tracks <- legacy_layout_to_tracks(layout_list)

  plot_obj <- plot_integrations(
    integrations = integrations,
    host = host,
    virus = stats::setNames(virus_length, virus_name),
    tracks = tracks,
    visual_ratio = visual_ratio,
    clear = clear,
    chrom_file = chrom_file
  )

  gi <- create_gi_from_table(input_file = input_file, cfg = plot_obj$cfg)

  invisible(list(cfg = plot_obj$cfg, gi = gi, data = plot_obj$plot_df))
}
