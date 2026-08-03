validate_virus_info <- function(virus_name, virus_length, host_chr) {
  if (missing(virus_name) || is.null(virus_name) || length(virus_name) != 1 ||
      is.na(virus_name) || !nzchar(virus_name)) {
    stop("virus_name must be a single non-empty sequence name.")
  }

  virus_name <- trimws(gsub("(?i)^chr", "", as.character(virus_name)))

  if (!is.numeric(virus_length) || length(virus_length) != 1 ||
      is.na(virus_length) || virus_length <= 0) {
    stop("virus_length must be a single positive number.", call. = FALSE)
  }

  virus_length <- as.numeric(virus_length)

  if (virus_name %in% host_chr) {
    stop("virus_name conflicts with a host chromosome name: ", virus_name, call. = FALSE)
  }

  list(
    virus_name = virus_name,
    virus_length = virus_length
  )
}

#' Create Plotting Configuration
#'
#' @param host Character scalar specifying a built-in host genome.
#' @param chrom_file Optional path to a host chromosome-size table.
#' @param virus_name Name of the virus sequence.
#' @param virus_length Length of the virus sequence in base pairs.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @return A plotting configuration list.
#' @keywords internal
create_config <- function(host = NULL,
                          chrom_file = NULL,
                          virus_name,
                          virus_length,
                          visual_ratio = 0.1) {
  host_df <- resolve_host_chrom_sizes(host = host, chrom_file = chrom_file)
  virus_info <- validate_virus_info(virus_name, virus_length, host_df$chr)
  virus_name <- virus_info$virus_name
  virus_length <- virus_info$virus_length

  if (!is.numeric(visual_ratio) || length(visual_ratio) != 1 ||
      is.na(visual_ratio) || visual_ratio <= 0 || visual_ratio >= 1) {
    stop("visual_ratio must be a single number between 0 and 1.")
  }

  chrom_df <- rbind(
    host_df,
    data.frame(chr = virus_name, start = 0, end = virus_length, stringsAsFactors = FALSE)
  )

  widths <- numeric(nrow(chrom_df))
  names(widths) <- chrom_df$chr
  widths[virus_name] <- visual_ratio

  host_mask <- chrom_df$chr != virus_name
  host_len <- chrom_df$end[host_mask] - chrom_df$start[host_mask]
  widths[host_mask] <- host_len / sum(host_len) * (1 - visual_ratio)

  list(
    data = chrom_df,
    widths = widths,
    virus_name = virus_name
  )
}

#' Convert an Integration Table to GInteractions
#'
#' @param input_file Path to an integration table.
#' @param cfg Plotting configuration returned by \code{create_config()}.
#' @return A \code{GInteractions} object.
#' @export
create_gi_from_table <- function(input_file, cfg) {
  if (!file.exists(input_file)) {
    stop("File not found: ", input_file)
  }

  raw <- as.data.frame(as_integrations(input_file), stringsAsFactors = FALSE)

  host_chr <- raw$host_chr
  virus_chr <- rep(cfg$virus_name, nrow(raw))
  host_pos <- as.numeric(raw$host_pos)
  virus_pos <- as.numeric(raw$virus_pos)

  keep <- !is.na(host_chr) & host_chr %in% cfg$data$chr & !is.na(host_pos) & !is.na(virus_pos)
  if (!all(keep)) {
    warning(sum(!keep), " rows were skipped because chromosome names or position values were invalid.")
  }

  raw <- raw[keep, , drop = FALSE]
  host_chr <- host_chr[keep]
  virus_chr <- virus_chr[keep]
  host_pos <- host_pos[keep]
  virus_pos <- virus_pos[keep]

  if (nrow(raw) == 0) {
    stop("No valid integration records were found in data.txt.")
  }

  gr_host <- GenomicRanges::GRanges(
    seqnames = host_chr,
    ranges = IRanges::IRanges(host_pos, width = 1)
  )
  gr_virus <- GenomicRanges::GRanges(
    seqnames = virus_chr,
    ranges = IRanges::IRanges(virus_pos, width = 1),
    strand = raw$virus_strand
  )

  suppressWarnings(
    InteractionSet::GInteractions(
      gr_host,
      gr_virus,
      mode = "strict",
      Depth = as.numeric(raw$support),
      Label = raw$sample,
      Source = raw$method,
      ViralStrand = raw$virus_strand
    )
  )
}

#' Draw Ideogram Track in Genomic Mode
#'
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param grid_col Optional named vector of sector colors.
#' @param border_col Border color for ideogram rectangles.
#' @param label_cex Chromosome label text size.
#' @param axis_label_cex Genomic axis label text size.
#' @param axis Which sectors should show genomic axis labels. One of
#'   \code{"virus"}, \code{"all"}, or \code{"none"}.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_ideogram <- function(height, cfg, grid_col = NULL, border_col = "white",
                          label_cex = 0.52, axis_label_cex = 0.25,
                          axis = c("virus", "all", "none")) {
  axis <- match.arg(axis)
  ideogram_df <- cfg$data[, c("chr", "start", "end"), drop = FALSE]
  grid_col <- resolve_ideogram_colors(cfg, grid_col = grid_col)

  circlize::circos.genomicTrackPlotRegion(
    data = ideogram_df,
    ylim = c(-0.2, 1.5),
    bg.border = NA,
    track.height = height,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      xlim <- circlize::CELL_META$cell.xlim
      fill_col <- grid_col[chr]

      circlize::circos.genomicRect(
        region = region,
        value = value,
        ybottom = 0,
        ytop = 1,
        col = fill_col,
        border = border_col
      )

      circlize::circos.text(
        mean(xlim), 1.25, chr,
        facing = "bending.inside",
        niceFacing = TRUE,
        adj = c(0.5, -0.2),
        cex = ideogram_label_cex(chr, xlim, cfg, label_cex)
      )

      if (show_ideogram_axis(chr, cfg, axis)) {
        axis_labels <- make_axis_labels(xlim, use_kb = chr == cfg$virus_name)
        circlize::circos.axis(
          h = 1.1,
          major.at = axis_labels$at,
          labels = axis_labels$labels,
          labels.cex = axis_label_cex,
          direction = "outside",
          major.tick.length = 0.02
        )
      }
    }
  )

  invisible(NULL)
}

#' Draw Scatter Track in Genomic Mode
#'
#' @param data Data frame converted from a \code{GInteractions} object.
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param track_label Optional track label to draw on the virus sector.
#' @param track_label_cex Track label text size.
#' @param track_label_col Track label color.
#' @param method_col Optional named vector of colors for different methods.
#' @param point_color Default point color.
#' @param baseline_col Baseline color for the scatter track.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_scatter <- function(data, height, cfg, track_label = NULL,
                         track_label_cex = 0.8,
                         track_label_col = "grey30",
                         method_col = NULL,
                         point_color = "blue", baseline_col = "grey90") {
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  max_depth <- if ("Depth" %in% colnames(data)) max(log10(data$Depth + 1), na.rm = TRUE) else 1
  if (!is.finite(max_depth) || max_depth <= 0) {
    max_depth <- 1
  }

  depth <- if ("Depth" %in% colnames(data)) data$Depth else rep(10, nrow(data))
  source_col <- if ("Source" %in% colnames(data)) as.character(data$Source) else rep("integration", nrow(data))
  col_map <- resolve_method_colors(source_col, method_col = method_col)
  pt_col <- col_map[source_col]
  pt_col[is.na(pt_col)] <- point_color
  pt_col <- grDevices::adjustcolor(pt_col, alpha.f = 0.85)
  pt_cex <- log10(depth + 1) / max_depth * 1.5 + 0.5
  cex_range <- range(pt_cex, na.rm = TRUE)
  if (!all(is.finite(cex_range)) || diff(cex_range) == 0) {
    pt_y <- rep(0.12, length(pt_cex))
  } else {
    pt_y <- 0.08 + (pt_cex - cex_range[1]) / diff(cex_range) * (0.16 - 0.08)
  }

  scatter_df <- data.frame(
    chr = as.character(data$seqnames1),
    start = data$start1,
    end = data$end1,
    y = pt_y,
    cex = pt_cex,
    pt_col = pt_col,
    stringsAsFactors = FALSE
  )
  scatter_df <- ensure_virus_track_row(
    df = scatter_df,
    virus_name = cfg$virus_name,
    virus_length = cfg$data$end[cfg$data$chr == cfg$virus_name][1],
    defaults = list(y = 0.72, cex = 0, pt_col = NA_character_)
  )

  circlize::circos.genomicTrackPlotRegion(
    data = scatter_df,
    ylim = c(0, 1),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0, 0), col = baseline_col)

      if (nrow(region) > 0) {
        circlize::circos.genomicPoints(
          region,
          value,
          numeric.column = "y",
          pch = 19,
          col = value$pt_col,
          cex = value$cex
        )
      }

      if (!is.null(track_label) && identical(chr, cfg$virus_name)) {
        draw_track_label(
          track_label,
          y = 0.72,
          cex = track_label_cex,
          col = track_label_col
        )
      }
    }
  )

  invisible(NULL)
}

#' Draw Virus Feature Track in Genomic Mode
#'
#' @param features Virus feature intervals from \code{virus_features()}.
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param fill Feature fill color or named vector by feature type.
#' @param label Logical. Whether to label features.
#' @param label_cex Feature label text size.
#' @param label_col Feature label color.
#' @param label_min_width Minimum feature width, in base pairs, required for
#'   drawing an internal feature label.
#' @param border_col Feature border color.
#' @return Invisibly returns \code{NULL}.
draw_virus_genes <- function(features, height, cfg, fill = NULL, label = TRUE,
                             label_cex = 0.32, label_col = "grey20",
                             label_min_width = 250,
                             border_col = "white") {
  virus_length <- cfg$data$end[cfg$data$chr == cfg$virus_name][1]
  feature_df <- layout_virus_features(features, virus_length = virus_length)
  feature_df$chr <- cfg$virus_name
  feature_df$ymin <- feature_df$level - 0.85
  feature_df$ymax <- feature_df$level - 0.15

  feature_fill <- if (is.null(fill)) {
    normalize_named_colors(grDevices::hcl.colors(length(unique(feature_df$type)), "Set 3"), unique(feature_df$type))
  } else {
    normalize_named_colors(fill, unique(feature_df$type))
  }
  feature_df$fill <- unname(feature_fill[feature_df$type])

  feature_region <- feature_df[, c("chr", "start", "end"), drop = FALSE]
  feature_region$start <- pmax(0, pmin(feature_region$start, virus_length))
  feature_region$end <- pmax(0, pmin(feature_region$end, virus_length))
  feature_df$label <- ifelse(
    feature_df$end - feature_df$start >= label_min_width,
    feature_df$feature,
    NA_character_
  )
  feature_value <- feature_df[, c("ymin", "ymax", "fill", "feature", "label"), drop = FALSE]

  circlize::circos.genomicTrackPlotRegion(
    data = cbind(feature_region, feature_value),
    ylim = c(0, max(feature_df$level)),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      if (!identical(chr, cfg$virus_name)) {
        return(NULL)
      }

      circlize::circos.genomicRect(
        region,
        value,
        ybottom = value$ymin,
        ytop = value$ymax,
        col = value$fill,
        border = border_col
      )

      if (isTRUE(label) && any(!is.na(value$label))) {
        label_region <- region[!is.na(value$label), , drop = FALSE]
        label_value <- value[!is.na(value$label), , drop = FALSE]
        circlize::circos.genomicText(
          label_region,
          label_value,
          y = (label_value$ymin + label_value$ymax) / 2,
          labels.column = "label",
          facing = "bending.inside",
          niceFacing = TRUE,
          cex = label_cex,
          col = label_col
        )
      }
    }
  )

  invisible(NULL)
}

show_ideogram_axis <- function(chr, cfg, axis) {
  if (identical(axis, "none")) {
    return(FALSE)
  }
  if (identical(axis, "all")) {
    return(TRUE)
  }

  identical(chr, cfg$virus_name)
}

ideogram_label_cex <- function(chr, xlim, cfg, label_cex) {
  if (identical(chr, cfg$virus_name)) {
    return(max(label_cex, 0.78))
  }

  chr_width <- diff(xlim)
  if (!is.finite(chr_width)) {
    return(label_cex)
  }
  if (chr_width < 75e6) {
    return(label_cex * 0.75)
  }
  if (chr_width < 125e6) {
    return(label_cex * 0.88)
  }

  label_cex
}

#' Draw Virus Density Track in Genomic Mode
#'
#' @param data Data frame converted from integration records.
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param bins Number of bins across the virus genome.
#' @param col Histogram fill color.
#' @param track_label Optional track label.
#' @param track_label_cex Track label text size.
#' @param track_label_col Track label color.
#' @return Invisibly returns \code{NULL}.
draw_virus_density <- function(data, height, cfg, bins = 50, col = "#4D4D4D",
                               track_label = NULL,
                               track_label_cex = 0.8,
                               track_label_col = "grey30") {
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  hist_df <- make_virus_histogram_data(data = data, cfg = cfg, bins = bins)
  max_count <- max(hist_df$count, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) {
    max_count <- 1
  }

  circlize::circos.genomicTrackPlotRegion(
    data = hist_df,
    ylim = c(0, max_count),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      if (!identical(chr, cfg$virus_name)) {
        return(NULL)
      }

      circlize::circos.genomicRect(
        region,
        value,
        ybottom = 0,
        ytop = value$count,
        col = col,
        border = NA
      )

      if (!is.null(track_label)) {
        draw_track_label(
          track_label,
          y = max_count * 0.8,
          cex = track_label_cex,
          col = track_label_col
        )
      }
    }
  )

  invisible(NULL)
}

make_virus_histogram_data <- function(data, cfg, bins = 50) {
  if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) || bins < 1) {
    stop("bins must be a positive integer.", call. = FALSE)
  }
  bins <- as.integer(bins)

  virus_length <- cfg$data$end[cfg$data$chr == cfg$virus_name][1]
  breaks <- seq(0, virus_length, length.out = bins + 1)
  pos <- data$virus_pos[!is.na(data$virus_pos)]
  counts <- tabulate(
    findInterval(pos, breaks, rightmost.closed = TRUE, all.inside = TRUE),
    nbins = bins
  )

  data.frame(
    chr = cfg$virus_name,
    start = breaks[-length(breaks)],
    end = breaks[-1],
    count = counts,
    stringsAsFactors = FALSE
  )
}

#' Draw Link Layer in Genomic Mode
#'
#' @param link_data Data frame converted from a \code{GInteractions} object.
#' @param cfg Plotting configuration list.
#' @param radius Optional link radius.
#' @param lwd Link line width.
#' @param method_col Optional named vector of colors for different methods.
#' @param default_col Fallback link color.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_link <- function(link_data, cfg, radius = NULL, lwd = 0.35, method_col = NULL,
                      default_col = "grey") {
  if (is.null(link_data) || nrow(link_data) == 0) {
    return(invisible(NULL))
  }

  source_col <- if ("Source" %in% colnames(link_data)) as.character(link_data$Source) else rep("integration", nrow(link_data))
  col_map <- resolve_method_colors(source_col, method_col = method_col)

  for (i in seq_len(nrow(link_data))) {
    link_col <- col_map[source_col[i]]
    if (is.na(link_col)) {
      link_col <- default_col
    }

    region1 <- data.frame(
      chr = as.character(link_data$seqnames1[i]),
      start = link_data$start1[i],
      end = link_data$end1[i]
    )
    region2 <- data.frame(
      chr = as.character(link_data$seqnames2[i]),
      start = link_data$start2[i],
      end = link_data$end2[i]
    )

    if (is.null(radius)) {
      circlize::circos.genomicLink(region1, region2, col = link_col, border = link_col, lwd = lwd)
    } else {
      circlize::circos.genomicLink(region1, region2, rou = radius, col = link_col, border = link_col, lwd = lwd)
    }
  }

  invisible(NULL)
}

#' Draw Histogram Track in Genomic Mode
#'
#' @param data Data frame converted from a \code{GInteractions} object.
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param track_label Optional track label to draw on the HBV sector.
#' @param bins Number of genomic bins per sector.
#' @param col Histogram fill color.
#' @return Invisibly returns \code{NULL}.
draw_histogram <- function(data, height, cfg, track_label = NULL, bins = NULL, col = "#4D4D4D") {
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  bins <- if (is.null(bins)) 50 else bins
  if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) || bins < 1) {
    stop("bins must be a positive integer.")
  }
  bins <- as.integer(bins)

  hist_df <- make_histogram_data(data = data, cfg = cfg, bins = bins)
  hist_df <- ensure_virus_track_row(
    df = hist_df,
    virus_name = cfg$virus_name,
    virus_length = cfg$data$end[cfg$data$chr == cfg$virus_name][1],
    defaults = list(count = 0)
  )
  if (nrow(hist_df) == 0) {
    return(invisible(NULL))
  }

  max_count <- max(hist_df$count, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) {
    max_count <- 1
  }

  circlize::circos.genomicTrackPlotRegion(
    data = hist_df,
    ylim = c(0, max_count),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      if (nrow(region) > 0) {
        circlize::circos.genomicRect(
          region,
          value,
          ybottom = 0,
          ytop = value$count,
          col = col,
          border = NA
        )
      }

      if (!is.null(track_label) && identical(chr, cfg$virus_name)) {
        draw_track_label(track_label, y = max_count * 0.8)
      }
    }
  )

  invisible(NULL)
}

make_histogram_data <- function(data, cfg, bins = 50) {
  out <- vector("list", nrow(cfg$data))

  for (i in seq_len(nrow(cfg$data))) {
    chr <- cfg$data$chr[i]
    chr_start <- cfg$data$start[i]
    chr_end <- cfg$data$end[i]
    chr_pos <- data$start1[as.character(data$seqnames1) == chr]
    chr_pos <- chr_pos[!is.na(chr_pos)]

    breaks <- seq(chr_start, chr_end, length.out = bins + 1)
    if (length(unique(breaks)) < 2) {
      next
    }

    counts <- tabulate(
      findInterval(chr_pos, breaks, rightmost.closed = TRUE, all.inside = TRUE),
      nbins = bins
    )

    out[[i]] <- data.frame(
      chr = chr,
      start = breaks[-length(breaks)],
      end = breaks[-1],
      count = counts,
      stringsAsFactors = FALSE
    )
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    return(data.frame(chr = character(), start = numeric(), end = numeric(), count = numeric()))
  }

  do.call(rbind, out)
}

draw_method_legend <- function(method_col) {
  if (is.null(method_col) || length(method_col) == 0) {
    return(invisible(NULL))
  }

  graphics::legend(
    "topright",
    legend = names(method_col),
    col = unname(method_col),
    pch = 19,
    title = "Methods",
    bty = "n",
    cex = 0.8,
    inset = 0.02
  )

  invisible(NULL)
}

get_method_legend_spec <- function(layout_list, plot_df) {
  for (task in layout_list) {
    if (task$type %in% c("scatter", "links")) {
      methods <- if ("Source" %in% colnames(plot_df)) as.character(plot_df$Source) else "integration"
      return(resolve_method_colors(methods, method_col = task$method_col))
    }
  }

  NULL
}

filter_plot_data <- function(plot_df, sample_label = NULL) {
  if (is.null(sample_label)) {
    return(plot_df)
  }

  plot_df[plot_df$Label %in% sample_label, , drop = FALSE]
}

ensure_virus_track_row <- function(df, virus_name, virus_length, defaults = list()) {
  if (is.null(virus_name) || is.na(virus_name) || !nzchar(virus_name)) {
    return(df)
  }

  if (virus_name %in% df$chr) {
    return(df)
  }

  extra <- data.frame(
    chr = virus_name,
    start = 0,
    end = virus_length,
    stringsAsFactors = FALSE
  )

  for (nm in names(defaults)) {
    extra[[nm]] <- defaults[[nm]]
  }

  missing_cols <- setdiff(colnames(df), colnames(extra))
  for (nm in missing_cols) {
    extra[[nm]] <- NA
  }

  extra <- extra[, colnames(df), drop = FALSE]
  rbind(df, extra)
}

draw_track_label <- function(label, y, x_frac = 0.5, cex = 0.8, col = "grey30") {
  xlim <- circlize::CELL_META$cell.xlim
  x <- xlim[1] + diff(xlim) * x_frac

  circlize::circos.text(
    x = x,
    y = y,
    labels = label,
    facing = "inside",
    niceFacing = TRUE,
    adj = c(0.5, 0.5),
    cex = cex,
    col = col
  )
}

resolve_ideogram_colors <- function(cfg, grid_col = NULL) {
  chr <- as.character(cfg$data$chr)
  if (is.null(grid_col)) {
    grid_col <- stats::setNames(grDevices::hcl.colors(length(chr), "Dark 3"), chr)
  } else {
    grid_col <- normalize_named_colors(grid_col, chr)
  }

  if (cfg$virus_name %in% names(grid_col)) {
    grid_col[cfg$virus_name] <- "#595959"
  }

  grid_col
}

resolve_method_colors <- function(methods, method_col = NULL) {
  methods <- sort(unique(as.character(methods)))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  if (length(methods) == 0) {
    methods <- "integration"
  }

  if (is.null(method_col)) {
    base_col <- c(
      Pos = "#C77C7C",
      PB = "#7B9ACC",
      integration = "#7FB77E"
    )
    method_col <- base_col[methods]
    missing <- is.na(method_col)
    if (any(missing)) {
      method_col[missing] <- grDevices::hcl.colors(sum(missing), "Dark 3")
    }
    return(stats::setNames(unname(method_col), methods))
  }

  normalize_named_colors(method_col, methods)
}

normalize_named_colors <- function(colors, keys) {
  if (is.null(colors)) {
    return(stats::setNames(rep(NA_character_, length(keys)), keys))
  }

  if (!is.null(names(colors))) {
    matched <- colors[keys]
    missing <- is.na(matched)
    if (any(missing)) {
      matched[missing] <- rep(unname(colors)[1], sum(missing))
    }
    return(stats::setNames(unname(matched), keys))
  }

  if (length(colors) == 1) {
    return(stats::setNames(rep(colors, length(keys)), keys))
  }

  stats::setNames(rep(colors, length.out = length(keys)), keys)
}

validate_layout_list <- function(layout_list) {
  if (is.null(layout_list) || length(layout_list) == 0) {
    stop("layout_list must contain at least one track definition.")
  }

  invalid_obj_idx <- which(!vapply(layout_list, is.list, logical(1)))
  if (length(invalid_obj_idx) > 0) {
    stop("Every layout item must be a list. Invalid index: ",
         paste(invalid_obj_idx, collapse = ", "))
  }

  invalid_idx <- which(vapply(layout_list, function(x) is.null(x$type) || !nzchar(x$type), logical(1)))
  if (length(invalid_idx) > 0) {
    stop("Every layout item must include a non-empty type. Invalid index: ",
         paste(invalid_idx, collapse = ", "))
  }

  valid_types <- c("ideogram", "scatter", "histogram", "links")
  invalid_type_idx <- which(vapply(layout_list, function(x) !x$type %in% valid_types, logical(1)))
  if (length(invalid_type_idx) > 0) {
    stop("Unsupported layout type at index: ",
         paste(invalid_type_idx, collapse = ", "),
         ". Supported types: ", paste(valid_types, collapse = ", "),
         call. = FALSE)
  }

  scatter_missing_label <- which(vapply(
    layout_list,
    function(x) identical(x$type, "scatter") && (is.null(x$sample_label) || !nzchar(x$sample_label)),
    logical(1)
  ))
  if (length(scatter_missing_label) > 0) {
    stop("Each scatter layout item must include a non-empty sample_label. Invalid index: ",
         paste(scatter_missing_label, collapse = ", "),
         call. = FALSE)
  }

  invisible(NULL)
}

#' @importFrom utils tail
#' @noRd
make_axis_labels <- function(xlim, use_kb = FALSE) {
  axis_unit <- if (use_kb) "kb" else "Mb"
  axis_scale <- if (use_kb) 1000 else 1e6
  axis_len <- diff(xlim) / axis_scale
  tick_step <- if (use_kb) {
    if (axis_len <= 10) 1 else if (axis_len <= 50) 10 else 50
  } else {
    50
  }

  label_values <- seq(0, axis_len, by = tick_step)
  major_values <- label_values
  major_labels <- paste0(label_values, " ", axis_unit)

  if (length(major_values) == 0 || tail(major_values, 1) < axis_len) {
    major_values <- c(major_values, axis_len)
    major_labels <- c(major_labels, "")
  }

  list(
    at = xlim[1] + major_values * axis_scale,
    labels = major_labels,
    unit = axis_unit
  )
}
