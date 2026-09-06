#' Create a virus gene track specification
#'
#' @param features Virus feature intervals. Accepts a data frame or
#'   \code{vi_virus_features} object with \code{feature}, \code{start},
#'   \code{end}, and optional \code{type} columns.
#' @param height Track height.
#' @param fill Feature fill color, unnamed palette, or named vector by feature
#'   name or feature type. By default, each distinct feature name gets its own
#'   color.
#' @param label Logical. Whether to label features.
#' @param label_cex Feature label text size.
#' @param label_col Feature label color.
#' @param label_min_width Minimum feature width, in base pairs, required for a
#'   label to be drawn.
#' @param border_col Feature border color. Defaults to \code{NA} so split
#'   circular features do not show a seam at the origin.
#' @return A \code{vi_virus_circos_track} object.
#' @export
virus_circos_track_genes <- function(features, height = 0.12, fill = NULL,
                                     label = TRUE, label_cex = 0.55,
                                     label_col = "grey20",
                                     label_min_width = 250,
                                     border_col = NA) {
  structure(
    list(
      type = "genes",
      height = height,
      features = features,
      fill = fill,
      label = label,
      label_cex = label_cex,
      label_col = label_col,
      label_min_width = label_min_width,
      border_col = border_col
    ),
    class = "vi_virus_circos_track"
  )
}

#' Create a virus read track specification
#'
#' @param height Track height.
#' @param bins Number of bins across the virus genome.
#' @param fill Histogram fill color.
#' @return A \code{vi_virus_circos_track} object.
#' @export
track_virus_reads <- function(height = 0.12, bins = 25, fill = "#A6CEE3") {
  bins <- validate_virus_circos_bins(bins)
  structure(
    list(
      type = "virus_reads",
      height = height,
      bins = bins,
      fill = fill,
      border_col = "white"
    ),
    class = "vi_virus_circos_track"
  )
}

#' Add a virus track to an existing virus-circos plot object
#'
#' @param x A \code{vi_virus_circos_plot} object.
#' @param track A \code{vi_virus_circos_track} object.
#' @return An updated \code{vi_virus_circos_plot} object.
#' @export
add_virus_circos_track <- function(x, track) {
  if (!inherits(x, "vi_virus_circos_plot")) {
    stop("x must be a vi_virus_circos_plot object.", call. = FALSE)
  }
  x$tracks[[length(x$tracks) + 1L]] <- validate_virus_circos_track(track)
  x$drawn <- FALSE
  x
}

#' Plot a virus genome with concentric circos tracks
#'
#' The plot initializes a single circular sector whose length equals the virus
#' genome length. Feature intervals are laid out from the outermost ring inward.
#' All intervals with the same feature name stay on one ring, which keeps split
#' records for a circular virus gene visually connected across the origin.
#' Overlapping feature groups are pushed to the next inner ring. Labels are
#' centered in the colored blocks when the interval is wide enough, including
#' across-origin centers for split circular features; short labels are omitted
#' and reported with a warning.
#'
#' Virus read counts are prepared from \code{integrations} before drawing and
#' then rendered as a single outer histogram track.
#'
#' @param virus_length Virus genome length in base pairs.
#' @param virus_name Virus sequence name used for the sector label.
#' @param features Optional virus feature intervals.
#' @param integrations Optional integration records used by virus read tracks.
#' @param tracks Optional list of \code{vi_virus_circos_track} objects.
#' @param fill Default feature fill color, unnamed palette, or named vector by
#'   feature name or feature type. When \code{NULL}, each distinct feature name
#'   gets its own color.
#' @param height Track height for the gene layer created from \code{features}.
#' @param label Logical. Whether to label gene blocks.
#' @param label_cex Gene label text size.
#' @param label_col Gene label color.
#' @param label_min_width Minimum feature width, in base pairs, required for a
#'   gene label to be drawn.
#' @param border_col Gene block border color. Defaults to \code{NA} so split
#'   circular features do not show a seam at the origin.
#' @param show_virus_name Logical. Whether to draw the virus name in the center
#'   of the circular plot.
#' @param track_margin Spacing between concentric tracks.
#' @param start.degree Starting angle for the circular layout.
#' @param gap.degree Gap between sectors. With one sector, keep this at zero.
#' @param clear Logical. Whether to clear the current circlize device first.
#' @param draw Logical. Whether to draw immediately.
#' @return A \code{vi_virus_circos_plot} object.
#' @export
plot_virus_circos <- function(virus_length, virus_name = "virus",
                              features = NULL, integrations = NULL,
                              tracks = NULL, fill = NULL,
                              height = 0.12, label = TRUE,
                              label_cex = 0.55,
                              label_col = "grey20",
                              label_min_width = 250,
                              border_col = NA,
                              show_virus_name = TRUE,
                              track_margin = c(0.004, 0.004),
                              start.degree = 90, gap.degree = 0,
                              clear = TRUE, draw = TRUE) {
  virus_length <- validate_virus_circos_length(virus_length)
  virus_name <- validate_virus_circos_name(virus_name)
  track_margin <- validate_virus_circos_margin(track_margin)
  if (!is.logical(show_virus_name) || length(show_virus_name) != 1L ||
      is.na(show_virus_name)) {
    stop("show_virus_name must be TRUE or FALSE.", call. = FALSE)
  }

  plot_tracks <- list()
  if (!is.null(tracks)) {
    if (!is.list(tracks)) {
      stop("tracks must be a list of vi_virus_circos_track objects.", call. = FALSE)
    }
    plot_tracks <- tracks
  }
  if (!is.null(features)) {
    plot_tracks[[length(plot_tracks) + 1L]] <- prepare_virus_circos_gene_track(
      virus_circos_track_genes(
        features = features,
        height = height,
        fill = fill,
        label = label,
        label_cex = label_cex,
        label_col = label_col,
        label_min_width = label_min_width,
        border_col = border_col
      ),
      virus_length = virus_length
    )
  }

  if (any(vapply(plot_tracks, function(track) identical(track$type, "virus_reads"),
                 logical(1)))) {
    if (is.null(integrations)) {
      stop("integrations must be supplied when using track_virus_reads().", call. = FALSE)
    }
    integrations <- as_integrations(integrations)
    validate_integration_virus_name(integrations, virus_name)
    plot_tracks <- lapply(
      plot_tracks,
      function(track) {
        if (identical(track$type, "virus_reads")) {
          prepare_virus_reads_track(
            track = track,
            integrations = integrations,
            virus_length = virus_length
          )
        } else {
          track
        }
      }
    )
  }

  x <- structure(
    list(
      virus_name = virus_name,
      virus_length = virus_length,
      track_margin = track_margin,
      start.degree = start.degree,
      gap.degree = gap.degree,
      show_virus_name = show_virus_name,
      tracks = lapply(plot_tracks, validate_virus_circos_track),
      drawn = FALSE
    ),
    class = "vi_virus_circos_plot"
  )

  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) {
    stop("draw must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(draw)) {
    x <- draw_virus_circos_plot(x, clear = clear)
    x$drawn <- TRUE
  }

  x
}

#' Draw a virus-circos plot object
#'
#' @param x A \code{vi_virus_circos_plot} object.
#' @param clear Logical. Whether to clear the current circlize device first.
#' @return The updated \code{vi_virus_circos_plot} object.
#' @export
draw_virus_circos_plot <- function(x, clear = TRUE) {
  if (!inherits(x, "vi_virus_circos_plot")) {
    stop("x must be a vi_virus_circos_plot object.", call. = FALSE)
  }

  if (isTRUE(clear)) {
    circlize::circos.clear()
  }

  circlize::circos.par(
    start.degree = x$start.degree,
    gap.degree = x$gap.degree,
    cell.padding = c(0, 0, 0, 0),
    track.margin = x$track_margin,
    points.overflow.warning = FALSE
  )

  circlize::circos.initialize(factors = x$virus_name,
                              xlim = c(0, x$virus_length))

  for (track in x$tracks) {
    draw_virus_circos_track(track, virus_name = x$virus_name,
                            virus_length = x$virus_length)
  }

  if (isTRUE(x$show_virus_name)) {
    draw_center_virus_name(x$virus_name)
  }

  x$drawn <- TRUE
  invisible(x)
}

#' Print a virus-circos plot
#'
#' @param x A \code{vi_virus_circos_plot} object.
#' @param ... Ignored.
#' @export
print.vi_virus_circos_plot <- function(x, ...) {
  if (isTRUE(x$drawn)) {
    invisible(x)
  } else {
    return(draw_virus_circos_plot(x))
  }
}

#' Draw a single virus-circos track
#'
#' @param track A \code{vi_virus_circos_track} object.
#' @param virus_name Virus sequence name.
#' @param virus_length Virus genome length in base pairs.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_virus_circos_track <- function(track, virus_name, virus_length) {
  track <- validate_virus_circos_track(track)
  if (track$type == "genes") {
    return(draw_virus_circos_gene_track(
      track = track,
      virus_name = virus_name,
      virus_length = virus_length
    ))
  }

  if (track$type == "virus_reads") {
    return(draw_virus_reads_track(
      track = track,
      virus_name = virus_name
    ))
  }

  stop("Unsupported virus track type: ", track$type, call. = FALSE)
}

validate_virus_circos_track <- function(track) {
  if (!inherits(track, "vi_virus_circos_track")) {
    stop("Each track must be created by virus_circos_track_genes() or track_virus_reads().",
         call. = FALSE)
  }

  if (!is.numeric(track$height) || length(track$height) != 1L || is.na(track$height) ||
      track$height <= 0) {
    stop("track height must be a single positive number.", call. = FALSE)
  }

  if (identical(track$type, "virus_reads")) {
    track$bins <- validate_virus_circos_bins(track$bins)
  }

  track
}

validate_virus_circos_length <- function(virus_length) {
  if (!is.numeric(virus_length) || length(virus_length) != 1L ||
      is.na(virus_length) || virus_length <= 0) {
    stop("virus_length must be a single positive number.", call. = FALSE)
  }
  as.numeric(virus_length)
}

validate_virus_circos_name <- function(virus_name) {
  if (!is.character(virus_name) || length(virus_name) != 1L ||
      is.na(virus_name) || !nzchar(virus_name)) {
    stop("virus_name must be a single non-empty string.", call. = FALSE)
  }
  trimws(gsub("(?i)^chr", "", virus_name))
}

validate_virus_circos_margin <- function(track_margin) {
  if (!is.numeric(track_margin) || length(track_margin) != 2L ||
      any(is.na(track_margin)) || any(track_margin < 0)) {
    stop("track_margin must contain two non-negative numbers.", call. = FALSE)
  }
  as.numeric(track_margin)
}

validate_virus_circos_bins <- function(bins) {
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) ||
      bins <= 0 || bins != as.integer(bins)) {
    stop("bins must be a single positive integer.", call. = FALSE)
  }
  as.integer(bins)
}

assign_virus_feature_levels <- function(features, virus_length) {
  virus_length <- validate_virus_circos_length(virus_length)
  df <- as.data.frame(
    normalize_virus_features(features, virus_length = virus_length),
    stringsAsFactors = FALSE
  )

  df$level <- 1L
  feature_order <- unique(df$feature)
  level_intervals <- list()
  for (feature_name in feature_order) {
    idx <- which(df$feature == feature_name)
    feature_intervals <- df[idx, c("start", "end"), drop = FALSE]
    placed <- FALSE
    for (level in seq_along(level_intervals)) {
      if (!virus_feature_group_overlaps(feature_intervals, level_intervals[[level]])) {
        df$level[idx] <- level
        level_intervals[[level]] <- rbind(level_intervals[[level]], feature_intervals)
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      level_intervals[[length(level_intervals) + 1L]] <- feature_intervals
      df$level[idx] <- length(level_intervals)
    }
  }

  df
}

virus_feature_group_overlaps <- function(feature_intervals, level_intervals) {
  if (is.null(level_intervals) || nrow(level_intervals) == 0L) {
    return(FALSE)
  }

  for (i in seq_len(nrow(feature_intervals))) {
    overlap <- feature_intervals$start[i] < level_intervals$end &
      feature_intervals$end[i] > level_intervals$start
    if (any(overlap, na.rm = TRUE)) {
      return(TRUE)
    }
  }

  FALSE
}

draw_virus_circos_gene_track <- function(track, virus_name, virus_length) {
  if (is.null(track$prepared_features)) {
    track <- prepare_virus_circos_gene_track(track, virus_length = virus_length)
  }
  feature_df <- track$prepared_features

  levels <- sort(unique(feature_df$level))
  for (level in levels) {
    level_df <- feature_df[feature_df$level == level, , drop = FALSE]
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = track$height,
      bg.border = NA,
      panel.fun = function(x, y) {
        chr <- circlize::CELL_META$sector.index
        if (!identical(chr, virus_name)) {
          return(NULL)
        }

        circlize::circos.rect(
          xleft = level_df$plot_start,
          ybottom = 0.15,
          xright = level_df$plot_end,
          ytop = 0.85,
          col = level_df$fill,
          border = track$border_col
        )

        if (isTRUE(track$label) && any(!is.na(level_df$draw_label))) {
          label_df <- level_df[!is.na(level_df$draw_label), , drop = FALSE]
          circlize::circos.text(
            x = label_df$label_x,
            y = 0.5,
            labels = label_df$draw_label,
            facing = "inside",
            niceFacing = TRUE,
            adj = c(0.5, 0.5),
            cex = track$label_cex,
            col = track$label_col
          )
        }
      }
    )
  }

  invisible(NULL)
}

prepare_virus_circos_gene_track <- function(track, virus_length) {
  track <- validate_virus_circos_track(track)
  feature_df <- assign_virus_feature_levels(track$features, virus_length = virus_length)
  feature_df <- set_virus_circos_plot_intervals(
    feature_df,
    virus_length = virus_length
  )
  feature_df$fill_key <- resolve_virus_circos_fill_key(feature_df, track$fill)
  fill_map <- resolve_virus_circos_fill(track$fill, unique(feature_df$fill_key))
  feature_df$fill <- unname(fill_map[feature_df$fill_key])
  feature_df$draw_label <- ifelse(
    feature_df$plot_end - feature_df$plot_start >= track$label_min_width,
    feature_df$feature,
    NA_character_
  )
  feature_df <- keep_one_virus_feature_label(feature_df, virus_length = virus_length)

  short_features <- setdiff(unique(feature_df$feature), stats::na.omit(feature_df$draw_label))
  if (length(short_features) > 0L && isTRUE(track$label)) {
    warning(
      "Labels omitted for short virus features: ",
      paste(short_features, collapse = ", "),
      call. = FALSE
    )
  }

  track$prepared_features <- feature_df
  track
}

set_virus_circos_plot_intervals <- function(feature_df, virus_length) {
  feature_df$plot_start <- pmax(0, pmin(feature_df$start, virus_length))
  feature_df$plot_end <- pmax(0, pmin(feature_df$end, virus_length))

  for (feature_name in unique(feature_df$feature)) {
    idx <- which(feature_df$feature == feature_name)
    if (length(idx) < 2L) {
      next
    }

    origin_idx <- idx[feature_df$start[idx] <= 1]
    terminal_candidates <- idx[feature_df$start[idx] > 1]
    if (length(origin_idx) == 0L || length(terminal_candidates) == 0L) {
      next
    }

    terminal_idx <- terminal_candidates[
      feature_df$end[terminal_candidates] == max(feature_df$end[terminal_candidates],
                                                 na.rm = TRUE)
    ]
    feature_df$plot_start[origin_idx] <- 0
    feature_df$plot_end[terminal_idx] <- virus_length
  }

  feature_df
}

keep_one_virus_feature_label <- function(feature_df, virus_length) {
  feature_df$label_x <- (feature_df$plot_start + feature_df$plot_end) / 2
  for (feature_name in unique(feature_df$feature)) {
    idx <- which(feature_df$feature == feature_name & !is.na(feature_df$draw_label))
    if (length(idx) == 0L) {
      next
    }

    wrapped_label <- circular_virus_feature_label_position(
      feature_df[idx, , drop = FALSE],
      virus_length = virus_length
    )
    if (!is.null(wrapped_label)) {
      keep <- wrapped_label$row
      feature_df$label_x[keep] <- wrapped_label$x
    } else {
      widths <- feature_df$plot_end[idx] - feature_df$plot_start[idx]
      keep <- idx[which.max(widths)]
    }

    feature_df$draw_label[setdiff(idx, keep)] <- NA_character_
  }

  feature_df
}

circular_virus_feature_label_position <- function(feature_df, virus_length) {
  origin_idx <- which(feature_df$plot_start <= 0)
  terminal_idx <- which(feature_df$plot_end >= virus_length)
  if (length(origin_idx) == 0L || length(terminal_idx) == 0L) {
    return(NULL)
  }

  terminal_start <- min(feature_df$plot_start[terminal_idx], na.rm = TRUE)
  origin_end <- max(feature_df$plot_end[origin_idx], na.rm = TRUE)
  width <- (virus_length - terminal_start) + origin_end
  if (!is.finite(width) || width <= 0) {
    return(NULL)
  }

  label_x <- (terminal_start + width / 2) %% virus_length
  contains_label <- feature_df$plot_start <= label_x & feature_df$plot_end >= label_x
  candidate_rows <- which(contains_label)
  if (length(candidate_rows) == 0L) {
    widths <- feature_df$plot_end - feature_df$plot_start
    row <- which.max(widths)
  } else {
    row <- candidate_rows[1L]
  }

  list(row = as.integer(rownames(feature_df)[row]), x = label_x)
}

draw_center_virus_name <- function(virus_name) {
  graphics::text(
    x = 0,
    y = 0,
    labels = virus_name,
    cex = 1,
    col = "grey20"
  )
  invisible(NULL)
}

prepare_virus_reads_track <- function(track, integrations, virus_length) {
  track <- validate_virus_circos_track(track)
  data <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)

  pos <- suppressWarnings(as.numeric(data$virus_pos))
  keep <- !is.na(pos) & pos > 0 & pos <= virus_length
  if (!all(keep)) {
    warning(sum(!keep), " virus read records were dropped because their positions were invalid.",
            call. = FALSE)
  }
  pos <- pos[keep]
  if (length(pos) == 0L) {
    track$data <- data.frame(start = numeric(), end = numeric(), value = numeric())
    return(track)
  }

  weights <- suppressWarnings(as.numeric(data$support_reads))[keep]
  weights[is.na(weights)] <- 0

  breaks <- seq(0, virus_length, length.out = as.integer(track$bins) + 1L)
  bins <- findInterval(pos, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  counts <- tapply(weights, bins, sum)
  counts <- counts[as.character(seq_len(length(breaks) - 1L))]
  counts <- as.numeric(counts)
  counts[is.na(counts)] <- 0

  track$data <- data.frame(
    start = breaks[-length(breaks)],
    end = breaks[-1],
    value = counts,
    stringsAsFactors = FALSE
  )

  track
}

draw_virus_reads_track <- function(track, virus_name) {
  data <- track$data
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0L) {
    return(invisible(NULL))
  }

  max_value <- max(data$value, na.rm = TRUE)
  if (!is.finite(max_value) || max_value <= 0) {
    max_value <- 1
  }

  circlize::circos.trackPlotRegion(
    ylim = c(0, max_value),
    track.height = track$height,
    bg.border = NA,
    panel.fun = function(x, y) {
      if (!identical(circlize::CELL_META$sector.index, virus_name)) {
        return(NULL)
      }

      circlize::circos.rect(
        xleft = data$start,
        ybottom = 0,
        xright = data$end,
        ytop = data$value,
        col = track$fill,
        border = track$border_col
      )
    }
  )

  invisible(NULL)
}

resolve_virus_circos_fill_key <- function(feature_df, fill) {
  if (!is.null(fill) && !is.null(names(fill))) {
    type_names <- unique(feature_df$type)
    feature_names <- unique(feature_df$feature)
    if (all(names(fill) %in% feature_names)) {
      return(feature_df$feature)
    }
    if (all(names(fill) %in% type_names)) {
      return(feature_df$type)
    }
    if (any(names(fill) %in% feature_names)) {
      return(feature_df$feature)
    }
    if (any(names(fill) %in% type_names)) {
      return(feature_df$type)
    }
  }

  feature_df$feature
}

resolve_virus_circos_fill <- function(fill, keys) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0L) {
    keys <- "feature"
  }

  if (is.null(fill)) {
    cols <- grDevices::hcl.colors(length(keys), "Set 3")
    return(stats::setNames(cols, keys))
  }

  if (!is.null(names(fill))) {
    mapped <- fill[keys]
    missing <- is.na(mapped)
    if (any(missing)) {
      mapped[missing] <- rep(unname(fill)[1L], sum(missing))
    }
    return(stats::setNames(unname(mapped), keys))
  }

  if (length(fill) == 1L) {
    return(stats::setNames(rep(fill, length(keys)), keys))
  }

  stats::setNames(rep(fill, length.out = length(keys)), keys)
}
