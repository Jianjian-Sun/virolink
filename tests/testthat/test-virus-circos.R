test_that("assign_virus_feature_levels pushes overlaps inward", {
  features <- data.frame(
    feature = c("gene1", "gene2", "gene3"),
    start = c(10, 40, 120),
    end = c(80, 100, 180),
    stringsAsFactors = FALSE
  )

  result <- assign_virus_feature_levels(features, virus_length = 200)

  expect_equal(result$level, c(1L, 2L, 1L))
})

test_that("split records for the same virus gene stay on one track", {
  features <- data.frame(
    feature = c("S", "P", "X", "C", "P", "S"),
    start = c(1, 1, 1376, 1816, 2309, 2850),
    end = c(837, 1625, 1840, 2454, 3182, 3182),
    type = "gene",
    stringsAsFactors = FALSE
  )

  result <- assign_virus_feature_levels(features, virus_length = 3215)

  expect_length(unique(result$level[result$feature == "S"]), 1L)
  expect_length(unique(result$level[result$feature == "P"]), 1L)
  expect_false(identical(
    unique(result$level[result$feature == "S"]),
    unique(result$level[result$feature == "P"])
  ))
})

test_that("split records for the same circular virus gene close the origin", {
  features <- data.frame(
    feature = c("S", "P", "X", "C", "P", "S"),
    start = c(1, 1, 1376, 1816, 2309, 2850),
    end = c(837, 1625, 1840, 2454, 3182, 3182),
    type = "gene",
    stringsAsFactors = FALSE
  )

  track <- prepare_virus_circos_gene_track(
    virus_circos_track_genes(features),
    virus_length = 3215
  )
  prepared <- track$prepared_features

  s_features <- prepared[prepared$feature == "S", , drop = FALSE]
  p_features <- prepared[prepared$feature == "P", , drop = FALSE]

  expect_true(any(s_features$plot_start == 0))
  expect_true(any(s_features$plot_end == 3215))
  expect_true(any(p_features$plot_start == 0))
  expect_true(any(p_features$plot_end == 3215))
})

test_that("plot_virus_circos consumes circular features normalized by virus_features helpers", {
  features <- data.frame(
    feature = "P",
    start = 2309,
    end = 4807,
    type = "polymerase",
    stringsAsFactors = FALSE
  )

  track <- prepare_virus_circos_gene_track(
    virus_circos_track_genes(features),
    virus_length = 3182
  )

  prepared <- track$prepared_features
  expect_equal(nrow(prepared), 2L)
  expect_true(any(prepared$start == 1 & prepared$end == 1625))
  expect_true(any(prepared$plot_start == 0))
  expect_true(any(prepared$plot_end == 3182))
})

test_that("split records for the same virus gene only keep one label", {
  features <- data.frame(
    feature = c("S", "P", "X", "C", "P", "S"),
    start = c(1, 1, 1376, 1816, 2309, 2850),
    end = c(837, 1625, 1840, 2454, 3182, 3182),
    type = "gene",
    stringsAsFactors = FALSE
  )

  track <- prepare_virus_circos_gene_track(
    virus_circos_track_genes(features, label_min_width = 120),
    virus_length = 3215
  )
  prepared <- track$prepared_features

  expect_equal(sum(!is.na(prepared$draw_label[prepared$feature == "S"])), 1L)
  expect_equal(sum(!is.na(prepared$draw_label[prepared$feature == "P"])), 1L)
})

test_that("split virus gene labels are centered across the circular origin", {
  features <- data.frame(
    feature = c("S", "S"),
    start = c(1, 2850),
    end = c(837, 3182),
    type = "gene",
    stringsAsFactors = FALSE
  )

  track <- prepare_virus_circos_gene_track(
    virus_circos_track_genes(features, label_min_width = 120),
    virus_length = 3215
  )
  label_row <- track$prepared_features[
    !is.na(track$prepared_features$draw_label),
    ,
    drop = FALSE
  ]

  expected_center <- (2850 + ((3215 - 2850) + 837) / 2) %% 3215
  expect_equal(nrow(label_row), 1L)
  expect_equal(label_row$label_x, expected_center)
})

test_that("origin closure is internal rather than a public parameter", {
  expect_false("close_origin" %in% names(formals(plot_virus_circos)))
  expect_false("close_origin" %in% names(formals(virus_circos_track_genes)))
})

test_that("default virus gene colors differ by feature name", {
  features <- data.frame(
    feature = c("S", "P", "X"),
    start = c(1, 900, 1800),
    end = c(500, 1400, 2300),
    type = c("surface", "polymerase", "regulatory"),
    stringsAsFactors = FALSE
  )

  track <- prepare_virus_circos_gene_track(
    virus_circos_track_genes(features),
    virus_length = 3215
  )

  expect_equal(track$prepared_features$fill_key, features$feature)
  expect_equal(length(unique(track$prepared_features$fill)), 3L)
})

test_that("plot_virus_circos uses readable compact defaults", {
  features <- data.frame(
    feature = c("S", "P"),
    start = c(1, 900),
    end = c(500, 1600),
    stringsAsFactors = FALSE
  )

  p <- plot_virus_circos(
    virus_length = 3215,
    features = features,
    draw = FALSE
  )

  expect_equal(p$track_margin, c(0.004, 0.004))
  expect_equal(p$tracks[[1]]$label_cex, 0.55)
  expect_true(is.na(p$tracks[[1]]$border_col))
  expect_true(p$show_virus_name)
})

test_that("virus read tracks use integrations passed to plot_virus_circos", {
  features <- data.frame(
    feature = "S",
    start = 1,
    end = 500,
    stringsAsFactors = FALSE
  )
  integrations <- read_integrations(
    system.file("extdata", "data.txt", package = "virolink")
  )

  expect_false("data" %in% names(formals(track_virus_reads)))
  expect_false("label" %in% names(formals(track_virus_reads)))

  p <- plot_virus_circos(
    virus_length = 3182,
    virus_name = "HBV",
    features = features,
    integrations = integrations,
    tracks = list(track_virus_reads()),
    draw = FALSE
  )

  expect_equal(p$tracks[[1]]$type, "virus_reads")
  expect_equal(p$tracks[[1]]$bins, 25L)
  expect_equal(p$tracks[[1]]$fill, "#A6CEE3")
  expect_gt(nrow(p$tracks[[1]]$data), 0L)
  expect_equal(p$tracks[[2]]$type, "genes")
})

test_that("virus read tracks require integrations at plot time", {
  features <- data.frame(
    feature = "S",
    start = 1,
    end = 500,
    stringsAsFactors = FALSE
  )

  expect_error(
    plot_virus_circos(
      virus_length = 1000,
      features = features,
      tracks = list(track_virus_reads()),
      draw = FALSE
    ),
    "integrations"
  )
})

test_that("plot_virus_circos can omit the center virus name", {
  p <- plot_virus_circos(
    virus_length = 3215,
    virus_name = "HBV",
    show_virus_name = FALSE,
    draw = FALSE
  )

  expect_false(p$show_virus_name)
})

test_that("plot_virus_circos returns a plot object and warns on short labels", {
  features <- data.frame(
    feature = c("long_gene", "short_gene"),
    start = c(10, 150),
    end = c(120, 165),
    type = c("gene", "gene"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    p <- plot_virus_circos(
      virus_length = 200,
      features = features,
      draw = FALSE,
      label_min_width = 40
    ),
    "short_gene"
  )

  expect_s3_class(p, "vi_virus_circos_plot")
  expect_false(p$drawn)
})

test_that("plot_virus_circos validates virus_length", {
  expect_error(
    plot_virus_circos(virus_length = 0, draw = FALSE),
    "positive"
  )
})

test_that("virus read tracks validate bins before drawing", {
  expect_error(
    track_virus_reads(bins = 0),
    "bins must be a single positive integer"
  )
  expect_error(
    track_virus_reads(bins = NA),
    "bins must be a single positive integer"
  )
})
