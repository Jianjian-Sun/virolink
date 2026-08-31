test_that("host features standardize common columns", {
  features <- host_features(data.frame(
    chrom = c("chr2", "2"),
    from = c(100, 300),
    to = c(200, 450),
    name = c("A", "B"),
    feature_type = c("gene", "enhancer"),
    stringsAsFactors = FALSE
  ))

  expect_s3_class(features, "vi_host_features")
  expect_equal(features$chr, c("2", "2"))
  expect_equal(colnames(features), c("chr", "start", "end", "feature", "type", "strand"))
})

test_that("locus context plots integrations and annotations", {
  integrations <- data.frame(
    host_chr = c("chr2", "chr2", "chr2", "chr3"),
    host_pos = c(1000, 1300, 9000, 1100),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30, 40),
    support_reads = c(5, 10, 8, 2),
    sample = c("A", "B", "C", "D"),
    virus_strand = c("+", "-", "+", "-"),
    library_type = c("DNA", "RNA", "DNA", "RNA")
  )
  annotations <- host_features(
    chr = c("2", "2"),
    start = c(800, 1250),
    end = c(1600, 2200),
    feature = c("geneA", "enhancerB"),
    type = c("gene", "regulatory")
  )

  p <- plot_locus_context(
    integrations,
    chr = "2",
    pos = 1200,
    window = 800,
    annotations = annotations,
    group_by = "library_type",
    size_by = "support_reads"
  )

  expect_s3_class(p, "ggplot")
  build <- ggplot2::ggplot_build(p)
  point_layers <- Filter(
    function(layer) "shape" %in% names(layer) && "colour" %in% names(layer),
    build$data
  )
  point_layer <- point_layers[[length(point_layers)]]
  expect_equal(nrow(point_layer), 2L)
  expect_gt(length(unique(point_layer$colour)), 1L)
  expect_gt(length(unique(point_layer$size)), 1L)
  expect_true(any(vapply(
    p$scales$scales,
    function(scale) "fill" %in% scale$aesthetics,
    logical(1)
  )))
  expect_true(any(vapply(
    build$data,
    function(layer) "fill" %in% names(layer) && any(layer$alpha == 0.45),
    logical(1)
  )))
})

test_that("locus context can hide the host feature type legend", {
  integrations <- data.frame(
    host_chr = c("chr2", "chr2"),
    host_pos = c(1000, 1300),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20),
    support_reads = c(5, 10),
    sample = c("A", "B"),
    virus_strand = c("+", "-"),
    library_type = c("DNA", "RNA")
  )
  annotations <- host_features(
    chr = "2",
    start = 800,
    end = 1600,
    feature = "geneA",
    type = "gene"
  )

  p <- plot_locus_context(
    integrations,
    chr = "2",
    pos = 1200,
    window = 800,
    annotations = annotations,
    size_by = NULL,
    show_feature_legend = FALSE
  )

  fill_scales <- Filter(
    function(scale) "fill" %in% scale$aesthetics,
    p$scales$scales
  )
  expect_length(fill_scales, 1L)
  expect_identical(fill_scales[[1]]$guide, "none")
})

test_that("locus context works without annotations", {
  integrations <- data.frame(
    host_chr = c("2", "2"),
    host_pos = c(1000, 1300),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20),
    support_reads = c(5, 10),
    sample = c("A", "B"),
    virus_strand = c("+", "-"),
    library_type = c("DNA", "RNA")
  )

  p <- plot_locus_context(
    integrations,
    chr = "chr2",
    pos = 1200,
    window = 500,
    size_by = NULL
  )

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("locus context validates inputs", {
  integrations <- data.frame(
    host_chr = "2",
    host_pos = 1000,
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = 10,
    support_reads = 5,
    sample = "A",
    virus_strand = "+",
    library_type = "DNA"
  )

  expect_error(
    plot_locus_context(integrations, chr = "2", pos = 1000, window = -1),
    "window must be a single positive number"
  )
  expect_error(
    plot_locus_context(integrations, chr = "2", pos = 1000, group_by = "missing"),
    "Unknown group_by column"
  )
  expect_error(
    plot_locus_context(integrations, chr = "2", pos = 1000, size_by = "library_type"),
    "size_by must refer to a numeric column"
  )
  expect_error(
    plot_locus_context(integrations, chr = "2", pos = 999999),
    "No integration records were found"
  )
  expect_error(
    plot_locus_context(integrations, chr = "2", pos = 1000, show_feature_legend = NA),
    "show_feature_legend must be TRUE or FALSE"
  )
})
