test_that("breakpoint map plots host and virus coordinates", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2"),
    host_pos = c(100, 500, 200),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30),
    support_reads = c(5, 10, 8),
    sample = c("A", "B", "A"),
    virus_strand = c("+", "-", "+"),
    library_type = c("DNA", "DNA", "RNA")
  )
  host <- data.frame(
    chr = c("1", "2"),
    start = c(0, 0),
    end = c(1000, 800)
  )
  features <- virus_features(
    feature = c("A", "B"),
    start = c(1, 25),
    end = c(15, 40),
    type = c("gene", "gene")
  )

  p <- plot_breakpoint_map(
    integrations,
    host = host,
    virus = c(HBV = 50),
    features = features,
    group_by = "library_type",
    size_by = "support_reads"
  )

  expect_s3_class(p, "ggplot")
  build <- ggplot2::ggplot_build(p)
  expect_gt(length(unique(build$data[[3]]$colour)), 1L)
  expect_gt(length(unique(build$data[[3]]$size)), 1L)
  expect_lte(max(build$data[[3]]$y), 50)
})

test_that("breakpoint map accepts virus genome objects", {
  integrations <- data.frame(
    host_chr = "1",
    host_pos = 100,
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = 10,
    support_reads = 5,
    sample = "A",
    virus_strand = "+",
    library_type = "DNA"
  )
  host <- data.frame(chr = "1", start = 0, end = 1000)
  virus <- virus_genome("HBV", 50)

  p <- plot_breakpoint_map(integrations, host = host, virus = virus)

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("breakpoint map validates mapped columns", {
  integrations <- data.frame(
    host_chr = "1",
    host_pos = 100,
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = 10,
    support_reads = 5,
    sample = "A",
    virus_strand = "+",
    library_type = "DNA"
  )
  host <- data.frame(chr = "1", start = 0, end = 1000)

  expect_error(
    plot_breakpoint_map(integrations, host = host, virus = c(HBV = 50), group_by = "missing"),
    "Unknown group_by column"
  )
  expect_error(
    plot_breakpoint_map(integrations, host = host, virus = c(HBV = 50), size_by = "library_type"),
    "size_by must refer to a numeric column"
  )
})

test_that("breakpoint map drops out-of-range coordinates", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2"),
    host_pos = c(100, 1500, 200),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 80),
    support_reads = c(5, 10, 8),
    sample = c("A", "B", "C"),
    virus_strand = c("+", "-", "+"),
    library_type = c("DNA", "DNA", "RNA")
  )
  host <- data.frame(chr = c("1", "2"), start = c(0, 0), end = c(1000, 800))

  expect_warning(
    p <- plot_breakpoint_map(integrations, host = host, virus = c(HBV = 50)),
    "outside the host or virus genome"
  )
  build <- ggplot2::ggplot_build(p)
  expect_equal(nrow(build$data[[2]]), 1L)
})
