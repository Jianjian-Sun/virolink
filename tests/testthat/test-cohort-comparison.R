test_that("cohort comparison plots event burden by group and fill", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2", "2"),
    host_pos = c(100, 200, 300, 400),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30, 40),
    support_reads = c(5, 10, 8, 2),
    sample = c("A", "A", "B", "B"),
    virus_strand = c("+", "-", "+", "-"),
    library_type = c("DNA", "RNA", "DNA", "RNA")
  )

  p <- plot_cohort_comparison(integrations, group_by = "sample", fill_by = "library_type")

  expect_s3_class(p, "ggplot")
  build <- ggplot2::ggplot_build(p)
  expect_equal(sum(build$data[[1]]$ymax - build$data[[1]]$ymin), 4)
  expect_gt(length(unique(build$data[[1]]$fill)), 1L)
})

test_that("cohort comparison summarizes support and normalizes proportions", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2", "2"),
    host_pos = c(100, 200, 300, 400),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30, 40),
    support_reads = c(5, 15, 8, 2),
    sample = c("A", "A", "B", "B"),
    virus_strand = c("+", "-", "+", "-"),
    library_type = c("DNA", "RNA", "DNA", "RNA")
  )

  p <- plot_cohort_comparison(
    integrations,
    group_by = "sample",
    fill_by = "library_type",
    measure = "support_reads",
    normalize = TRUE
  )
  build <- ggplot2::ggplot_build(p)
  group_totals <- stats::aggregate(
    build$data[[1]]$ymax - build$data[[1]]$ymin,
    by = list(x = build$data[[1]]$x),
    FUN = sum
  )

  expect_equal(group_totals$x, c(1, 2))
  expect_equal(round(group_totals[[2]], 8), c(1, 1))
})

test_that("cohort comparison works without fill_by", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2"),
    host_pos = c(100, 200, 300),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30),
    support_reads = c(5, 10, 8),
    sample = c("A", "A", "B"),
    virus_strand = c("+", "-", "+"),
    library_type = c("DNA", "RNA", "DNA")
  )

  p <- plot_cohort_comparison(integrations, fill_by = NULL)

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("cohort comparison works without sample when grouping is omitted", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2"),
    host_pos = c(100, 200, 300),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30),
    support_reads = c(5, 10, 8),
    virus_strand = c("+", "-", "+"),
    library_type = c("DNA", "RNA", "DNA")
  )

  p <- plot_cohort_comparison(integrations)

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("cohort comparison collapses low-frequency fill categories", {
  integrations <- data.frame(
    host_chr = rep("1", 5),
    host_pos = seq(100, 500, by = 100),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = seq(10, 50, by = 10),
    support_reads = rep(1, 5),
    sample = rep("A", 5),
    virus_strand = rep("+", 5),
    library_type = c("A", "A", "B", "C", "D")
  )

  p <- plot_cohort_comparison(integrations, fill_by = "library_type", top_n = 2)
  build <- ggplot2::ggplot_build(p)

  expect_equal(nrow(build$data[[1]]), 3L)
})

test_that("cohort comparison validates inputs", {
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

  expect_error(
    plot_cohort_comparison(integrations, group_by = "missing"),
    "Unknown group_by column"
  )
  expect_error(
    plot_cohort_comparison(integrations, fill_by = NULL, top_n = 2),
    "top_n can only be used"
  )
  expect_error(
    plot_cohort_comparison(integrations, normalize = NA),
    "normalize must be TRUE or FALSE"
  )
})
