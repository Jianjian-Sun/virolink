test_that("integration signature plots candidate clonality clusters", {
  integrations <- data.frame(
    host_chr = c("1", "1", "1", "2"),
    host_pos = c(1000, 1100, 10000, 500),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 20, 30, 40),
    support_reads = c(5, 10, 8, 2),
    sample = c("A", "B", "A", "B"),
    virus_strand = c("+", "-", "+", "-"),
    library_type = c("DNA", "RNA", "DNA", "RNA")
  )

  p <- plot_integration_signature(
    integrations,
    mode = "clonality",
    group_by = "sample",
    cluster_window = 200,
    measure = "events"
  )

  expect_s3_class(p, "ggplot")
  build <- ggplot2::ggplot_build(p)
  expect_equal(sum(build$data[[1]]$ymax - build$data[[1]]$ymin), 4)
  expect_gt(length(unique(build$data[[1]]$fill)), 1L)
})

test_that("host-virus clustering can split host-near virus-distant records", {
  integrations <- data.frame(
    host_chr = c("1", "1"),
    host_pos = c(1000, 1100),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(10, 5000),
    support_reads = c(5, 10),
    sample = c("A", "A"),
    virus_strand = c("+", "-"),
    library_type = c("DNA", "RNA")
  )

  host_only <- plot_integration_signature(
    integrations,
    mode = "clonality",
    group_by = NULL,
    cluster_by = "host",
    cluster_window = 200
  )
  host_virus <- plot_integration_signature(
    integrations,
    mode = "clonality",
    group_by = NULL,
    cluster_by = "host_virus",
    cluster_window = 200
  )

  expect_equal(nrow(ggplot2::ggplot_build(host_only)$data[[1]]), 1L)
  expect_equal(nrow(ggplot2::ggplot_build(host_virus)$data[[1]]), 2L)
})

test_that("integration signature plots strand orientation", {
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

  p <- plot_integration_signature(
    integrations,
    mode = "strand",
    group_by = "sample",
    measure = "support_reads"
  )

  expect_s3_class(p, "ggplot")
  build <- ggplot2::ggplot_build(p)
  expect_equal(sum(build$data[[1]]$y), 23)
})

test_that("integration signature validates inputs", {
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
    plot_integration_signature(integrations, group_by = "missing"),
    "Unknown group_by column"
  )
  expect_error(
    plot_integration_signature(integrations, cluster_window = 0),
    "cluster_window must be a single positive number"
  )
  expect_error(
    plot_integration_signature(integrations, top_n = 1.5),
    "top_n must be a single positive integer"
  )
})
