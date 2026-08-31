test_that("landscape plots sites and density", {
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

  sites <- plot_integration_landscape(integrations, host = host)
  density <- plot_integration_landscape(
    integrations,
    host = host,
    mode = "density",
    group_by = "sample",
    bin_size = 250
  )

  expect_s3_class(sites, "ggplot")
  expect_s3_class(density, "ggplot")
  expect_no_error(ggplot2::ggplot_build(sites))

  density_build <- ggplot2::ggplot_build(density)
  expect_gt(length(unique(density_build$data[[2]]$fill)), 1L)
  expect_lte(max(density_build$data[[2]]$xmax), sum(host$end - host$start))
})

test_that("landscape validates grouping columns", {
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
    plot_integration_landscape(integrations, host = host, group_by = "missing"),
    "Unknown group_by column"
  )
})
