test_that("interactive explorer returns an HTML component", {
  testthat::skip_if_not_installed("ggiraph")
  testthat::skip_if_not_installed("DT")
  testthat::skip_if_not_installed("htmltools")

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
  host <- data.frame(chr = c("1", "2"), start = c(0, 0), end = c(1000, 800))

  explorer <- plot_interactive_explorer(
    integrations,
    host = host,
    virus = c(HBV = 50),
    group_by = "library_type",
    size_by = "support_reads"
  )

  expect_s3_class(explorer, "shiny.tag.list")
  expect_equal(length(explorer), 2L)
})

test_that("interactive explorer supports selected table columns", {
  testthat::skip_if_not_installed("ggiraph")
  testthat::skip_if_not_installed("DT")
  testthat::skip_if_not_installed("htmltools")

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

  explorer <- plot_interactive_explorer(
    integrations,
    host = host,
    virus = c(HBV = 50),
    table_columns = c("sample", "host_chr", "host_pos"),
    group_by = NULL,
    size_by = NULL
  )

  expect_s3_class(explorer, "shiny.tag.list")
})

test_that("interactive explorer validates inputs", {
  testthat::skip_if_not_installed("ggiraph")
  testthat::skip_if_not_installed("DT")
  testthat::skip_if_not_installed("htmltools")

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
    plot_interactive_explorer(integrations, host = host, virus = c(HBV = 50), group_by = "missing"),
    "Unknown group_by column"
  )
  expect_error(
    plot_interactive_explorer(integrations, host = host, virus = c(HBV = 50), size_by = "library_type"),
    "size_by must refer to a numeric column"
  )
  expect_error(
    plot_interactive_explorer(
      integrations,
      host = host,
      virus = c(HBV = 50),
      table_columns = "missing"
    ),
    "Unknown table_columns"
  )
})
