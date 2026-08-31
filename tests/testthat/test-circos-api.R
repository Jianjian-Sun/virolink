test_that("plot_integrations can build first and draw once", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2", "2"),
    host_pos = c(100, 600, 300, 800),
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = c(50, 150, 450, 850),
    sample = c("A", "B", "A", "B"),
    support_reads = c(10, 20, 15, 30),
    virus_strand = c("+", "-", "+", "-"),
    library_type = c("single_end", "paired_end", "single_end", "paired_end")
  )
  host <- host_genome(data.frame(
    chr = c("1", "2"),
    start = c(0, 0),
    end = c(1000, 1000)
  ))
  virus <- virus_genome(c(HBV = 1000))
  features <- virus_features(
    feature = c("geneA", "short"),
    start = c(50, 500),
    end = c(450, 560),
    type = c("gene", "gene")
  )

  p <- plot_integrations(
    integrations,
    host = host,
    virus = virus,
    tracks = list(
      track_ideogram(),
      track_sites(color = "library_type", label = "Sites"),
      track_virus_genes(features, label_min_width = 100),
      track_virus_density(label = "Virus density"),
      track_links(color = "library_type")
    ),
    draw = FALSE
  )

  expect_s3_class(p, "vi_integration_plot")
  expect_false(p$drawn)
  expect_true("library_type" %in% colnames(p$plot_df))
  expect_named(
    integration_group_legend_spec(p$tracks, p$plot_df),
    c("paired_end", "single_end")
  )

  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 900, height = 700)
  on.exit({
    circlize::circos.clear()
    grDevices::dev.off()
  }, add = TRUE)

  expect_no_error(p <- draw_integration_plot(p))
  expect_true(p$drawn)
  expect_no_error(print(p))
})

test_that("plot_integrations rejects a virus name mismatch", {
  integrations <- data.frame(
    host_chr = "1",
    host_pos = 100,
    host_strand = "*",
    virus_chr = "HPV16",
    virus_pos = 10,
    virus_strand = "+",
    support_reads = 5
  )
  host <- host_genome(data.frame(chr = "1", start = 0, end = 1000))
  virus <- virus_genome(c(HBV = 1000))

  expect_error(
    plot_integrations(integrations, host = host, virus = virus, draw = FALSE),
    "virus_chr.*HBV"
  )
})

test_that("create_gi_from_table keeps strands and explicit metadata", {
  integrations <- data.frame(
    host_chr = "1",
    host_pos = 100,
    host_strand = "+",
    virus_chr = "HBV",
    virus_pos = 10,
    virus_strand = "-",
    support_reads = 5,
    sample = "A",
    library_type = "single_end"
  )
  file <- tempfile(fileext = ".tsv")
  utils::write.table(integrations, file, sep = "\t", row.names = FALSE, quote = FALSE)
  cfg <- create_config(
    host = data.frame(chr = "1", start = 0, end = 1000),
    virus_name = "HBV",
    virus_length = 50
  )

  gi <- create_gi_from_table(file, cfg)

  expect_s4_class(gi, "GInteractions")
  anchors <- InteractionSet::anchors(gi)
  expect_equal(as.character(GenomicRanges::strand(anchors$first)), "+")
  expect_equal(as.character(GenomicRanges::strand(anchors$second)), "-")
  expect_equal(S4Vectors::mcols(gi)$sample, "A")
  expect_equal(S4Vectors::mcols(gi)$library_type, "single_end")
})

test_that("create_gi_from_table treats unknown host strand as plus for plotting", {
  integrations <- data.frame(
    host_chr = "1",
    host_pos = 100,
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = 10,
    virus_strand = "*",
    support_reads = 5
  )
  file <- tempfile(fileext = ".tsv")
  utils::write.table(integrations, file, sep = "\t", row.names = FALSE, quote = FALSE)
  cfg <- create_config(
    host = data.frame(chr = "1", start = 0, end = 1000),
    virus_name = "HBV",
    virus_length = 50
  )

  gi <- create_gi_from_table(file, cfg)

  anchors <- InteractionSet::anchors(gi)
  expect_equal(as.character(GenomicRanges::strand(anchors$first)), "+")
  expect_equal(as.character(GenomicRanges::strand(anchors$second)), "*")
})

test_that("create_gi_from_table rejects a virus name mismatch", {
  integrations <- data.frame(
    host_chr = "1",
    host_pos = 100,
    host_strand = "*",
    virus_chr = "HPV16",
    virus_pos = 10,
    virus_strand = "*",
    support_reads = 5
  )
  file <- tempfile(fileext = ".tsv")
  utils::write.table(integrations, file, sep = "\t", row.names = FALSE, quote = FALSE)
  cfg <- create_config(
    host = data.frame(chr = "1", start = 0, end = 1000),
    virus_name = "HBV",
    virus_length = 50
  )

  expect_error(
    create_gi_from_table(file, cfg),
    "virus_chr.*HBV"
  )
})
