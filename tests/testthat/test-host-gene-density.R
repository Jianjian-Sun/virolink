test_that("host gene density is prepared from a local GTF annotation", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(
    "1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\"; gene_type \"protein_coding\";",
    "1\ttest\tgene\t800000\t800100\t.\t-\t.\tgene_id \"g2\"; gene_name \"B\"; gene_type \"lncRNA\";",
    "1\ttest\tgene\t1200000\t1200100\t.\t+\t.\tgene_id \"g3\"; gene_name \"C\"; gene_type \"protein_coding\";"
  ), gtf)

  host <- host_genome(
    sequences = data.frame(chr = "1", start = 0, end = 2000000),
    annotation = gtf
  )
  virus <- virus_genome(name = "HBV", length = 1000)
  integrations <- data.frame(
    host_chr = "1", host_pos = 100, host_strand = "*",
    virus_chr = "HBV", virus_pos = 10, virus_strand = "+", support_reads = 5
  )

  p <- plot_integrations(
    integrations = integrations,
    host = host,
    virus = virus,
    tracks = list(track_host_gene_density()),
    draw = FALSE
  )

  density_track <- p$tracks[[1]]
  expect_identical(density_track$type, "host_gene_density")
  expect_equal(density_track$data$count, c(2, 1))
  expect_identical(density_track$data$chr, c("1", "1"))
  expect_true(inherits(p$host$genes, "data.frame"))
})

test_that("host gene density draws only on host sectors", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines("1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id \"g1\";", gtf)
  host <- host_genome(data.frame(chr = "1", start = 0, end = 1000000), annotation = gtf)
  virus <- virus_genome(name = "HBV", length = 1000)
  integrations <- data.frame(host_chr = "1", host_pos = 100, host_strand = "*",
                             virus_chr = "HBV", virus_pos = 10, virus_strand = "+", support_reads = 5)
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 600, height = 600)
  on.exit({ circlize::circos.clear(); grDevices::dev.off() }, add = TRUE)

  expect_no_error(plot_integrations(
    integrations, host, virus,
    tracks = list(track_ideogram(), track_host_gene_density(), track_sites(), track_links())
  ))
})
