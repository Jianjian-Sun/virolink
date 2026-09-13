test_that("host_genome stores canonical sequences and automatic colors", {
  host <- host_genome(data.frame(
    chr = c("chr2", "chr1"),
    start = c(0, 0),
    end = c(200, 100)
  ))

  expect_s3_class(host, "vi_host_genome")
  expect_named(host, c("sequences", "genes", "sources", "metadata", "cache", "colors"))
  expect_equal(host$sequences$chr, c("2", "1"))
  expect_null(host$genes)
  expect_identical(names(host$colors), c("2", "1"))
  expect_equal(length(host$colors), 2)
})

test_that("host_genome reads a local FASTA index", {
  fai <- tempfile(fileext = ".fai")
  writeLines(c(paste("scaffold_1", 1200, 0, 80, 81, sep = "\t"), paste("scaffold_2", 900, 0, 80, 81, sep = "\t")), fai)

  host <- host_genome(fai)

  expect_equal(host$sequences$chr, c("scaffold_1", "scaffold_2"))
  expect_equal(host$sequences$end, c(1200, 900))
})




test_that("virus_genome keeps an accession separate from the plotting sequence name", {
  virus <- virus_genome(
    name = "HPV16",
    length = 7904,
    accession = "NC_003977.2",
    features = virus_features(feature = "E6", start = 104, end = 559, type = "gene")
  )

  expect_identical(virus$name, "HPV16")
  expect_identical(virus$sequences$chr, "HPV16")
  expect_identical(virus$metadata$accession, "NC_003977.2")
  expect_identical(virus$metadata$source, "NCBI")
  expect_identical(virus$features$chr, "HPV16")
})


test_that("UCSC human primary chromosomes use natural genomic order", {
  input <- data.frame(
    chr = c("1", "10", "2", "X", "Y", "22", "3"),
    start = 0,
    end = 100,
    stringsAsFactors = FALSE
  )

  ordered <- filter_primary_chromosomes(input, host = "hg38", assembly = "hg38")

  expect_identical(ordered$chr, c("1", "2", "3", "10", "22", "X", "Y"))
})

test_that("ideogram defaults to host and virus genomic axes", {
  expect_identical(track_ideogram()$axis, "all")
})
