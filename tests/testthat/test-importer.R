caller_fixture <- function(file) {
  candidates <- c(
    file.path("inst", "extdata", file),
    file.path("..", "..", "inst", "extdata", file),
    system.file("extdata", file, package = "virolink")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("Could not find caller fixture: ", file, call. = FALSE)
  }
  path
}

test_that("CTAT-VIF output imports to vi_integrations", {
  path <- caller_fixture("CTAT-VIF_test.GRCH38.insertion_candidates.tsv")

  integrations <- read_ctat_vif_integrations(path)

  expect_s3_class(integrations, "vi_integrations")
  expect_true(all(c(
    "host_chr", "host_pos", "host_strand", "virus_chr", "virus_pos",
    "virus_strand", "support_reads", "raw_record_id", "breakpoint_type",
    "prelim_support_reads", "split_reads", "span_reads", "total_reads"
  ) %in% colnames(integrations)))
  expect_false("sample" %in% colnames(integrations))
  expect_false("method" %in% colnames(integrations))
  expect_false("caller" %in% colnames(integrations))
  expect_equal(integrations$host_chr[1], "11")
  expect_equal(integrations$host_pos[1], 117344940)
  expect_equal(integrations$host_strand[1], "+")
  expect_equal(integrations$virus_chr[1], "HPV16")
  expect_equal(integrations$virus_pos[1], 1215)
  expect_equal(integrations$virus_strand[1], "+")
  expect_equal(integrations$support_reads[1], 112)
  expect_equal(integrations$raw_record_id[1], "chr11~117344940~+~HPV16~1215~+")
  expect_equal(integrations$breakpoint_type[1], "Split")
  expect_equal(integrations$prelim_support_reads[1], 117)
  expect_equal(integrations$split_reads[1], 33)
  expect_equal(integrations$span_reads[1], 79)
  expect_equal(integrations$total_reads[1], 112)
})

test_that("CTAT-VIF importer only adds sample when explicitly provided", {
  path <- caller_fixture("CTAT-VIF_test.GRCH38.insertion_candidates.tsv")

  integrations <- read_ctat_vif_integrations(path, sample = "tumor")

  expect_s3_class(integrations, "vi_integrations")
  expect_equal(unique(integrations$sample), "tumor")
})

test_that("VirusFinder2 output imports to vi_integrations", {
  skip_if_not_installed("readxl")
  path <- caller_fixture("VirusFinder2_integration-sites.tsv")

  integrations <- read_virusfinder2_integrations(path)

  expect_s3_class(integrations, "vi_integrations")
  expect_true(all(c(
    "host_chr", "host_pos", "host_strand", "virus_chr", "virus_pos",
    "virus_strand", "support_reads"
  ) %in% colnames(integrations)))
  expect_false("sample" %in% colnames(integrations))
  expect_false("method" %in% colnames(integrations))
  expect_false("caller" %in% colnames(integrations))
  expect_equal(integrations$host_chr[1], "1")
  expect_equal(integrations$host_pos[1], 24020709)
  expect_equal(integrations$host_strand[1], "+")
  expect_equal(integrations$virus_chr[1], "HBV")
  expect_equal(integrations$virus_pos[1], 1)
  expect_equal(integrations$virus_strand[1], "+")
  expect_equal(integrations$support_reads[1], 25)
})

test_that("Virus-Clip output imports to vi_integrations and keeps confidence", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    paste(
      c(
        "Chromosome 1", "Position 1", "Strand 1",
        "Chromosome 2", "Position 2", "Strand 2",
        "#Supporting reads (pair+softclipping)", "confidence"
      ),
      collapse = "\t"
    ),
    paste(
      c("chr1", "24020709", "+", "HBV", "1", "+", "25", "high"),
      collapse = "\t"
    )
  ), path)

  integrations <- read_virus_clip_integrations(path)

  expect_s3_class(integrations, "vi_integrations")
  expect_true(all(c(
    "host_chr", "host_pos", "host_strand", "virus_chr", "virus_pos",
    "virus_strand", "support_reads", "confidence"
  ) %in% colnames(integrations)))
  expect_equal(integrations$confidence[1], "high")
  expect_equal(integrations$support_reads[1], 25)
})

test_that("caller importer can auto-detect bundled caller outputs", {
  skip_if_not_installed("readxl")
  ctat_path <- caller_fixture("CTAT-VIF_test.GRCH38.insertion_candidates.tsv")
  vf_path <- caller_fixture("VirusFinder2_integration-sites.tsv")

  expect_s3_class(read_caller_integrations(ctat_path), "vi_integrations")
  expect_s3_class(read_caller_integrations(vf_path), "vi_integrations")
})

test_that("caller importer auto-detects Virus-Clip output", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    paste(
      c(
        "Chromosome 1", "Position 1", "Strand 1",
        "Chromosome 2", "Position 2", "Strand 2",
        "#Supporting reads (pair+softclipping)", "confidence"
      ),
      collapse = "\t"
    ),
    paste(
      c("chr1", "24020709", "+", "HBV", "1", "+", "25", "high"),
      collapse = "\t"
    )
  ), path)

  integrations <- read_caller_integrations(path, format = "auto")

  expect_s3_class(integrations, "vi_integrations")
  expect_equal(integrations$confidence[1], "high")
})

test_that("BSVF output imports to vi_integrations", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    paste(
      c(
        "Chr", "breakpoint", "virus-start", "virus-end",
        "virusstrand", "how-many-reads-support", "cluster-name"
      ),
      collapse = "\t"
    ),
    paste(
      c("chr1", "24020709", "200", "300", "+", "18", "cluster1"),
      collapse = "\t"
    )
  ), path)

  integrations <- read_bsvf_integrations(path)

  expect_s3_class(integrations, "vi_integrations")
  expect_true(all(c(
    "host_chr", "host_pos", "host_strand", "virus_chr", "virus_pos",
    "virus_strand", "support_reads", "virus_start", "virus_end",
    "raw_record_id", "cluster_name"
  ) %in% colnames(integrations)))
  expect_equal(integrations$host_chr[1], "1")
  expect_equal(integrations$host_pos[1], 24020709)
  expect_equal(integrations$host_strand[1], "*")
  expect_equal(integrations$virus_chr[1], "virus")
  expect_equal(integrations$virus_pos[1], 250)
  expect_equal(integrations$virus_strand[1], "+")
  expect_equal(integrations$support_reads[1], 18)
  expect_equal(integrations$cluster_name[1], "cluster1")
})

test_that("caller importer auto-detects BSVF output", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    paste(
      c(
        "Chr", "breakpoint", "virus-start", "virus-end",
        "virusstrand", "how-many-reads-support", "cluster-name"
      ),
      collapse = "\t"
    ),
    paste(
      c("chr1", "24020709", "200", "300", "+", "18", "cluster1"),
      collapse = "\t"
    )
  ), path)

  integrations <- read_caller_integrations(path, format = "auto")

  expect_s3_class(integrations, "vi_integrations")
  expect_equal(integrations$virus_pos[1], 250)
})
