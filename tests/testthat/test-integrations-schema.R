test_that("vi_integrations core schema only requires integration-event fields", {
  integrations <- as_integrations(data.frame(
    host_chr = "chr1",
    host_pos = 100,
    host_strand = "*",
    virus_chr = "HBV",
    virus_pos = 10,
    virus_strand = ".",
    support_reads = 5
  ))

  expect_s3_class(integrations, "vi_integrations")
  expect_equal(integrations$host_chr, "1")
  expect_equal(integrations$host_strand, "*")
  expect_equal(integrations$virus_chr, "HBV")
  expect_equal(integrations$virus_strand, "*")
  expect_equal(integrations$support_reads, 5)
  expect_false("sample" %in% colnames(integrations))
  expect_false("method" %in% colnames(integrations))
})

test_that("vi_integrations uses star for unknown strands", {
  integrations <- as_integrations(data.frame(
    host_chr = c("1", "1"),
    host_pos = c(100, 200),
    host_strand = c(".", "unknown"),
    virus_chr = "HBV",
    virus_pos = c(10, 20),
    virus_strand = c(".", "?"),
    support_reads = c(5, 6)
  ))

  expect_equal(integrations$host_strand, c("*", "*"))
  expect_equal(integrations$virus_strand, c("*", "*"))
})

test_that("vi_integrations requires support_reads", {
  expect_error(
    as_integrations(data.frame(
      host_chr = "1",
      host_pos = 100,
      host_strand = "+",
      virus_chr = "HBV",
      virus_pos = 10,
      virus_strand = "-"
    )),
    "Integration data must contain columns"
  )
})

test_that("vi_integrations does not accept legacy alias-only columns", {
  expect_error(
    as_integrations(data.frame(
      chr = "1",
      host_loc = 100,
      viral_chr = "HBV",
      viral_loc = 10,
      viral_strand = "-",
      reads = 5
    )),
    "Integration data must contain columns"
  )
})
