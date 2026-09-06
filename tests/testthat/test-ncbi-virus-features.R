test_that("GFF3 records are converted to vi_virus_features", {
  gff3 <- paste(
    "##gff-version 3",
    "NC_001526.4\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-E6;gene=E6;product=early protein",
    "NC_001526.4\tRefSeq\tCDS\t120\t800\t.\t+\t0\tID=cds-E6;Parent=gene-E6;gene=E6;product=E6 protein",
    "NC_001526.4\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-E6;Parent=gene-E6",
    sep = "\n"
  )

  result <- parse_ncbi_gff3(gff3, accession = "NC_001526.4")

  expect_s3_class(result, "vi_virus_features")
  expect_equal(names(result), c("feature", "start", "end", "type"))
  expect_equal(result$feature, c("E6", "rna-E6", "E6"))
  expect_equal(result$type, c("gene", "mRNA", "CDS"))
  expect_equal(result$start, c(100, 100, 120))
  expect_equal(result$end, c(900, 900, 800))
})

test_that("NCBI helper constructs a reproducible GFF3 URL", {
  expect_match(
    ncbi_gff3_url("NC_001526.4"),
    "db=nuccore.*id=NC_001526.4.*rettype=gff3.*retmode=text"
  )
})

test_that("virus_info is the public annotation helper", {
  expect_true(is.function(virus_info))
  expect_false("features" %in% names(formals(virus_info)))
})

test_that("virus_features filters downloaded annotation records", {
  gff3 <- paste(
    "##gff-version 3",
    "NC_001526.4\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-E6;gene=E6",
    "NC_001526.4\tRefSeq\tCDS\t120\t800\t.\t+\t0\tID=cds-E6;gene=E6",
    sep = "\n"
  )
  annotation <- list(accession = "NC_001526.4", gff3 = gff3)

  result <- virus_features(annotation, type = "gene")

  expect_s3_class(result, "vi_virus_features")
  expect_equal(result$type, "gene")
  expect_equal(result$feature, "E6")
})

test_that("virus_features expands circular virus features", {
  gff3 <- paste(
    "##gff-version 3",
    "##sequence-region NC_003977.2 1 3215",
    "NC_003977.2\tRefSeq\tgene\t2307\t1623\t.\t+\t.\tID=gene-P;gene=P",
    sep = "\n"
  )
  annotation <- list(
    accession = "NC_003977.2",
    gff3 = gff3,
    metadata = data.frame(seqid = "NC_003977.2", start = 1, end = 3215, stringsAsFactors = FALSE)
  )

  result <- virus_features(annotation, type = "gene")

  expect_s3_class(result, "vi_virus_features")
  expect_equal(nrow(result), 2L)
  expect_equal(result$feature, c("P", "P"))
  expect_equal(result$start, c(1, 2307))
  expect_equal(result$end, c(1623, 3215))
})

test_that("virus_features filters annotation records before standardization", {
  gff3 <- paste(
    "##gff-version 3",
    "##sequence-region NC_003977.2 1 3215",
    "NC_003977.2\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-X;gene=X",
    "NC_003977.2\tRefSeq\tCDS\t2307\t1623\t.\t+\t0\tID=cds-P;gene=P",
    sep = "\n"
  )
  annotation <- list(
    accession = "NC_003977.2",
    gff3 = gff3,
    metadata = data.frame(seqid = "NC_003977.2", start = 1, end = 3215, stringsAsFactors = FALSE)
  )

  expect_warning(
    result <- virus_features(annotation, type = "gene"),
    NA
  )
  expect_equal(result$feature, "X")
  expect_equal(result$type, "gene")
})

test_that("virus_features filters out unrelated invalid records before standardization", {
  gff3 <- paste(
    "##gff-version 3",
    "##sequence-region NC_003977.2 1 3182",
    "NC_003977.2\tRefSeq\tgene\t1376\t1840\t.\t+\t.\tID=gene-X;gene=X",
    "NC_003977.2\tRefSeq\tCDS\t500\t500\t.\t+\t0\tID=cds-invalid;gene=bad",
    sep = "\n"
  )
  annotation <- list(
    accession = "NC_003977.2",
    gff3 = gff3,
    metadata = data.frame(seqid = "NC_003977.2", start = 1, end = 3182, stringsAsFactors = FALSE)
  )

  expect_warning(
    result <- virus_features(annotation, type = "gene"),
    NA
  )
  expect_equal(result$feature, "X")
})

test_that("virus_features validates selected annotation records after filtering", {
  gff3 <- paste(
    "##gff-version 3",
    "##sequence-region NC_003977.2 1 3182",
    "NC_003977.2\tRefSeq\tgene\t1376\t1840\t.\t+\t.\tID=gene-X;gene=X",
    "NC_003977.2\tRefSeq\tCDS\t.\t500\t.\t+\t0\tID=cds-invalid",
    sep = "\n"
  )
  annotation <- list(
    accession = "NC_003977.2",
    gff3 = gff3,
    metadata = data.frame(seqid = "NC_003977.2", start = 1, end = 3182, stringsAsFactors = FALSE)
  )

  expect_warning(
    result <- virus_features(annotation, type = "gene"),
    NA
  )
  expect_equal(result$feature, "X")
})

test_that("virus_features expands circular features that extend past genome length", {
  gff3 <- paste(
    "##gff-version 3",
    "##sequence-region NC_003977.2 1 3182",
    "NC_003977.2\tRefSeq\tgene\t2309\t4807\t.\t+\t.\tID=gene-P;gene=P",
    sep = "\n"
  )
  annotation <- list(
    accession = "NC_003977.2",
    gff3 = gff3,
    metadata = data.frame(seqid = "NC_003977.2", start = 1, end = 3182, stringsAsFactors = FALSE)
  )

  result <- virus_features(annotation, type = "gene")

  expect_equal(nrow(result), 2L)
  expect_equal(result$feature, c("P", "P"))
  expect_equal(result$start, c(1, 2309))
  expect_equal(result$end, c(1625, 3182))
})

test_that("normalize_virus_features expands circular user features with a virus length", {
  features <- data.frame(
    feature = c("P", "S"),
    start = c(2309, 2850),
    end = c(4807, 500),
    type = c("polymerase", "surface"),
    stringsAsFactors = FALSE
  )

  result <- normalize_virus_features(features, virus_length = 3182)

  expect_s3_class(result, "vi_virus_features")
  expect_equal(result$feature, c("S", "P", "P", "S"))
  expect_equal(result$start, c(1, 1, 2309, 2850))
  expect_equal(result$end, c(500, 1625, 3182, 3182))
})
