library(GenomicRanges)

context("Data import functions")

test_that("BED file is imported to GRanges correctly", {
    bed_file <- system.file("extdata", "Jurkat_CTCF_chr1.narrowPeak", package="diffloop")
    gr <- bedToGRanges(bed_file)
    expect_equal(length(gr), 5509)
    expect_equal(as.character(seqnames(gr)[1]), "chr1")
    expect_equal(start(gr)[5509], 249240118)
    expect_equal(gr$id[100], "Jurkat_CTCF_peak_100")
})

test_that("RDA looks good and CCDs work as expected", {
    rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
    load(rda)
    lo <- subsetLoops(loops.small, c(1,2,5,6,7,8,9,27,69))
    ccd <- callCCDs(lo, petWeights = TRUE, lowCoveragePercentile = 0.5)
    expect_equal(length(ccd), 3)
    ccd <- callCCDs(lo, petWeights = FALSE, lowCoveragePercentile = 0.5)
    expect_equal(length(ccd), 2)
})