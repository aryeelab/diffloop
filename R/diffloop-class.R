#' @include diffloop.R
NULL

#' A class to represent ChIA-PET interaction data.
#'
#' @slot anchors A GRanges object describing loop anchor locations
#' @slot loops A matrix. Each row is an interaction between two anchors
#' @slot counts A matrix with the number paired-end reads per loop per sample
#' @slot colData A data.frame with features (colums) for each sample (rows)
#' @import GenomicRanges
#' @export
loopdata <- setClass("loopdata", slots = c(anchors = "GRanges", 
    loops = "matrix", counts = "matrix", colData = "data.frame"))

#' A class to represent ChIA-PET interaction data and an edgeR fit. 
#'
#' @slot loopdata A loopdata with anchors, loops, counts, and colData
#' @slot fit An edgeR fit from running the loopFit function
#' @import edgeR
#' @export
loopfit <- setClass("loopfit", slots = c(loopdata = "loopdata", 
    fit = "DGEGLM"))

#' A class to represent ChIA-PET interaction data and an association
#' test for the loops.
#'
#' @slot loopdata A loopdata object with anchors, loops, counts, and colData
#' @slot results A data.frame of association values per loop in the loopdata
#' @export
looptest <- setClass("looptest", slots = c(loopdata = "loopdata", 
    results = "data.frame"))


