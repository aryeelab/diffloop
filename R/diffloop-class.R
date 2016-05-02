#' @include diffloop.R
NULL

#' A class to represent ChIA-PET interaction data and annotations
#'
#' @slot anchors A GRanges object describing loop anchor locations
#' @slot interactions A matrix. Each row is an interaction between two anchors
#' @slot counts A matrix with the number paired-end reads per loop per sample
#' @slot colData A data.frame with features (colums) for each sample (rows)
#' @slot rowData A data.frame with features (columns) for each loop (rows)
#' 
#' @import GenomicRanges
#' @export
loops <- setClass("loops", slots = c(anchors = "GRanges", 
    interactions = "matrix", counts = "matrix", colData = "data.frame", rowData = "data.frame"))

#' A class to represent ChIA-PET interaction data and an edgeR fit. 
#'
#' @slot loops A loops object with anchors, interactions, counts, colData,and rowData
#' @slot fit An edgeR fit from running the loopFit function
#' @import edgeR
#' @export
loopfit <- setClass("loopfit", slots = c(loops = "loops", 
    fit = "DGEGLM"))
