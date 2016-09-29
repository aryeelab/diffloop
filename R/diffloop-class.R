#' @include diffloop.R
NULL

setClassUnion("nim", c("numeric", "integer", "matrix"))

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
    interactions = "nim", counts = "nim",
    colData = "data.frame", rowData = "data.frame"))

