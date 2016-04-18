#' @include diffloop-class.R
NULL

#' chr1:36000000-36300000 loops
#'
#' A loopdata object containing unique 108 loops with 27 anchors
#' for 6 samples and corresponding colData
#'
#' @format A small loopdata object 
#' \describe{
#'   \item{anchors}{GRanges object of anchor locations}
#'   \item{loops}{indexes of interactions }
#'   \item{samples}{Two replicates each of jurkat, naive, and primed cells}
#'   \item{colData}{Groups identifying cell type and unnormalized sizeFactors}
#'   ...
#' }
#' @return A loopdata object
#' @source subsetRegion(full,GRanges(c('1'),IRanges(c(36000000),c(36300000))))
"jpn_chr1reg"