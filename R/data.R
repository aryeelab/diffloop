#' @include diffloop-class.R
NULL

#' chr1:36000000-36300000 loops
#'
#' A loops object containing unique 108 loops with 27 anchors
#' for 6 samples and corresponding colData/rowData
#'
#' @format A small loops object 
#' \describe{
#'   \item{anchors}{GRanges object of anchor locations}
#'   \item{loops}{indexes of interactions }
#'   \item{samples}{Two replicates each of jurkat, naive, and primed cells}
#'   \item{colData}{Groups identifying cell type and unnormalized sizeFactors}
#'   \item{rowData}{Base initialization with only loopWidth values}
#'   ...
#' }
#' @return A loops object
#' @source subsetRegion(loops,GRanges(c('1'),IRanges(c(36000000),c(36300000))))
"loops.small"

#' Human protein coding genes
#'
#' A GRanges object with the human protein-coding genes
#'
#' @format A GRanges object 
#' \describe{
#'   \item{seqnames}{Chromosomes without "chr"}
#'   \item{ranges}{start/end loci}
#'   \item{strand}{not specified ('*' everywhere)}
#'   \item{id}{Gene Name}
#'   ...
#' }
#' @return A GRanges object
#' @source biomaRt July 2015 stable build
"human.genes"

#' Human/mouse exon locations
#'
#' A dataframe used for plotting annotation for human
#' and mouse. Each loaded .rda has the same variable
#' called "geneinfo" (so don't co-load these), but 
#' the files differ by an m orh
#'
#' @format A GRanges object 
#' \describe{
#'   \item{chrom}{Chromosomes without "chr"}
#'   \item{start}{exon start location}
#'   \item{stop}{exon end location}
#'   \item{gene}{Gene Name}
#'   \item{score}{dummy column there for sushi}
#'   \item{strand}{+1 or -1 to indicate side of DNA}
#'   ...
#' }
#' @return A data.frame
#' @source biomaRt July 2015 stable build
"geneinfo"
