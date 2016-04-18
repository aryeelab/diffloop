#' diffloop: A package for differential DNA loop calling from ChIA-PET data
#'
#' The diffloop package contains a suite of tools and S4 data objects to
#' efficiently facilitate the analysis of ChIA-PET datasets. Key features
#' include differential loop calling, visualization of looping in regions,
#' quality-control metrics, and principal component analysis across 
#' experiments. 
#' 
#' @section diffloop classes:
#' Three classes mostly comprise the methodology in \code{diffloop}. First,
#' \code{loopdata} is a basic structure that contains one or more ChIA-PET
#' experiments, \code{loopfit} links an \code{edgeR} fit to a \code{loopdata}
#' and currently has little functionality except for generating a 
#' \code{looptest} object, which links summary statistics per loop with a 
#' \code{loopdata} object.
#'   
#' @docType package
#' @name diffloop
#' @import methods
NULL
