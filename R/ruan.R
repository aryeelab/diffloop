#' @include sugar.R
NULL

#' Compute Chromatin Contact Domains (CCDs)
#'
#' \code{callCCDs} determines regions 
#'
#' This funciton returns a GRanges object of regions determined to be Chromatin Contact
#' Domains as defined in the Tang et al. 2015 paper from the Ruan group. Users
#' can choose to weight the loops by the total number PETs (across all samples)
#' or not and what percent. For details of this method, see page 12 of the supplement
#' of PMID:26686651. Make sure there are only loops within a chromosome before calling this.
#'
#' @param lo A loops object
#' @param petWeights Boolean to weight loop coverage by number of PETs. Default = FALSE
#' @param lowCoveragePercentile Percentile of low coverage to be dropped. Default = 0.05
#'
#' @return A GRanges object of called Chromatin Contact domains
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' lo <- subsetLoops(loops.small, c(1,2,5,6,7,8,9,27,69))
#' ccd <- callCCDs(lo, petWeights = TRUE, lowCoveragePercentile = 0.5)
#' 
#' @importFrom stats quantile
#' 
#' @export
setGeneric(name = "callCCDs", def = function(lo, petWeights = FALSE,lowCoveragePercentile = 0.05) standardGeneric("callCCDs"))

#' @rdname callCCDs
setMethod(f = "callCCDs", def = function(lo, petWeights = FALSE, lowCoveragePercentile = 0.05) {
    sdf <- summary(lo)
    loopSpan <- GRanges(data.frame(chr = sdf$chr_1, start = sdf$start_1, end = sdf$end_2))
    if(petWeights) {
        mcols(loopSpan)$tpc <- rowSums(lo@counts)
        ans <- GRanges(GenomicRanges::coverage(loopSpan, weight = "tpc"))
    } else {
        ans <- GRanges(GenomicRanges::coverage(loopSpan))
    }
    low <- quantile(mcols(ans)$score, probs = c(lowCoveragePercentile))
    return(reduce(ans[mcols(ans)$score >= low]))
})