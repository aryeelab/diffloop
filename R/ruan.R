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


#' Retain loops that have anchors in two specified regions
#'
#' \code{subsetRegionAB} returns a loops object where one anchor 
#' maps to regionA and the other maps to region B
#'
#'
#' @param dlo A loops object
#' @param regionA A GRanges object
#' @param regionB A GRanges object
#'
#' @return A loops object
#'
#' @examples
#' # Return the width for loops 
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' regA <- GRanges(c('1'),IRanges(c(36000000),c(36100000)))
#' regB <- GRanges(c('1'),IRanges(c(36200000),c(36300000)))
#' splits <- subsetRegionAB(dlo, regA, regB)
#'
#' @export
setGeneric(name = "subsetRegionAB", def = function(dlo, regionA, regionB) standardGeneric("subsetRegionAB"))

#' @rdname subsetRegionAB
setMethod(f = "subsetRegionAB", signature = c("loops", "GRanges", "GRanges"),
          definition = function(dlo, regionA, regionB) {
    sdf <- summary(dlo)
    lAnchor <- makeGRangesFromDataFrame(sdf, seqnames.field = "chr_1", start.field = "start_1", end.field = "end_1")
    rAnchor <- makeGRangesFromDataFrame(sdf, seqnames.field = "chr_2", start.field = "start_2", end.field = "end_2")
    numKeep <- union(intersect(queryHits(findOverlaps(lAnchor, regionA)), queryHits(findOverlaps(rAnchor, regionB))),
                     intersect(queryHits(findOverlaps(rAnchor, regionA)), queryHits(findOverlaps(lAnchor, regionB))))
    
    return(subsetLoops(dlo, numKeep))
})

#' Compute boundary scores for genomic loci in between anchors
#'
#' \code{computeBoundaryScores} determines the
#'  boundary scores corresponding to the Genomic region either
#' in bewteen pairs of anchors. To achieve this, the number of PETs
#' for a set of samples (default is all = 0) is summed over a window
#' (default 1MB) on the left (A) and the right (B) of gap. Thus sum 
#' of the number of PETs in these windows is devided by the number 
#' of PETs that span the two loci, plus 1 (C). A larger value
#' corresponds to a stronger boundary. 
#' 
#' Warning: this function is slow; there is a progress bar outputted to give 
#' an anticipated runntime. 
#'
#' @param dlo A loops object
#' @param samples = 0 Vector indexing which samples should be used. 0 is all 
#' @param windowSize = 500000 BP length on left and right of putative boundary to define A/B
#'
#' @return A GRanges object with genomic loci and boundary scores in the mcols
#'
#' @examples
#' # Return the width for loops 
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' BS <- computeBoundaryScores(loops.small, samples = 0, windowSize = 500000)
#'
#' @importFrom pbapply pbsapply
#' @export
setGeneric(name = "computeBoundaryScores", def = function(dlo, samples = 0, windowSize = 500000) standardGeneric("computeBoundaryScores"))

#' @rdname computeBoundaryScores
setMethod(f = "computeBoundaryScores", signature = c("loops", "ANY", "ANY"),
          definition = function(dlo, samples = 0, windowSize = 500000) {
    stopifnot(min(samples) >= 0, max(samples) <= dim(dlo)[3]) # Not a valid entry for samples
    stopifnot(windowSize > 0)
    if(samples[1] == 0) samples <- c(1:as.numeric(dim(dlo)[3]))
    gaps <- data.frame(GenomicRanges::gaps(dlo@anchors))
    ns <- c("seqnames", "start", "end")
    aRegs <- GRanges(setNames(data.frame(gaps[,1], gaps$start - windowSize, gaps$start), ns))
    bRegs <- GRanges(setNames(data.frame(gaps[,1], gaps$end , gaps$end + windowSize), ns))

    ovA <- findOverlaps(aRegs, dlo@anchors)
    ovB <- findOverlaps(bRegs, dlo@anchors)
    qA <- queryHits(ovA)
    qB <- queryHits(ovB)
    
    tru <- intersect(qA, qB) # requires anchors in each to remove ends of chromosomes
    scores <- pbsapply(tru, function(idx){
        Aanc <- subjectHits(ovA[qA == idx,])
        Banc <- subjectHits(ovB[qB == idx,])
        A <- sum(dlo@counts[dlo@interactions[,1] %in% Aanc & dlo@interactions[,2] %in% Aanc,samples])
        B <- sum(dlo@counts[dlo@interactions[,1] %in% Banc & dlo@interactions[,2] %in% Banc,samples])
        C <- sum(sum(dlo@counts[dlo@interactions[,1] %in% Aanc & dlo@interactions[,2] %in% Banc,samples]),
                 sum(dlo@counts[dlo@interactions[,1] %in% Banc & dlo@interactions[,2] %in% Aanc,samples]))
        (A + B)/(C + 1)
    })
    regs <- GenomicRanges::gaps(dlo@anchors)[tru]
    mcols(regs)$BoundaryScores <- scores
    return(regs)
})
