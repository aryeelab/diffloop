#' @include pipeline_io.R
NULL

#' Per-sample loop quantities
#'
#' \code{numLoops} counts number of loops for each sample based on the
#' index of nloops and returns a data.frame
#'
#' This function shows the number of unique loops with at least
#' \code{nloops} in counts. Can be used to quickly visualize relative
#' sequencing depth between samples
#'
#' @param dlo A loopdata object
#' @param nloops A numeric vector of loops to be considered
#'
#' @return A data.frame
#'
#' @examples
#' # Determine what samples have loops with 1-20 counts
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' nLoops <- numLoops(jpn_chr1reg, 1:20)
#' 
#' # Determine what samples loops with 1-10 counts by default
#' nLoops <- numLoops(jpn_chr1reg)
#' @export
setGeneric(name = "numLoops", def = function(dlo, nloops = 1:10) standardGeneric("numLoops"))
.numLoops <- function(dlo, nloops) {
    nres <- sapply(nloops, function(n) colSums(dlo@counts > n))
    return(as.data.frame(cbind(nloops, t(nres))))
}

#' @rdname numLoops
setMethod(f = "numLoops", signature = c("loopdata", "numeric"), 
    definition = function(dlo, nloops) {
        .numLoops(dlo, nloops)
    })

#' @rdname numLoops
setMethod(f = "numLoops", signature = c("loopdata", "missing"), 
    definition = function(dlo, nloops) {
        .numLoops(dlo, 1:10)
    })


#' Remove self loops
#'
#' \code{removeSelfLoops} removes instances where a loop is observed
#' between the same anchor
#'
#' This function removes loops from the \code{loops} slot that reference
#' the same index of the \code{anchors} slot.
#'
#' @param dlo A loopdata/looptest object
#'
#' @return A loopdata/looptest object
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_unique <- removeSelfLoops(jpn_chr1reg)
#' 
#' @export
setGeneric(name = "removeSelfLoops", def = function(dlo) standardGeneric("removeSelfLoops"))

#' @rdname removeSelfLoops
setMethod(f = "removeSelfLoops", signature = c("loopdata"), definition = function(dlo) {
    return(subsetLoops(dlo, dlo@loops[, 1] != dlo@loops[, 2]))
})

#' @rdname removeSelfLoops
setMethod(f = "removeSelfLoops", signature = c("looptest"), definition = function(dlo) {
    idx <- dlo@loopdata@loops[, 1] != dlo@loopdata@loops[, 2]
    slot(dlo, "loopdata", check=TRUE) <- subsetLoops(dlo@loopdata, idx)
    slot(dlo, "results", check=TRUE) <- dlo@results[idx,]
    return(dlo)
})

# Return Boolean Vector for loops in loopdata if both anchors
# are unique Internal method
setGeneric(name = "uniqueLoops", def = function(dlo) standardGeneric("uniqueLoops"))
setMethod(f = "uniqueLoops", signature = c("loopdata"), definition = function(dlo) {
    return(((!is.na(dlo@loops[, 1]) & !is.na(dlo@loops[, 2])) & 
        (dlo@loops[, 1] != dlo@loops[, 2])))
})

# Determine loop type as either unique (with called anchors),
# self, one unique anchor = single, or no unique anchors
# Internal method
setGeneric(name = "classifyLoops", def = function(dlo) standardGeneric("classifyLoops"))
setMethod(f = "classifyLoops", signature = c("loopdata"), definition = function(dlo) {
    x <- dlo@loops
    type <- rep(NA, nrow(x))
    unique <- (!is.na(x[, 1]) & !is.na(x[, 2])) & (x[, 1] != 
        x[, 2])
    self <- (!is.na(x[, 1]) & !is.na(x[, 2])) & (x[, 1] == x[, 
        2])
    single <- (!is.na(x[, 1]) & is.na(x[, 2])) | (is.na(x[, 1]) & 
        !is.na(x[, 2]))
    none <- (is.na(x[, 1]) & is.na(x[, 2]))
    type[unique] <- "unique"
    type[self] <- "self"
    type[single] <- "single"
    type[none] <- "none"
    type
})

#' Types of loops
#'
#' \code{loopMetrics} counts number of loops for each sample and returns
#' whether they are single, self, unique, or none
#'
#' This function shows the number of loops for each sample based on four
#' types. Single refers to having only one anchor for a the loop whereas
#' none has no unique anchors. If using the \code{loopdataMake} pipeline, only
#' self and unique loops will be observed when running this function
#'
#' @param dlo A loopdata object
#'
#' @return A data.frame
#'
#' @examples
#' # Return loop metrics for number of each type for each sample 
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' loopMetrics(jpn_chr1reg)
#' 
#' @export
setGeneric(name = "loopMetrics", def = function(dlo) standardGeneric("loopMetrics"))

#' @rdname loopMetrics
setMethod(f = "loopMetrics", signature = c("loopdata"), definition = function(dlo) {
    return(do.call("rbind", by(dlo@counts, classifyLoops(dlo), 
        colSums, na.rm = TRUE)))
})


#' Loop widths
#'
#' \code{loopWidth} returns the width of a loop, which is defined as the
#' distance between the anchors containing a loop
#'
#' This function returns a positive integer value of the number of basepairs
#' that separate two loops. If they are on separate chromosomes, it still
#' returns a value, but it will be non-sensical, so consider subsetting to
#' only intrachromosomal loops.
#'
#' @param dlo A loopdata object
#'
#' @return An integer vector
#'
#' @examples
#' # Return the width for loops 
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' w <- loopWidth(jpn_chr1reg)
#'
#' @export
setGeneric(name = "loopWidth", def = function(dlo) standardGeneric("loopWidth"))

#' @rdname loopWidth
setMethod(f = "loopWidth", signature = c("loopdata"), definition = function(dlo) {
    ret <- rep(NA, nrow(dlo@loops))
    uniqueloops <- dlo@loops
    w <- start(dlo@anchors[uniqueloops[, 2]]) - end(dlo@anchors[uniqueloops[, 
        1]]) + 1
    ret <- abs(w)
    return(ret)
})

#' Subset loops
#'
#' \code{subsetLoops} restricts the loops and counts matrix to only those 
#' specified by \code{idxa}, either numerically or logically
#'
#' This function returns a loopdata object where the loops are retained only
#' if they meet a logical criteria or are included in the numeric vector of
#' \code{idxa}. Only the anchors that reference a loop in the subsetted 
#' loopdata object are retained. 
#'
#' @param dlo A loopdata object
#' @param idxa A numeric vector or logical vector
#'
#' @return A loopdata object
#'
#' @examples
#' # Return the first 10 loops
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' #' ten <- subsetLoops(jpn_chr1reg, 1:10)
#' 
#' # Subset loops with widths greater than 10000
#' big <- subsetLoops(jpn_chr1reg, loopWidth(jpn_chr1reg) >= 10000)
#'
#' @export
setGeneric(name = "subsetLoops", def = function(dlo, idxa) standardGeneric("subsetLoops"))

.subsetLoops <- function(dlo, idxa) {
    # Keep Unique Interactions and their counts
    slot(dlo, "loops", check = TRUE) <- dlo@loops[idxa, , drop = FALSE]
    slot(dlo, "counts", check = TRUE) <- dlo@counts[idxa, , drop = FALSE]
    return(cleanup(dlo))
}

#' @rdname subsetLoops
setMethod(f = "subsetLoops", signature = c("loopdata", "logical"), 
    definition = function(dlo, idxa) {
        if(all(idxa))
            return(dlo)
        else
            return(.subsetLoops(dlo, idxa))
    })

#' @rdname subsetLoops
setMethod(f = "subsetLoops", signature = c("loopdata", "numeric"), 
    definition = function(dlo, idxa) {
        .subsetLoops(dlo, idxa)
    })

#' Filter loops
#'
#' \code{filterLoops} filters out loops that aren't wide, aren't prevalent
#' within samples or prevalent between samples
#'
#' Function that restricts loops in a loopdata object. \code{width} specifies 
#' the minimum width between anchors. Default is zero. \code{nreplicates}
#' restricts loops to at least this specified amount of counts is present
#' in at least one sample. Instead of \code{nreplicates} being present in only
#' one sample, \code{nsamples} specifies how many individual samples that a
#' loop must have \code{nreplicates} in to be included after filtering. 
#'
#' @param dlo A loopdata object
#' @param width Minimum loop width 
#' @param nreplicates Minimum number of counts per loop
#' @param nsamples Minimum number of samples per loop per counts
#'
#' @return A loopdata object
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' # Restrict loops to > 5kb width
#' filtered.jpn1 <- filterLoops(jpn_chr1reg, 5000, 0, 0)
#' # Restrict loops to > 5kb width and have >= 3 replicates in >= 1 sample
#' filtered.jpn2 <- filterLoops(jpn_chr1reg, 5000, 3, 1)
#' # Restrict loops to > 10kb width and have >= 3 replicates in >= 2 samples
#' filtered.jpn3 <- filterLoops(jpn_chr1reg, 10000, 3, 2)

#' @export
setGeneric(name = "filterLoops", def = function(dlo, width = 0, 
    nreplicates = 0, nsamples = 1) standardGeneric("filterLoops"))

#' @rdname filterLoops
setMethod(f = "filterLoops", definition = function(dlo, width = 0, 
    nreplicates = 0, nsamples = 1) {
    n <- dim(dlo@loops)[1]
    idxw <- rep(TRUE, n)
    idxr <- rep(TRUE, n)
    idxn <- rep(TRUE, n)
    idxa <- rep(TRUE, n)
    
    if (width != 0) 
        idxw <- (loopWidth(dlo) >= width)
    if (nreplicates != 0) 
        idxr <- apply(dlo@counts, 1, function(t) any(t >= nreplicates))
    if (nsamples != 1 | nsamples != 0) 
        idxn <- apply(dlo@counts, 1, function(t) sum(as.numeric(t >= 
            nreplicates)) >= nsamples)
    idxa <- (idxw & idxr & idxn)
    ndlo <- subsetLoops(dlo, idxa)
    if (width > 0)
        ndlo <- removeSelfLoops(ndlo)
    return(ndlo)
})

#' Determine genes contained within loops
#'
#' \code{loopGenes} determines all gene bodies partially or fully contained
#' in a loop.
#'
#' Function that annotates all loops. 'NA' if looping between chromosomes.
#' Otherwise, all gene names that are contained within a loop. 'None' if no
#' genes are in the loop body. If there are multiple, the function returns a
#' comma separated list. The length of the object returned by this function
#' should be the same length as the number of rows in the \code{loops} slot.
#'
#' @param dlo A loopdata object
#' @param genesGR A GRanges object of genes with mcol 'id'
#'
#' @return A matrix of comma separated gene names
#'
#' @examples
#' # Determine the genes housed in the loops from our example
#' genes <- getHumanGenes()
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_chr1reg <- loopGenes(jpn_chr1reg,genes)
#' 
#' @import GenomicRanges
#' @importFrom IRanges IRanges 
#' @export
setGeneric(name = "loopGenes", def = function(dlo, genesGR) standardGeneric("loopGenes"))

#' @rdname loopGenes
setMethod(f = "loopGenes", signature = c("loopdata", "GRanges"), 
    definition = function(dlo, genesGR) {
        values <- apply(dlo@loops, 1, function(interaction) {
            chr1 <- as.character(dlo@anchors[interaction[1]]@seqnames)
            chr2 <- as.character(dlo@anchors[interaction[2]]@seqnames)
            if (chr1 != chr2) {
                "NA"
            } else {
                big <- max(dlo@anchors[interaction[2]]@ranges)
                small <- min((dlo@anchors[interaction[1]]@ranges))
                loop <- GRanges(seqnames = c(chr1), ranges = IRanges(start = c(small), 
                  end = c(big)))
                
                # Supposes that the first mcol of the genesGR are the gene
                # names
                geneHits <- mcols(genesGR[findOverlaps(loop, 
                  genesGR)@from])[, 1]
                if (length(geneHits) == 0) {
                  "NONE"
                } else {
                  paste(unique(geneHits), collapse = ", ")
                }
            }
        })
        values <- as.matrix(values)
        colnames(values) <- c("GenesInLoop")
        return(values)
    })

#' Loops within chromosomes
#'
#' \code{intrachromosomal} restricts loops to those where anchors are observed
#' on the same chromosomes
#'
#' This function subsets the \code{loopdata} object into only those loops that
#' have both anchors on the same chromosome
#'
#' @param dlo A loopdata object
#'
#' @return A loopdata object where all loops are on the same chromosome.
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' 
#' # Compute number of loops on same chromosome
#' dim(intrachromosomal(jpn_chr1reg))
#' samechromo <- intrachromosomal(jpn_chr1reg)
#' 
#' @export
setGeneric(name = "intrachromosomal", def = function(dlo) standardGeneric("intrachromosomal"))

#' @rdname intrachromosomal
setMethod(f = "intrachromosomal", signature = c("loopdata"), 
    definition = function(dlo) {
        idx <- as.character(dlo@anchors[dlo@loops[, 1]]@seqnames) == 
            as.character(dlo@anchors[dlo@loops[, 2]]@seqnames)
        return(subsetLoops(dlo, idx))
    })


#' Loops between chromosomes
#'
#' \code{interchromosomal} restricts loops to those where anchors are observed
#' on different chromosomes 
#'
#' This function subsets the \code{loopdata} object into only those loops that
#' have anchors on different chromosomes
#'
#' @param dlo A loopdata object
#' @return A loopdata object with all loops on different chromosomes
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' 
#' # Compute number of loops on same chromosome
#' dim(intrachromosomal(jpn_chr1reg))
#' samechromo <- intrachromosomal(jpn_chr1reg)
#' 
#' # Compute number of loops on same chromosome
#' # dim(interchromosomal(jpn_chr1reg))
#' # This will throw and error since the toy only has intrachromosomal loops
#' 
#' @export
setGeneric(name = "interchromosomal", def = function(dlo) standardGeneric("interchromosomal"))

#' @rdname interchromosomal
setMethod(f = "interchromosomal", signature = c("loopdata"), 
    definition = function(dlo) {
        idx <- as.character(dlo@anchors[dlo@loops[, 1]]@seqnames) != 
            as.character(dlo@anchors[dlo@loops[, 2]]@seqnames)
        return(subsetLoops(dlo, idx))
        
    })
