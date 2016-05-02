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
#' @param dlo A loops object
#' @param nloops A numeric vector of counts to be considered
#'
#' @return A data.frame
#'
#' @examples
#' # Determine what samples have loops with 1-20 counts
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' nLoops <- numLoops(loops.small, 1:20)
#' 
#' # Determine what samples loops with 1-10 counts by default
#' nLoops <- numLoops(loops.small)
#' @export
setGeneric(name = "numLoops", def = function(dlo, nloops = 1:10) standardGeneric("numLoops"))
.numLoops <- function(dlo, nloops) {
    nres <- sapply(nloops, function(n) colSums(dlo@counts > n))
    return(as.data.frame(cbind(nloops, t(nres))))
}

#' @rdname numLoops
setMethod(f = "numLoops", signature = c("loops", "numeric"), 
    definition = function(dlo, nloops) {
        .numLoops(dlo, nloops)
    })

#' @rdname numLoops
setMethod(f = "numLoops", signature = c("loops", "missing"), 
    definition = function(dlo, nloops) {
        .numLoops(dlo, 1:10)
    })


#' Remove self loops
#'
#' \code{removeSelfLoops} removes instances where a loop is observed
#' between the same anchor
#'
#' This function removes loops from the \code{interactions} slot that reference
#' the same index of the \code{anchors} slot.
#'
#' @param dlo A loops object
#'
#' @return A loops object
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' jpn_unique <- removeSelfLoops(loops.small)
#' 
#' @export
setGeneric(name = "removeSelfLoops", def = function(dlo) standardGeneric("removeSelfLoops"))

#' @rdname removeSelfLoops
setMethod(f = "removeSelfLoops", signature = c("loops"), definition = function(dlo) {
    return(subsetLoops(dlo, dlo@interactions[, 1] != dlo@interactions[, 
        2]))
})

# Return Boolean Vector for loops in loops if both anchors
# are unique... Internal method
setGeneric(name = "uniqueLoops", def = function(dlo) standardGeneric("uniqueLoops"))
setMethod(f = "uniqueLoops", signature = c("loops"), definition = function(dlo) {
    return(((!is.na(dlo@interactions[, 1]) & !is.na(dlo@interactions[, 
        2])) & (dlo@interactions[, 1] != dlo@interactions[, 2])))
})

# Determine loop type as either unique (with called anchors),
# self, one unique anchor = single, or no unique anchors
# Internal method
setGeneric(name = "classifyLoops", def = function(dlo) standardGeneric("classifyLoops"))
setMethod(f = "classifyLoops", signature = c("loops"), definition = function(dlo) {
    x <- dlo@interactions
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
#' none has no unique anchors. If using the \code{loopsMake} pipeline, only
#' self and unique loops will be observed when running this function
#'
#' @param dlo A loops object
#'
#' @return A data.frame
#'
#' @examples
#' # Return loop metrics for number of each type for each sample 
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' loopMetrics(loops.small)
#' 
#' @export
setGeneric(name = "loopMetrics", def = function(dlo) standardGeneric("loopMetrics"))

#' @rdname loopMetrics
setMethod(f = "loopMetrics", signature = c("loops"), definition = function(dlo) {
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
#' only intrachromosomal loops. Also, self-loops will return a postive number
#' that is the inter-anchor width. These loops should be handled using the 
#' removeSelfLoops() function.
#'
#' @param dlo A loops object
#'
#' @return An integer vector
#'
#' @examples
#' # Return the width for loops 
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' w <- loopWidth(loops.small)
#'
#' @export
setGeneric(name = "loopWidth", def = function(dlo) standardGeneric("loopWidth"))

#' @rdname loopWidth
setMethod(f = "loopWidth", signature = c("loops"), definition = function(dlo) {
    w <- start(dlo@anchors[dlo@interactions[, 2]]) - end(dlo@anchors[dlo@interactions[, 
        1]]) + 1
    return(w)
})

#' Subset loops
#'
#' \code{subsetLoops} restricts the loops and counts matrix to only those 
#' specified by \code{idxa}, either numerically or logically
#'
#' This function returns a loops object where the loops are retained only
#' if they meet a logical criteria or are included in the numeric vector of
#' \code{idxa}. Only the anchors that reference a loop in the subsetted 
#' loops object are retained. 
#'
#' @param dlo A loops object
#' @param idxa A numeric vector or logical vector
#'
#' @return A loops object
#'
#' @examples
#' # Return the first 10 loops
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' #' ten <- subsetLoops(loops.small, 1:10)
#' 
#' # Subset loops with widths greater than 10000
#' big <- subsetLoops(loops.small, loopWidth(loops.small) >= 10000)
#'
#' @export
setGeneric(name = "subsetLoops", def = function(dlo, idxa) standardGeneric("subsetLoops"))

.subsetLoops <- function(dlo, idxa) {
    # Keep Unique Interactions and their counts
    slot(dlo, "interactions", check = TRUE) <- as.matrix(dlo@interactions[idxa, , drop = FALSE])
    slot(dlo, "counts", check = TRUE) <- as.matrix(dlo@counts[idxa, , drop = FALSE])
    nRowData <- dlo@rowData[idxa, , drop = FALSE]
    row.names(nRowData) <- NULL
    slot(dlo, "rowData", check = TRUE) <- nRowData
    return(cleanup(dlo))
}

#' @rdname subsetLoops
setMethod(f = "subsetLoops", signature = c("loops", "logical"), 
    definition = function(dlo, idxa) {
        if (all(idxa)) 
            return(dlo) else return(.subsetLoops(dlo, idxa))
    })

#' @rdname subsetLoops
setMethod(f = "subsetLoops", signature = c("loops", "numeric"), 
    definition = function(dlo, idxa) {
        .subsetLoops(dlo, idxa)
    })

#' Filter loops
#'
#' \code{filterLoops} filters out loops that aren't wide, aren't prevalent
#' within samples or prevalent between samples
#'
#' Function that restricts loops in a loops object. \code{width} specifies 
#' the minimum width between anchors. Default is zero. \code{nreplicates}
#' restricts loops to at least this specified amount of counts is present
#' in at least one sample. Instead of \code{nreplicates} being present in only
#' one sample, \code{nsamples} specifies how many individual samples that a
#' loop must have \code{nreplicates} in to be included after filtering. 
#'
#' @param dlo A loops object
#' @param width Minimum loop width 
#' @param nreplicates Minimum number of counts per loop
#' @param nsamples Minimum number of samples per loop per counts
#'
#' @return A loops object
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' # Restrict loops to > 5kb width
#' filtered.jpn1 <- filterLoops(loops.small, 5000, 0, 0)
#' # Restrict loops to > 5kb width and have >= 3 replicates in >= 1 sample
#' filtered.jpn2 <- filterLoops(loops.small, 5000, 3, 1)
#' # Restrict loops to > 10kb width and have >= 3 replicates in >= 2 samples
#' filtered.jpn3 <- filterLoops(loops.small, 10000, 3, 2)

#' @export
setGeneric(name = "filterLoops", def = function(dlo, width = 0, 
    nreplicates = 0, nsamples = 1) standardGeneric("filterLoops"))

#' @rdname filterLoops
setMethod(f = "filterLoops", definition = function(dlo, width = 0, 
    nreplicates = 0, nsamples = 1) {
    n <- dim(dlo@interactions)[1]
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
#' @param dlo A loops object
#' @param genesGR A GRanges object of genes with mcol 'id'
#'
#' @return A matrix of comma separated gene names
#'
#' @examples
#' # Determine the genes housed in the loops from our example
#' genes <- getHumanGenes()
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' loops.small <- loopGenes(loops.small,genes)
#' 
#' @import GenomicRanges
#' @importFrom IRanges IRanges 
#' @export
setGeneric(name = "loopGenes", def = function(dlo, genesGR) standardGeneric("loopGenes"))

#' @rdname loopGenes
setMethod(f = "loopGenes", signature = c("loops", "GRanges"), 
    definition = function(dlo, genesGR) {
        values <- apply(dlo@interactions, 1, function(interaction) {
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
#' \code{intrachromosomal} restricts interactions to those where anchors are observed
#' on the same chromosomes
#'
#' This function subsets the \code{loops} object into only those interactions that
#' have both anchors on the same chromosome
#'
#' @param dlo A loops object
#'
#' @return A loops object where all loops are on the same chromosome.
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' 
#' # Compute number of interactions on same chromosome
#' dim(intrachromosomal(loops.small))
#' samechromo <- intrachromosomal(loops.small)
#' 
#' @export
setGeneric(name = "intrachromosomal", def = function(dlo) standardGeneric("intrachromosomal"))

#' @rdname intrachromosomal
setMethod(f = "intrachromosomal", signature = c("loops"), definition = function(dlo) {
    idx <- as.character(dlo@anchors[dlo@interactions[, 1]]@seqnames) == 
        as.character(dlo@anchors[dlo@interactions[, 2]]@seqnames)
    return(subsetLoops(dlo, idx))
})


#' Loops between chromosomes
#'
#' \code{interchromosomal} restricts loops to those where anchors are observed
#' on different chromosomes 
#'
#' This function subsets the \code{loops} object into only those loops that
#' have anchors on different chromosomes
#'
#' @param dlo A loops object
#' @return A loops object with all loops on different chromosomes
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' 
#' # Compute number of interactions on same chromosome
#' dim(intrachromosomal(loops.small))
#' samechromo <- intrachromosomal(loops.small)
#' 
#' # Compute number of interactions on same chromosome
#' # dim(interchromosomal(loops.small))
#' # This will throw and error since the toy only has intrachromosomal loops
#' 
#' @export
setGeneric(name = "interchromosomal", def = function(dlo) standardGeneric("interchromosomal"))

#' @rdname interchromosomal
setMethod(f = "interchromosomal", signature = c("loops"), definition = function(dlo) {
    idx <- as.character(dlo@anchors[dlo@interactions[, 1]]@seqnames) != 
        as.character(dlo@anchors[dlo@interactions[, 2]]@seqnames)
    return(subsetLoops(dlo, idx))
    
})
