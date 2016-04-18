#' @include data.R
NULL

#' Link the anchors and loops back together 
#'
#' \code{summarize} takes a \code{loopdata/looptest} object and breaks the 
#' loopdata structure resulting in a \code{data.frame}. 
#'
#' This function returns a \code{data.frame} where the left and right anchors 
#' are visualized together along with the loop width, individual counts, and
#' any anchor meta-data that has been annotated into the anchors GRanges
#' object. When a \code{looptest} object is the input, the additional rows 
#' in the \code{results} slot is added on but the individual counts are 
#' excluded.
#'
#' @param .data A loopdata object to be summarized
#'
#' @return A data.frame
#'
#' @examples
#' # Summarizing the first ten loops in \code{jpn_chr1reg}
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' summarydf <- summarize(jpn_chr1reg[1:10,])
#' # Summarizing the loops and significance results between naive and primed
#' summarylt <- summarize(quickAssoc(jpn_chr1reg[,1:4])[1:10,])
#' @import plyr

#' @export
setMethod(f = "summarize", signature = c("loopdata"), definition = function(.data) {
    dlo <- .data
    # Grab all the left anchors in order of loop occurence
    leftAnchor2 <- as.data.frame(dlo@anchors[dlo@loops[, 1]])
    leftAnchor2 <- subset(leftAnchor2, select = -c(width, strand))
    colnames(leftAnchor2) <- paste(colnames(leftAnchor2), "1", 
        sep = "_")
    colnames(leftAnchor2)[1] <- "chr_1"
    
    # Grab all the right anchors in order of loop occurence
    rightAnchor2 <- as.data.frame(dlo@anchors[dlo@loops[, 2]])
    rightAnchor2 <- subset(rightAnchor2, select = -c(width, strand))
    colnames(rightAnchor2) <- paste(colnames(rightAnchor2), "2", 
        sep = "_")
    colnames(rightAnchor2)[1] <- "chr_2"
    
    # Add the loop features width and counts (per sample)
    loopwidth <- matrix(loopWidth(dlo), ncol = 1)
    colnames(loopwidth) <- "loopwidth"
    cbind(leftAnchor2, rightAnchor2, loopwidth, dlo@counts)
})

#' @rdname summarize-loopdata-method 
setMethod(f = "summarize", signature = c("looptest"), definition = function(.data) {
    dlo <- .data@loopdata
    # Grab all the left anchors in order of loop occurence
    leftAnchor2 <- as.data.frame(dlo@anchors[dlo@loops[, 1]])
    leftAnchor2 <- subset(leftAnchor2, select = -c(width, strand))
    colnames(leftAnchor2) <- paste(colnames(leftAnchor2), "1", 
        sep = "_")
    colnames(leftAnchor2)[1] <- "chr_1"
    
    # Grab all the right anchors in order of loop occurence
    rightAnchor2 <- as.data.frame(dlo@anchors[dlo@loops[, 2]])
    rightAnchor2 <- subset(rightAnchor2, select = -c(width, strand))
    colnames(rightAnchor2) <- paste(colnames(rightAnchor2), "2", 
        sep = "_")
    colnames(rightAnchor2)[1] <- "chr_2"
    
    # Add the loop features width and counts (per sample)
    loopwidth <- matrix(loopWidth(dlo), ncol = 1)
    colnames(loopwidth) <- "loopwidth"
    cbind(leftAnchor2, rightAnchor2, loopwidth, .data@results)
})

# Function that removes all anchors not being referenced in
# loops matrix and updates indices. For internal use only.
setGeneric(name = "cleanup", def = function(dlo) standardGeneric("cleanup"))
setMethod(f = "cleanup", signature = c("loopdata"), definition = function(dlo) {
    if (dim(dlo@counts)[1] == 0) {
        stop("Attempting to subset to empty looptest/loopdata object!")
    }
    # Grab indicies of anchors being referenced in loops
    idf <- data.frame(dlo@loops[, 1], dlo@loops[, 2])
    sdf <- stack(idf)
    udf <- sort(unique(sdf[, "values"]))
    
    # Keep only those anchors that are being used
    newAnchors <- dlo@anchors[udf]
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
    intsdf <- as.data.frame(dlo@loops)
    
    # Update loops indices
    leftmatch <- t(sapply(intsdf$left, function(x) mapping[mapping[, 
        1] == x, ]))
    rightmatch <- t(sapply(intsdf$right, function(x) mapping[mapping[, 
        1] == x, ]))
    
    # Format new indices matrix
    totalupdate <- cbind(unlist(leftmatch[, 2]), unlist(rightmatch[, 
        2]))
    upints <- as.matrix(totalupdate)
    colnames(upints) <- c("left", "right")
    
    # Update values
    slot(dlo, "loops", check = TRUE) <- upints
    slot(dlo, "anchors", check = TRUE) <- newAnchors
    return(dlo)
})

#' Combine nearby anchors into one peak
#'
#' \code{mergeAnchors} combines anchors that are within a user-defined radius
#'
#' This function takes a loopdata object and combines nearby anchors, up to
#' a distance specified by the \code{mergegap}. This likely will cause self
#' loops to form (loop where the left and right anchor are the same), which
#' can either be removed (by default) or retained with \code{selfloops}
#'
#' @param dlo A loopdata object whose anchors will be merged
#' @param mergegap An integer value of the bp between anchors to be merged
#' @param selfloops A logical value to either retain (T) or remove (F) 
#' resulting self-loops after merging anchors
#'
#' @return A loopdata object
#'
#' @examples
#' # Merge anchors within 1kb of each other, keeping self loops 
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' m1kb <- mergeAnchors(jpn_chr1reg, 1000, FALSE)
#'
#' # Merge anchors within 1kb of each other, removing self loops by default
#' m1kb_unique <- mergeAnchors(jpn_chr1reg, 1000)

#' @import GenomicRanges
#' @import reshape2
#' @export
setGeneric(name = "mergeAnchors", def = function(dlo, mergegap, 
    selfloops = FALSE) standardGeneric("mergeAnchors"))

.mergeAnchors <- function(dlo, mergegap, selfloops) {
    # Join the anchors
    newAnchors <- reduce(dlo@anchors, min.gapwidth = mergegap)
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
    intsdf <- as.data.frame(dlo@loops)
    
    # Update loops indices
    leftmatch <- t(sapply(intsdf$left, function(x) mapping[mapping[, 
        1] == x, ]))
    rightmatch <- t(sapply(intsdf$right, function(x) mapping[mapping[, 
        1] == x, ]))
    totalupdate <- cbind(unlist(leftmatch[, 2]), unlist(rightmatch[, 
        2]))
    upints <- as.matrix(totalupdate)
    colnames(upints) <- c("left", "right")
    
    # Link counts and loops
    df <- data.frame(upints, dlo@counts)
    dnames <- colnames(df)
    dM <- melt(df, id.vars = c("left", "right"))
    print(dM)
    updatedLink <- suppressWarnings(dcast(dM, left + right ~ 
        variable, sum))
    intz <- matrix(c(updatedLink$left, updatedLink$right), ncol = 2)
    colnames(intz) <- c("left", "right")
    
    countz <- data.matrix(updatedLink[, -1:-2])
    mergedObject <- loopdata(anchors = newAnchors, loops = intz, 
        counts = countz)
    
    if (selfloops) {
        return(mergedObject)
    } else {
        return(subsetLoops(mergedObject, mergedObject@loops[, 
            1] != mergedObject@loops[, 2]))
    }
}

#' @rdname mergeAnchors
setMethod(f = "mergeAnchors", signature = c("loopdata", "numeric", 
    "missing"), definition = function(dlo, mergegap, selfloops) {
    .mergeAnchors(dlo, mergegap, FALSE)
})

#' @rdname mergeAnchors
setMethod(f = "mergeAnchors", signature = c("loopdata", "numeric", 
    "logical"), definition = function(dlo, mergegap, selfloops) {
    .mergeAnchors(dlo, mergegap, selfloops)
})

#' Extract region from loopdata object
#'
#' \code{subsetRegion} takes a \code{loopdata} object and a \code{GRanges}
#' object and returns a \code{loopdata} object where both anchors map inside
#' the \code{GRanges} coordinates
#'
#' This function returns a data.frame where the left and right anchors are 
#' visualized together along with the loop width, individual counts, and
#' any anchor meta-data that has been annotated into the anchors \code{GRanges}
#' object.
#'
#' @param dlo A loopdata object to be summarized
#' @param region A GRanges object containing region of interest 
#' 
#' @return A loopdata object
#'
#' @examples
#' # Grab region chr1:36000000-36100000
#' library(GenomicRanges)
#' regA <- GRanges(c('1'),IRanges(c(36000000),c(36100000)))
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_chr1reg.small <- subsetRegion(jpn_chr1reg, regA)
#' @import GenomicRanges
#' @export
setGeneric(name = "subsetRegion", def = function(dlo, region) standardGeneric("subsetRegion"))

#' @rdname subsetRegion
setMethod(f = "subsetRegion", signature = c("loopdata", "GRanges"), 
    definition = function(dlo, region) {
        # Keep only those anchors that are being used
        newAnchors <- dlo@anchors[findOverlaps(region, dlo@anchors)@to, 
            ]
        
        # Create mapping from old indices to new indices
        mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
        intsdf <- as.data.frame(dlo@loops)
        
        # Update loops indices
        leftmatch <- t(sapply(intsdf$left, function(x) mapping[mapping[, 
            1] == x, ]))
        rightmatch <- t(sapply(intsdf$right, function(x) mapping[mapping[, 
            1] == x, ]))
        lm <- suppressWarnings(as.numeric(as.character(leftmatch[, 
            2])))
        rm <- suppressWarnings(as.numeric(as.character(rightmatch[, 
            2])))
        
        # Format new indices matrix
        totalupdate <- cbind(unlist(lm), unlist(rm))
        cc <- complete.cases(totalupdate)
        newloops <- matrix(totalupdate[cc], ncol = 2)
        colnames(newloops) <- c("left", "right")
        
        # Grab counts indicies; removes lines that don't map to
        # anything via making them NAs and then removing them
        newcounts <- matrix(dlo@counts[cc], ncol = ncol(dlo@counts))
        colnames(newcounts) <- colnames(dlo@counts)
        
        # Update values
        slot(dlo, "anchors", check = TRUE) <- newAnchors
        slot(dlo, "loops", check = TRUE) <- newloops
        slot(dlo, "counts", check = TRUE) <- newcounts
        return(dlo)
    })

#' Get number of anchors in each sample
#'
#' \code{numAnchors} takes a \code{loopdata} object and a summarizes the
#' number of anchors that support all the loops (count >= 1) in the object
#'
#' This function returns a data.frame where the column names specify the
#' sample in the original \code{loopdata} object and the only row shows
#' the number of anchors used to support that sample
#'
#' @param x A loopdata object to be summarized
#' 
#' @return A data.frame of each sample and the number of anchors
#'
#' @examples
#' # Show number of anchors each sample is supported by
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' numAnchors(jpn_chr1reg)
#' @export
setGeneric(name = "numAnchors", def = function(x) standardGeneric("numAnchors"))

#' @rdname numAnchors
setMethod(f = "numAnchors", signature = c("loopdata"), definition = function(x) {
    nAnchors <- sapply(1:as.numeric(dim(x)[3]), function(t) {
        length(unique(stack(as.data.frame(x@loops[x@counts[, 
            t] == "0", ])))$values)
    })
    nAnchors <- as.data.frame(t(nAnchors))
    colnames(nAnchors) <- colnames(x@counts)
    return(nAnchors)
})