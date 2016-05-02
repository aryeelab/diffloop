#' @include data.R
NULL

#' Link the anchors and interactions back together 
#'
#' \code{summary} takes a \code{loops} object and breaks the 
#' loop data structure resulting in a \code{data.frame}. 
#'
#' This function returns a \code{data.frame} where the left and right anchors 
#' are visualized together along with the loop width, individual counts, and
#' any anchor meta-data that has been annotated into the anchors GRanges
#' object as well as any rowData varianble
#'
#' @param object A loops object to be summarized
#'
#' @return A data.frame
#'
#' @examples
#' # Summarizing the first ten loops in \code{loops.small}
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' summarydf <- summary(loops.small[1:10,])
#' # Summarizing the loops and significance results between naive and primed
#' summarylt <- summary(quickAssoc(loops.small[,1:4])[1:10,])
#' @import plyr

#' @export
setMethod(f = "summary", signature = c("loops"), definition = function(object) {
    dlo <- object
    
    # Grab all the left anchors in order of loop occurence
    leftAnchor2 <- as.data.frame(dlo@anchors[dlo@interactions[, 
        1]])
    leftAnchor2 <- subset(leftAnchor2, select = -c(width, strand))
    colnames(leftAnchor2) <- paste(colnames(leftAnchor2), "1", 
        sep = "_")
    colnames(leftAnchor2)[1] <- "chr_1"
    
    # Grab all the right anchors in order of loop occurence
    rightAnchor2 <- as.data.frame(dlo@anchors[dlo@interactions[, 
        2]])
    rightAnchor2 <- subset(rightAnchor2, select = -c(width, strand))
    colnames(rightAnchor2) <- paste(colnames(rightAnchor2), "2", 
        sep = "_")
    colnames(rightAnchor2)[1] <- "chr_2"
    
    # Add the loop features
    cbind(leftAnchor2, rightAnchor2, dlo@counts, object@rowData)
})

# Function that removes all anchors not being referenced in
# interactions matrix and updates indices. For internal use
# only.
setGeneric(name = "cleanup", def = function(dlo) standardGeneric("cleanup"))
setMethod(f = "cleanup", signature = c("loops"), definition = function(dlo) {
    if (dim(dlo@counts)[1] == 0) {
        stop("Attempting to subset to empty loops object!")
    }
    
    # Grab indicies of anchors being referenced in loops
    idf <- data.frame(dlo@interactions[, 1], dlo@interactions[, 
        2])
    sdf <- stack(idf)
    udf <- sort(unique(sdf[, "values"]))
    
    # Keep only those anchors that are being used
    newAnchors <- dlo@anchors[udf]
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
    intsdf <- as.data.frame(dlo@interactions)
    
    # Update interactions indices
    leftmatch <- t(sapply(intsdf$left, function(x) mapping[mapping[, 
        1] == x, ]))
    rightmatch <- t(sapply(intsdf$right, function(x) mapping[mapping[, 
        1] == x, ]))
    
    # Format new indices matrix
    totalupdate <- cbind(unlist(leftmatch[, 2]), unlist(rightmatch[, 
        2]))
    upints <- as.matrix(totalupdate)
    colnames(upints) <- c("left", "right")
    rownames(upints) <- NULL
    
    # Update values
    slot(dlo, "interactions", check = TRUE) <- upints
    slot(dlo, "anchors", check = TRUE) <- newAnchors
    return(dlo)
})

#' Combine nearby anchors into one peak
#'
#' \code{mergeAnchors} combines anchors that are within a user-defined radius
#'
#' This function takes a loops object and combines nearby anchors, up to
#' a distance specified by the \code{mergegap}. This likely will cause self
#' loops to form (loop where the left and right anchor are the same), which
#' can either be removed (by default) or retained with \code{selfloops}
#'
#' @param dlo A loops object whose anchors will be merged
#' @param mergegap An integer value of the bp between anchors to be merged
#' @param selfloops A logical value to either retain (T) or remove (F) 
#' resulting self-loops after merging anchors
#'
#' @return A loops object
#'
#' @examples
#' # Merge anchors within 1kb of each other, keeping self loops 
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' m1kb <- mergeAnchors(loops.small, 1000, FALSE)
#'
#' # Merge anchors within 1kb of each other, removing self loops by default
#' m1kb_unique <- mergeAnchors(loops.small, 1000)

#' @import GenomicRanges
#' @import reshape2
#' @export
#' 
setGeneric(name = "mergeAnchors", def = function(dlo, mergegap, 
    selfloops = FALSE) standardGeneric("mergeAnchors"))

.mergeAnchors <- function(dlo, mergegap, selfloops) {
    # Join the anchors
    newAnchors <- reduce(dlo@anchors, min.gapwidth = mergegap)
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
    intsdf <- as.data.frame(dlo@interactions)
    
    # Update interactions indices
    leftmatch <- t(sapply(intsdf$left, function(x) mapping[mapping[, 
        1] == x, ]))
    rightmatch <- t(sapply(intsdf$right, function(x) mapping[mapping[, 
        1] == x, ]))
    totalupdate <- cbind(unlist(leftmatch[, 2]), unlist(rightmatch[, 
        2]))
    upints <- as.matrix(totalupdate)
    colnames(upints) <- c("left", "right")
    
    # Link counts and interactions
    df <- data.frame(upints, dlo@counts)
    dnames <- colnames(df)
    dM <- melt(df, id.vars = c("left", "right"))
    print(dM)
    updatedLink <- suppressWarnings(dcast(dM, left + right ~ 
        variable, sum))
    intz <- matrix(c(updatedLink$left, updatedLink$right), ncol = 2)
    colnames(intz) <- c("left", "right")
    
    countz <- data.matrix(updatedLink[, -1:-2])
    mergedObject <- loops(anchors = newAnchors, interactions = intz, 
        counts = countz)
    
    if (selfloops) {
        return(mergedObject)
    } else {
        return(subsetLoops(mergedObject, mergedObject@interactions[, 
            1] != mergedObject@interactions[, 2]))
    }
}

#' @rdname mergeAnchors
setMethod(f = "mergeAnchors", signature = c("loops", "numeric", 
    "missing"), definition = function(dlo, mergegap, selfloops) {
    .mergeAnchors(dlo, mergegap, FALSE)
})

#' @rdname mergeAnchors
setMethod(f = "mergeAnchors", signature = c("loops", "numeric", 
    "logical"), definition = function(dlo, mergegap, selfloops) {
    .mergeAnchors(dlo, mergegap, selfloops)
})


#' Extract region from loops object
#'
#' \code{subsetRegion} takes a \code{loops} object and a \code{GRanges}
#' object and returns a \code{loops} object where both anchors map inside
#' the \code{GRanges} coordinates by default. Once can specify where only
#' one anchor is in the region instead.
#'
#' By default, \code{nachors = 2}, meaning both anchors need to be in the 
#' region for the loop to be preserved when extracting. However, by specifying
#' a numeric 1, interactions with either the left or right anchor will be extracted.
#' Loops with both anchors in the region will be excluded (exclusive \code{or}).
#' To get an inclusive \code{or}, take the union of subsetting both with 1 and 2.
#'
#' @param dlo A loops object to be subsetted
#' @param region A GRanges object containing region of interest 
#' @param nanchors Number of anchors to be contained in GRanges object. Default 2
#' 
#' @return A loops object
#'
#' @examples
#' # Grab region chr1:36000000-36100000
#' library(GenomicRanges)
#' regA <- GRanges(c('1'),IRanges(c(36000000),c(36100000)))
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' # Both anchors in region
#' loops.small.two <- subsetRegion(loops.small, regA)
#' #Only one anchor in region
#' loops.small.one <- subsetRegion(loops.small, regA, 1)
#' #Either one or two anchors in region
#' loops.small.both <- union(loops.small.one, loops.small.two)
#' @import GenomicRanges
#' @export
setGeneric(name = "subsetRegion", def = function(dlo, region, 
    nanchors) standardGeneric("subsetRegion"))

#' @rdname subsetRegion
setMethod(f = "subsetRegion", signature = c("loops", "GRanges", 
    "numeric"), definition = function(dlo, region, nanchors) {
    if (nanchors == "1") {
        .subsetRegion1(dlo, region)
    } else if (nanchors == "2") {
        .subsetRegion2(dlo, region)
    } else {
        print("Please specify either 1 or 2 anchors in region")
    }
})

#' @rdname subsetRegion
setMethod(f = "subsetRegion", signature = c("loops", "GRanges", 
    "missing"), definition = function(dlo, region, nanchors) {
    .subsetRegion2(dlo, region)
})

.subsetRegion1 <- function(dlo, region) {
    # Keep only those anchors that are being used
    newAnchors <- dlo@anchors[findOverlaps(region, dlo@anchors)@to, 
        ]
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
    intsdf <- as.data.frame(dlo@interactions)
    
    # Update interactions indices
    leftmatch <- t(sapply(intsdf$left, function(x) mapping[mapping[, 
        1] == x, ]))
    rightmatch <- t(sapply(intsdf$right, function(x) mapping[mapping[, 
        1] == x, ]))
    lm <- suppressWarnings(as.numeric(as.character(leftmatch[, 
        2])))
    rm <- suppressWarnings(as.numeric(as.character(rightmatch[, 
        2])))
    
    keepTheseLoops <- xor(!is.na(lm), !is.na(rm))  # only one anchor in region
    keepTheseAnchors <- unique(unlist(intsdf[keepTheseLoops, 
        ], use.names = FALSE))
    
    # Go back through the motions
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, dlo@anchors[keepTheseAnchors, 
        ]))
    intsdf <- as.data.frame(dlo@interactions)
    
    # Update interactions indices
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
    newinteractions <- matrix(totalupdate[cc], ncol = 2)
    colnames(newinteractions) <- c("left", "right")
    
    # Subset rowData
    newRowData <- as.data.frame(dlo@rowData[cc, ])
    colnames(newRowData) <- colnames(dlo@rowData)
    
    # Grab counts indicies; removes lines that don't map to
    # anything via making them NAs and then removing them
    newcounts <- matrix(dlo@counts[cc], ncol = ncol(dlo@counts))
    colnames(newcounts) <- colnames(dlo@counts)
    
    # Update values
    slot(dlo, "anchors", check = TRUE) <- dlo@anchors[keepTheseAnchors, 
        ]
    slot(dlo, "interactions", check = TRUE) <- newinteractions
    slot(dlo, "counts", check = TRUE) <- newcounts
    slot(dlo, "rowData", check = TRUE) <- newRowData
    return(dlo)
}

.subsetRegion2 <- function(dlo, region) {
    # Keep only those anchors that are being used
    newAnchors <- dlo@anchors[findOverlaps(region, dlo@anchors)@to, 
        ]
    
    # Create mapping from old indices to new indices
    mapping <- as.data.frame(findOverlaps(dlo@anchors, newAnchors))
    intsdf <- as.data.frame(dlo@interactions)
    
    # Update interactions indices
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
    newinteractions <- matrix(totalupdate[cc], ncol = 2)
    colnames(newinteractions) <- c("left", "right")
    
    # Grab counts indicies; removes lines that don't map to
    # anything via making them NAs and then removing them
    newcounts <- matrix(dlo@counts[cc], ncol = ncol(dlo@counts))
    colnames(newcounts) <- colnames(dlo@counts)
    
    # Subset rowData
    newRowData <- as.data.frame(dlo@rowData[cc, ])
    colnames(newRowData) <- colnames(dlo@rowData)
    
    # Update values
    slot(dlo, "anchors", check = TRUE) <- newAnchors
    slot(dlo, "interactions", check = TRUE) <- newinteractions
    slot(dlo, "counts", check = TRUE) <- newcounts
    slot(dlo, "rowData", check = TRUE) <- newRowData
    return(dlo)
}

#' Get number of anchors in each sample
#'
#' \code{numAnchors} takes a \code{loops} object and a summarizes the
#' number of anchors that support all the interactions (count >= 1) in the object
#'
#' This function returns a data.frame where the column names specify the
#' sample in the original \code{loops} object and the only row shows
#' the number of anchors used to support that sample
#'
#' @param x A loops object to be summarized
#' 
#' @return A data.frame of each sample and the number of anchors
#'
#' @examples
#' # Show number of anchors each sample is supported by
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' numAnchors(loops.small)
#' @export
setGeneric(name = "numAnchors", def = function(x) standardGeneric("numAnchors"))

#' @rdname numAnchors
setMethod(f = "numAnchors", signature = c("loops"), definition = function(x) {
    nAnchors <- sapply(1:as.numeric(dim(x)[3]), function(t) {
        length(unique(stack(as.data.frame(x@interactions[x@counts[, 
            t] == "0", ])))$values)
    })
    nAnchors <- as.data.frame(t(nAnchors))
    colnames(nAnchors) <- colnames(x@counts)
    return(nAnchors)
})