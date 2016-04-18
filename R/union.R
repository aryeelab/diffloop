#' @include sugar.R
NULL

#' Combine two loopdata objects
#'
#' \code{union} combines two loopdata objects' loops and anchors and 
#' populates the \code{colData} matrix where available
#'
#' This function returns a single loopdata object that has all the 
#' anchors and loops contained in the two loopdata objects that were
#' part of the input. However, when the two objects have different samples,
#' the counts matrix will contain missing values (e.g. when loop counts
#' in x are not in y, those values are unknown). While the number of loops,
#' colData, and anchors should be correct, we need to correct the counts
#' using a subsetting function.
#'
#' @param x A loopdata object 
#' @param y A loopdata object
#'
#' @return A loopdata obect
#'
#' @examples
#' # divide and recombine samples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' naive <- jpn_chr1reg[,1:2]
#' primed <- jpn_chr1reg[,3:4]
#' np <- union(naive, primed)
#' # Subset from full to get correct counts
#' c.np <- loopdataSubset(np, jpn_chr1reg)

#' @import plyr
#' @export
setMethod(f = "union", signature = c("loopdata", "loopdata"), 
    definition = function(x, y) {
        o1df <- summarize(x)
        o2df <- summarize(y)
        
        # Get combined anchors
        a1 <- data.frame(o1df$chr_1, o1df$start_1, o1df$end_1)
        a2 <- data.frame(o2df$chr_1, o2df$start_1, o2df$end_1)
        a3 <- data.frame(o1df$chr_2, o1df$start_2, o1df$end_2)
        a4 <- data.frame(o2df$chr_2, o2df$start_2, o2df$end_2)
        
        colnames(a1) <- NULL
        colnames(a2) <- NULL
        colnames(a3) <- NULL
        colnames(a4) <- NULL
        
        la <- list(a1, a2, a3, a4)
        
        a0 <- ldply(la, data.frame)
        anchors <- makeGRangesFromDataFrame(a0, seqnames.field = "X1", 
            start.field = "X2", end.field = "X3")
        anchors <- reduce(anchors)
        
        # Index loops in df1
        int1 <- apply(data.frame(a1, a3), 1, function(t) {
            i1 <- which(as.list(seqnames(anchors)) == t[1] & 
                as.integer(start(ranges(anchors))) == as.integer(t[2]) & 
                as.integer(end(ranges(anchors))) == as.integer(t[3]))
            i2 <- which(as.list(seqnames(anchors)) == t[4] & 
                as.integer(start(ranges(anchors))) == as.integer(t[5]) & 
                as.integer(end(ranges(anchors))) == as.integer(t[6]))
            cbind(i1, i2)
        })
        
        # Index loops in df2
        int2 <- apply(data.frame(a2, a4), 1, function(t) {
            i1 <- which(as.list(seqnames(anchors)) == t[1] & 
                as.integer(start(ranges(anchors))) == as.integer(t[2]) & 
                as.integer(end(ranges(anchors))) == as.integer(t[3]))
            i2 <- which(as.list(seqnames(anchors)) == t[4] & 
                as.integer(start(ranges(anchors))) == as.integer(t[5]) & 
                as.integer(end(ranges(anchors))) == as.integer(t[6]))
            cbind(i1, i2)
        })
        
        in1 <- t(int1)
        in2 <- t(int2)
        
        colnames(in1) <- c("left", "right")
        colnames(in2) <- c("left", "right")
        
        # Get counts
        sam1 <- colnames(x@counts)
        sam2 <- colnames(y@counts)
        cs1 <- sapply(sam1, function(s) {
            indx <- grep(s, colnames(o1df))
            o1df[, indx]
        })
        cs2 <- sapply(sam2, function(s) {
            indx <- grep(s, colnames(o2df))
            o2df[, indx]
        })
        
        d1 <- melt(data.frame(in1, cs1), id.vars = c("left", 
            "right"))
        d2 <- melt(data.frame(in2, cs2), id.vars = c("left", 
            "right"))
        
        # Grab unique rows only, then reshape
        base <- unique(ldply(list(d1, d2), data.frame))
        bigTab <- suppressWarnings(dcast(base, left + right ~ 
            variable, max))
        newloops <- matrix(c(bigTab$left, bigTab$right), ncol = 2)
        colnames(newloops) <- c("left", "right")
        newCounts <- data.matrix(bigTab[, -1:-2])
        
        # Update colData
        
        unsorted <- rbind(x@colData, y@colData)
        newcolData <- unsorted[match(colnames(newCounts), rownames(unsorted)), 
            ]
        
        cat("Check for NAs; Subset this object from  more comprehensive object\n")
        return(loopdata(anchors = anchors, loops = newloops, 
            counts = newCounts, colData = newcolData))
    })


#' Subset two difloop objects
#'
#' \code{loopdataSubset} takes the loops and anchors present in dlo1
#' and uses the counts and samples from dlo2. 
#'
#' This function plays nice with \code{union} to ensure counts are correct
#' after taking the union of two loopdata objects. The subset function simply
#' returns the anchors and loops of dlo1 and the counts and colData of dlo2.
#'
#' @param dlo1 A loopdata object 
#' @param dlo2 A loopdata object
#'
#' @return A loopdata obect
#'
#' @examples
#' # divide and recombine samples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' naive <- jpn_chr1reg[,1:2]
#' primed <- jpn_chr1reg[,3:4]
#' np <- union(naive, primed)
#' # Subset from full to get correct counts
#' c.np <- loopdataSubset(np, jpn_chr1reg)

#' @export
setGeneric(name = "loopdataSubset", def = function(dlo1, dlo2) standardGeneric("loopdataSubset"))

#' @rdname loopdataSubset
setMethod(f = "loopdataSubset", signature = c("loopdata", "loopdata"), 
    definition = function(dlo1, dlo2) {
        va <- as.integer(findOverlaps(dlo1@anchors, dlo2@anchors)@from)
        idd <- apply(dlo2@loops, 1, function(t) {
            ta <- (is.element(as.integer(t[1]), va))
            tb <- (is.element(as.integer(t[2]), va))
            ta & tb
        })
        return(subsetLoops(dlo2, idd))
    })
