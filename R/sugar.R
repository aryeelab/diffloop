#' @include loopFunctions.R
NULL

#' Extract parts of a loopdata object
#' 
#' @param x A loopdata object for subsetting
#' @param i Loops to be subsetted
#' @param j Samples to be subsetted 
#' @param drop Other non-essential parameters needed for sub
#' 
#' @return A loopdata object
#' 
setMethod("[", signature(x = "loopdata", i = "numeric", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
    upints <- x@loops[i, ]
    slot(x, "loops", check = TRUE) <- upints
    ucounts <- x@counts[i, ]
    slot(x, "counts", check = TRUE) <- ucounts
    x <- cleanup(x)
    
    slot(x, "counts", check = TRUE) <- x@counts[, j]
    nonZero <- apply(x@counts, MARGIN = 1, function(t) !all(t == 
        0))
    x <- subsetLoops(x, nonZero)
    slot(x, "colData", check = TRUE) <- x@colData[j, ]
    return(x)
})

#' @rdname sub-loopdata-numeric-numeric-missing-method 
setMethod("[", signature(x = "loopdata", i = "missing", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
    slot(x, "counts", check = TRUE) <- x@counts[, j]
    nonZero <- apply(x@counts, MARGIN = 1, function(t) !all(t == 
        0))
    x <- subsetLoops(x, nonZero)
    slot(x, "colData", check = TRUE) <- x@colData[j, ]
    return(x)
})

#' @rdname sub-loopdata-numeric-numeric-missing-method 
setMethod("[", signature(x = "loopdata", i = "numeric", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
    upints <- x@loops[i, ]
    slot(x, "loops", check = TRUE) <- upints
    ucounts <- x@counts[i, ]
    slot(x, "counts", check = TRUE) <- ucounts
    return(cleanup(x))
})

#' Extract parts of looptest
#' 
#' Can only subsample loops, so only the i parameter can be handled
#'
#' @param x A looptest object for subsetting
#' @param i Loops to be subsetted
#' @param j NULL 
#' @param drop NULL
#' 
#' @return A looptest object
#' 
setMethod("[", signature(x = "looptest", i = "numeric", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
    slot(x, "loopdata", check = TRUE) <- x@loopdata[i, ]
    slot(x, "results", check = TRUE) <- x@results[i, ]
    return(x)
})

#' Extract first part of loopdata
#'
#' @param x A looptest object
#' @param n Number of lines to view
#' @param ... Other non-essential params
#' 
#' @return A loopdata object
#' 
setMethod("head", signature = "loopdata", function(x, n = 6, 
    ...) {
    return(loopdata(anchors = head(x@anchors, n), loops = head(x@loops, 
        n), counts = head(x@counts, n), colData = x@colData))
})

#' Extract first part of looptest
#' 
#' @param x A looptest object
#' @param n Number of lines to view
#' @param ... Other non-essential params
#' 
#' @return A looptest object
#' 
setMethod("head", signature = "looptest", function(x, n = 6, 
    ...) {
    return(looptest(loopdata = head(x@loopdata), results = head(x@results)))
})

#' Extract last part of loopdata
#' 
#' @param x A loopdata object
#' @param n Number of lines to view
#' @param ... Other non-essential params
#' 
#' @return A loopdata object
#' 
setMethod("tail", signature = "loopdata", function(x, n = 6, 
    ...) {
    return(loopdata(anchors = tail(x@anchors, n), loops = tail(x@loops, 
        n), counts = tail(x@counts, n), colData = x@colData))
})

#' Extract last part of looptest
#' 
#' @param x A looptest object
#' @param n Number of lines to view
#' @param ... Other non-essential params
#' 
#' @return A looptest object
#'
setMethod("tail", signature = "looptest", function(x, n = 6, 
    ...) {
    return(looptest(loopdata = tail(x@loopdata), results = tail(x@results)))
})

#' See dimensions of loopdata
#' 
#' @param x A loopdata object
#' 
#' @return A data.frame of dimensions of the loopdata object,
#' including number of anchors, loops, samples, and column data
#' attributes 
#'
setMethod("dim", signature = "loopdata", function(x) {
    anchors <- length(x@anchors)
    loops <- nrow(x@loops)
    samples <- ncol(x@counts)
    colData <- ncol(x@colData)
    return(data.frame(cbind(anchors, loops, samples, colData)))
})

#' See dimensions of looptest
#' 
#' @param x A looptest object
#' 
#' @return A data.frame of dimensons of the looptest object,
#' including number of anchors, loops, samples, column data
#' attributes, and columns in the results subobject
#'
setMethod("dim", signature = "looptest", function(x) {
    anchors <- length(x@loopdata@anchors)
    loops <- nrow(x@loopdata@loops)
    samples <- ncol(x@loopdata@counts)
    colData <- ncol(x@loopdata@colData)
    results <- dim(x@results)[2]
    return(data.frame(cbind(anchors, loops, samples, colData, 
        results)))
})

#' Remove 'chr' from GRanges seqnames
#'
#' \code{rmchr} takes a loopdata object or GRanges object and 
#' simply removes the 'chr' from seqnames, if is present
#'
#' Often times, performing functions on GRanges objects can go awry 
#' if the seqnames are systematically different. A common example of
#' this is when some GRanges objects has the format of 'chr1' 
#' while the other has '1'. We can remove 'chr' from the first object 
#'
#' @param dlo A loopdata object or GRanges object
#'
#' @return An identical loopdata/GRanges object except 'chr' removed
#'
#' @examples
#' library(GenomicRanges)
#' regA <- GRanges(c('1'),IRanges(c(36200000),c(36300000)))
#' addchr(regA)
#' regA
#' rmchr(regA)
#' regA

#' @import GenomicRanges
#' 
#' @export
setGeneric(name = "rmchr", def = function(dlo) standardGeneric("rmchr"))

#' @rdname rmchr
setMethod(f = "rmchr", signature = c("loopdata"), definition = function(dlo) {
    seqlevels(dlo@anchors) <- gsub("^chr(.*)$", "\\1", seqlevels(dlo@anchors))
    return(dlo)
})

#' @rdname rmchr
setMethod(f = "rmchr", signature = c("GRanges"), definition = function(dlo) {
    seqlevels(dlo) <- gsub("^chr(.*)$", "\\1", seqlevels(dlo))
    return(dlo)
})

#' Add 'chr' to GRanges seqnames
#'
#' \code{addchr} takes a loopdata object or GRanges object and 
#' simply adds 'chr' to seqnames
#'
#' Often times, performing functions on GRanges objects can go awry 
#' if the seqnames are systematically different. A common example of
#' this is when some GRanges objects has the format of 'chr1' 
#' while the other has '1'. We can add 'chr' to the first object 
#'
#' @param dlo A loopdata object or GRanges object
#'
#' @return An identical loopdata object or GRanges object 'chr' added 
#'
#' @examples
#' library(GenomicRanges)
#' regA <- GRanges(c('1'),ranges=IRanges(c(36200000),c(36300000)))
#' addchr(regA)
#' regA
#' rmchr(regA)
#' regA
#' @import GenomicRanges
#' @export
setGeneric(name = "addchr", def = function(dlo) standardGeneric("addchr"))

#' @rdname addchr
setMethod(f = "addchr", signature = c("loopdata"), definition = function(dlo) {
    seqlevels(dlo@anchors) <- gsub("^(.*)$", "chr\\1", seqlevels(dlo@anchors))
    return(dlo)
})

#' @rdname addchr
setMethod(f = "addchr", signature = c("GRanges"), definition = function(dlo) {
    seqlevels(dlo) <- gsub("^(.*)$", "chr\\1", seqlevels(dlo))
    return(dlo)
})

#' Read a file and make a GRanges object
#'
#' \code{bedToGRanges} takes a string corresponding to a file 
#' and creates a GRanges object, retaining meta-data
#'
#' Useful function to read in a .bed file to create a GRanges object
#' where the meta-data is presevered. Useful for later functions like
#' \code{annotateAnchors} 
#'
#' @param file A string specifiying .bed file location
#'
#' @return A GRanges object
#'
#' @examples
#' #Read in CTCF Jurkat peaks in
#' ctcf_j <- system.file('extdata','Jurkat_CTCF_chr1.narrowPeak',package = 'diffloop')
#' ctcf <- bedToGRanges(ctcf_j)
#' 
#' 
#' @import GenomicRanges
#' @importFrom utils read.table stack
#' @export
setGeneric(name = "bedToGRanges", def = function(file) standardGeneric("bedToGRanges"))

#' @rdname bedToGRanges 
setMethod(f = "bedToGRanges", signature = c("character"), definition = function(file) {
    df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
    if (length(df) < 3) 
        stop("Not a valid bed file- < three columns")
    
    if (length(df) > 6) 
        df <- df[, -c(7:length(df))]
    header <- c("chr", "start", "end", "id", "score", "strand")
    names(df) <- header[1:length(names(df))]
    
    if ("strand" %in% colnames(df)) 
        df$strand <- gsub(pattern = "[^+-]+", replacement = "*", 
            x = df$strand)
    
    # check to see if chromosome names are preceded with 'chr'
    if (any(apply(df, 2, nchar)[, 1] < 4)) {
        df$chr <- sub("^", "chr", df$chr)
    }
    if (length(df) == 3) {
        gr <- with(df, GRanges(chr, IRanges(start, end)))
    } else if (length(df) == 4) {
        gr <- with(df, GRanges(chr, IRanges(start, end), id = id))
    } else if (length(df) == 5) {
        gr <- with(df, GRanges(chr, IRanges(start, end), id = id, 
            score = score))
    } else {
        gr <- with(df, GRanges(chr, IRanges(start, end), id = id, 
            score = score, strand = strand))
    }
    return(gr)
})

#' Pad a GRanges object
#'
#' \code{padGRanges} takes a GRanges object and adds or substracts 
#' distance based on user-defined input. Upstream and downstream
#' consider strand information when available. Specify only either
#' pad or upstream/downstream when using
#'
#' @param gro A granges object
#' @param upstream Distance in BP added upstream
#' @param downstream Distance in BP added downstream
#' @param pad Distance in BP added 
#'
#' @return A GRanges object with adjusted start and end values
#'
#' @examples
#' #Read in CTCF Jurkat peaks in
#' ctcf_j <- system.file('extdata','Jurkat_CTCF_chr1.narrowPeak',package = 'diffloop')
#' ctcf <- bedToGRanges(ctcf_j)
#' ctcf.pad <- padGRanges(ctcf, pad = 1000)
#' 
#' 
#' @import GenomicRanges
#' @export
setGeneric(name = "padGRanges", def = function(gro, upstream = 0, 
    downstream = 0, pad = 0) standardGeneric("padGRanges"))

#' @rdname padGRanges 
setMethod(f = "padGRanges", signature = c("GRanges", "ANY", "ANY", 
    "ANY"), definition = function(gro, upstream = 0, downstream = 0, 
    pad = 0) {
    g.df <- data.frame(gro)
    g.df[g.df$strand == "-", ]$start <- g.df[g.df$strand == "-", 
        ]$start - pad - downstream
    g.df[g.df$strand == "-", ]$end <- g.df[g.df$strand == "-", 
        ]$end + pad + upstream
    pos <- g.df$strand == "+" | g.df$strand == "*"
    g.df[pos, ]$start <- g.df[pos, ]$start - pad - upstream
    g.df[pos, ]$end <- g.df[pos, ]$end + pad + downstream
    return(GRanges(g.df))
})
