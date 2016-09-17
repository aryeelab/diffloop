#' @include loopFunctions.R
NULL

#' Extract parts of a loops object
#' 
#' @param x A loops object for subsetting
#' @param i Loops to be subsetted
#' @param j Samples to be subsetted 
#' @param drop Other non-essential parameters needed for sub
#' 
#' @return A loops object
#' 
setMethod("[", signature(x = "loops", i = "numeric", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
    upints <- x@interactions[i, ]
    slot(x, "interactions", check = TRUE) <- upints
    ucounts <- x@counts[i, ]
    slot(x, "counts", check = TRUE) <- ucounts
    newRowData <- as.data.frame(x@rowData[i, ])
    colnames(newRowData) <- colnames(x@rowData)
    slot(x, "rowData", check = TRUE) <- newRowData
    x <- cleanup(x)
    
    slot(x, "counts", check = TRUE) <- x@counts[, j]
    nonZero <- apply(x@counts, MARGIN = 1, function(t) !all(t == 
        0))
    x <- subsetLoops(x, nonZero)
    slot(x, "colData", check = TRUE) <- x@colData[j, ]
    return(x)
})

#' @rdname sub-loops-numeric-numeric-missing-method 
setMethod("[", signature(x = "loops", i = "missing", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
    ncounts <- as.matrix(x@counts[, j], ncol = j)
    colnames(ncounts) <- colnames(x@counts)[j]
    slot(x, "counts", check = TRUE) <- ncounts
    nonZero <- apply(x@counts, MARGIN = 1, function(t) !all(t == 0))
    x <- subsetLoops(x, nonZero)
    slot(x, "colData", check = TRUE) <- x@colData[j, ]
    return(x)
})

#' @rdname sub-loops-numeric-numeric-missing-method 
setMethod("[", signature(x = "loops", i = "numeric", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
    upints <- x@interactions[i, ]
    slot(x, "interactions", check = TRUE) <- as.matrix(upints)
    ucounts <- x@counts[i, ]
    slot(x, "counts", check = TRUE) <- ucounts
    newRowData <- as.data.frame(x@rowData[i, ])
    colnames(newRowData) <- colnames(x@rowData)
    slot(x, "rowData", check = TRUE) <- newRowData
    return(cleanup(x))
})

#' Extract first part of loops object
#'
#' @param x A loops object
#' @param n Number of lines to view
#' @param ... Other non-essential params
#' 
#' @return A loops object
#' 
setMethod("head", signature = "loops", function(x, n = 6, ...) {
    
    anchors <- head(x@anchors, n)
    interactions <- head(x@interactions, n)
    counts <- head(x@counts, n)
    colData <- head(x@colData, n)
    rowData <- head(x@rowData, n)
    
    dlo <- loops()
    slot(dlo, "anchors", check = TRUE) <- anchors
    slot(dlo, "interactions", check = TRUE) <- interactions
    slot(dlo, "counts", check = TRUE) <- counts
    slot(dlo, "colData", check = TRUE) <- colData
    slot(dlo, "rowData", check = TRUE) <- rowData
    return(dlo)
})


#' Extract last part of loops object
#' 
#' @param x A loops object
#' @param n Number of lines to view
#' @param ... Other non-essential params
#' 
#' @return A loops object
#' 
setMethod("tail", signature = "loops", function(x, n = 6, ...) {
    
    anchors <- tail(x@anchors, n)
    interactions <- tail(x@interactions, n)
    counts <- tail(x@counts, n)
    colData <- tail(x@colData, n)
    rowData <- tail(x@rowData, n)
    
    dlo <- loops()
    slot(dlo, "anchors", check = TRUE) <- anchors
    slot(dlo, "interactions", check = TRUE) <- interactions
    slot(dlo, "counts", check = TRUE) <- counts
    slot(dlo, "colData", check = TRUE) <- colData
    slot(dlo, "rowData", check = TRUE) <- rowData
    
    return(dlo)
})

#' See dimensions of loops object
#' 
#' @param x A loops object
#' 
#' @return A data.frame of dimensions of the loops object,
#' including number of anchors, interactions, samples, and column data
#' attributes 
#'
setMethod("dim", signature = "loops", function(x) {
    anchors <- length(x@anchors)
    interactions <- nrow(x@interactions)
    samples <- ncol(x@counts)
    colData <- ncol(x@colData)
    rowData <- ncol(x@rowData)
    return(data.frame(cbind(anchors, interactions, samples, colData, rowData)))
})

#' Remove 'chr' from GRanges seqnames
#'
#' \code{rmchr} takes a loops object or GRanges object and 
#' simply removes the 'chr' from seqnames, if is present
#'
#' Often times, performing functions on GRanges objects can go awry 
#' if the seqnames are systematically different. A common example of
#' this is when some GRanges objects has the format of 'chr1' 
#' while the other has '1'. We can remove 'chr' from the first object 
#'
#' @param dlo A loops object or GRanges object
#'
#' @return An identical loops/GRanges object except 'chr' removed
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
setMethod(f = "rmchr", signature = c("loops"), definition = function(dlo) {
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
#' \code{addchr} takes a loops object or GRanges object and 
#' simply adds 'chr' to seqnames
#'
#' Often times, performing functions on GRanges objects can go awry 
#' if the seqnames are systematically different. A common example of
#' this is when some GRanges objects has the format of 'chr1' 
#' while the other has '1'. We can add 'chr' to the first object 
#'
#' @param dlo A loops object or GRanges object
#'
#' @return An identical loops object or GRanges object 'chr' added 
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
setMethod(f = "addchr", signature = c("loops"), definition = function(dlo) {
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
    g.df[g.df$strand == "-", ]$start <- g.df[g.df$strand == "-", ]$start - pad - downstream
    g.df[g.df$strand == "-", ]$end <- g.df[g.df$strand == "-",  ]$end + pad + upstream
    pos <- g.df$strand == "+" | g.df$strand == "*"
    g.df[pos, ]$start <- g.df[pos, ]$start - pad - upstream
    g.df[pos, ]$end <- g.df[pos, ]$end + pad + downstream
    return(GRanges(g.df))
})
