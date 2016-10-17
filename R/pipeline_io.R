#' @include core.R
NULL

#' Read preprocessed ChIA-PET data from dnaloop
#'
#' \code{loopsMake} reads in a data directory created by the
#' \code{dnaloop} preprocessing pipeline and returns a loops object
#'
#' This function reads in preprocessed ChIA-PET data produced by the
#' \code{dnaloop} preprocessing pipeline. The \code{samples} argument specifies
#' which samples are read. If \code{samples} is not specified all samples will
#' be read. The \code{type} option restricts loops whether they are on the same 
#' 'intra' or different 'inter' chormosome. Default is 'all'. 
#' 
#' IMPORTANT: Assumes the delimiter is a space, not a tab on the files. 
#'
#' @param beddir A string. The preprocessed data directory
#' @param samples A character vector. Optional list of samples to read in
#' @param mergegap An integer value of the radius to merge anchors; default 0
#' @param type Specificies 'intra', 'inter', or 'all' looping. Default 'all'
#' 
#' @return A loops object
#'
#' @examples
#' # Reading in all samples, no mergegap, all loops
#' bd<- system.file('extdata', 'esc_jurkat', package='diffloopdata')
#' # loops <- loopsMake(bd) #standard call
#'
#' # Reading in a subset of samples, 1kb mergegap, only intrachromosomal
#' # looping
#' samples <- c('naive_esc_1', 'naive_esc_2')
#' # naive.intra <- loopsMake(bd, samples, 1000, 'inter')
#'

#' @import foreach
#' @import GenomicRanges
#' @importFrom utils read.table stack
#' @importFrom dplyr group_by full_join
#' @import readr

#' @export
setGeneric(name = "loopsMake", def = function(beddir, samples = NA, 
    mergegap = 0, type = "all") standardGeneric("loopsMake"))

#' @import GenomicRanges
.loopsMake <- function(beddir, samples, mergegap, type) {
    
    ct <- list(col_character(), col_integer(), col_integer(), 
        col_character(), col_integer(), col_integer(), col_character(), 
        col_integer())
    
    restrictPets <- function(bt) {
        if (type == "intra") {
            return(bt[bt$X1 == bt$X4, ])
        } else if (type == "inter") {
            return(bt[bt$X1 != bt$X4, ])
        } else {
            return(bt)
        }
    }
    
    # Iterate through files to set up anchors
    anchorsraw <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "loop_counts.bedpe", sep = "."))
        bt <- read_delim(fullfile, " ", col_types = ct, col_names = FALSE)
        bt <- restrictPets(bt)
        plyr::rbind.fill(bt[, 1:3], setNames(bt[, 4:6], names(bt[, 1:3])))
    }
    
    anchors <- reduce(makeGRangesFromDataFrame(do.call(rbind, 
        anchorsraw), ignore.strand = TRUE, seqnames.field = "X1", 
        start.field = "X2", end.field = "X3"), min.gapwidth = mergegap)
    
    # Map individual reads to anchors
    getpets <- function(left, right, counts, sample) {
        ovl <- findOverlaps(left, anchors)
        leftanchor <- rep(NA, length(left))
        leftanchor[queryHits(ovl)] <- subjectHits(ovl)
        
        ovl <- findOverlaps(right, anchors)
        rightanchor <- rep(NA, length(right))
        rightanchor[queryHits(ovl)] <- subjectHits(ovl)
        
        df <- data.frame(left = leftanchor, right = rightanchor)
        g <- as.data.frame(dplyr::group_by(df, left, right))
        d <- cbind(g, counts)
        colnames(d) <- c("left", "right", "counts")
        dag <- aggregate(counts ~ left + right, FUN = sum, data=d)
        colnames(dag) <- c("left", "right", sample)
        dag
    }
    
    petlist <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "loop_counts.bedpe", 
            sep = "."))
        bt <- read_delim(fullfile, " ", col_types = ct, col_names = FALSE)
        bt <- restrictPets(bt)
        getpets(makeGRangesFromDataFrame(bt[, 1:3], ignore.strand = TRUE, 
            seqnames.field = "X1", start.field = "X2", end.field = "X3"), 
            makeGRangesFromDataFrame(bt[, 4:6], ignore.strand = TRUE, 
                seqnames.field = "X4", start.field = "X5", end.field = "X6"), 
            bt[, 8], sample)
    }
    .full_join <- function(a, b) {
        as.data.frame(dplyr::full_join(a, b, by = c("left", "right")))
    }
    
    # Map Counts
    pets <- Reduce(.full_join, petlist)
    iraw <- pets[, c("left", "right")]
    iraw <- t(apply(iraw, 1, function(x) {
        if (x[1] < x[2]) {
            x
        } else {
            x[c(2, 1)]
        }
    }))
    interactions <- iraw[order(iraw[, 1], iraw[, 2]), ]
    colnames(interactions) <- c("left", "right")
    counts <- as.matrix(pets[, -c(1:2)])[order(iraw[, 1], iraw[, 
        2]), ]
    counts[is.na(counts)] <- 0
    counts <- as.matrix(counts, ncol = length(samples))
    colnames(counts) <- samples
    
    # Initialize rowData slot (with loop widths)
    w <- (start(anchors[interactions[, 2]]) + end(anchors[interactions[, 2]]))/2 -
         (start(anchors[interactions[, 1]]) + end(anchors[interactions[, 1]]))/2
    w[w < 0] <- 0
    rowData <- as.data.frame(as.integer(w))
    colnames(rowData) <- c("loopWidth")
    
    # Remove rownames from matrices
    row.names(interactions) <- NULL
    row.names(counts) <- NULL
    
    # Initialize colData slot
    groups <- rep("group1", length(samples))
    if(length(samples) == 1){
        sizeFactor <- 1
    } else {
        lc <- log2(counts)
        keep <- rowSums(counts > 0) == ncol(lc)
        lc <- lc[keep, ]
        target <- 2^rowMeans(lc)
        sizeFactor <- colMedians(sweep(2^lc, 1, target, FUN = "/"), na.rm = TRUE)
    }
    dfcd <- data.frame(sizeFactor, groups)
    rownames(dfcd) <- samples
    
    # Create loops object
    dlo <- loops()
    slot(dlo, "anchors", check = TRUE) <- anchors
    slot(dlo, "interactions", check = TRUE) <- interactions
    slot(dlo, "counts", check = TRUE) <- counts
    slot(dlo, "colData", check = TRUE) <- dfcd
    slot(dlo, "rowData", check = TRUE) <- rowData
    
    return(dlo)
}

#' @rdname loopsMake
setMethod(f = "loopsMake", def = function(beddir, samples, mergegap = 0, 
    type = "all") {
    if (sum(is.na(samples)) > 0) {
        samples <- dir(beddir, pattern = ".loop_counts.bedpe")
        samples <- sub(".loop_counts.bedpe", "", samples)
    }
    .loopsMake(beddir, samples, mergegap, type)
})

#' Read preprocessed ChIA-PET data from mango
#'
#' \code{loopsMake.mango} reads in a data directory created by the
#' \code{mango} preprocessing pipeline and returns a loops object
#'
#' This function reads in preprocessed ChIA-PET data produced by the
#' \code{mango} preprocessing pipeline. The \code{samples} argument specifies
#' which samples are read. If \code{samples} is not specified all samples will
#' be read. The \code{ext} specifies which type of file to look for, either
#' \code{all} or \code{fdr}, with \code{all} being the default. Under the
#' default, all samples with the extension \code{.fdr.mango} will be processed.
#' Finally, the \code{FDR} parameter (default = 1) specifies the minimum threshold
#' for loops to be added to the greater \code{loops} object. Currently, we do not
#' support importing the verbose output, so the verbose parameter when executing
#' mango should be set to FALSE or the user will have ot parse the file before
#' reading into diffloop using \code{awk}, \code{cut}, or something similar. 
#'
#' @param beddir A string. The preprocessed data directory; Required
#' @param samples A character vector. Optional list of samples to read in
#' @param mergegap An integer value of the radius to merge anchors; default 500
#' @param ext Specificies 'all' or 'fdr' file format; default 'all'
#'
#' @return A loops object where 'chr' is removed from the anchors.
#'
#' @examples
#' # UPDATE THIS
#' bd<- system.file('extdata', 'esc_jurkat', package='diffloopdata')

#' @import foreach
#' @import GenomicRanges
#' @importFrom utils read.table stack
#' @importFrom dplyr group_by full_join
#' @import readr

#' @export
setGeneric(name = "loopsMake.mango", def = function(beddir, samples = NA, 
    mergegap = 500, ext = "all") standardGeneric("loopsMake.mango"))

#' @import GenomicRanges
.loopsMake.mango <- function(beddir, samples, mergegap, ext) {
    
    ct <- list(col_character(), col_integer(), col_integer(), 
        col_character(), col_integer(), col_integer(), col_integer(), 
        col_number())
    
    # Iterate through files to set up anchors
    anchorsraw <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "interactions", ext, "mango", sep = "."))
        bt <- read_delim(fullfile, "\t", col_types = ct, col_names = FALSE)
        plyr::rbind.fill(bt[, 1:3], setNames(bt[, 4:6], names(bt[, 1:3])))
    }
    
    anchors <- reduce(makeGRangesFromDataFrame(do.call(rbind, 
        anchorsraw), ignore.strand = TRUE, seqnames.field = "X1", 
        start.field = "X2", end.field = "X3"), min.gapwidth = mergegap)
    
    # Map individual reads to anchors
    getpets <- function(left, right, counts, sample) {
        ovl <- findOverlaps(left, anchors)
        leftanchor <- rep(NA, length(left))
        leftanchor[queryHits(ovl)] <- subjectHits(ovl)
        
        ovl <- findOverlaps(right, anchors)
        rightanchor <- rep(NA, length(right))
        rightanchor[queryHits(ovl)] <- subjectHits(ovl)
        
        df <- data.frame(left = leftanchor, right = rightanchor)
        g <- as.data.frame(dplyr::group_by(df, left, right))
        d <- cbind(g, counts)
        colnames(d) <- c("left", "right", "counts")
        dag <- aggregate(counts ~ left + right, FUN = sum, data=d)
        colnames(dag) <- c("left", "right", sample)
        dag
    }
    
    petlist <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "interactions", ext, "mango", sep = "."))
        bt <- read_delim(fullfile, "\t", col_types = ct, col_names = FALSE)
        getpets(makeGRangesFromDataFrame(bt[, 1:3], ignore.strand = TRUE, 
            seqnames.field = "X1", start.field = "X2", end.field = "X3"), 
            makeGRangesFromDataFrame(bt[, 4:6], ignore.strand = TRUE, 
                seqnames.field = "X4", start.field = "X5", end.field = "X6"), 
            bt[, 7], sample)
    }
    .full_join <- function(a, b) {
        as.data.frame(dplyr::full_join(a, b, by = c("left", "right")))
    }
    
    # Map Counts
    pets <- Reduce(.full_join, petlist)
    iraw <- pets[, c("left", "right")]
    iraw <- t(apply(iraw, 1, function(x) {
        if (x[1] < x[2]) {
            x
        } else {
            x[c(2, 1)]
        }
    }))
    interactions <- iraw[order(iraw[, 1], iraw[, 2]), ]
    colnames(interactions) <- c("left", "right")
    counts <- as.matrix(pets[, -c(1:2)])[order(iraw[, 1], iraw[, 2]), ]
    counts[is.na(counts)] <- 0
    counts <- as.matrix(counts, ncol = length(samples))
    colnames(counts) <- samples
    
    # Initialize rowData slot (with loop widths)
    w <- (start(anchors[interactions[, 2]]) + end(anchors[interactions[, 2]]))/2 -
         (start(anchors[interactions[, 1]]) + end(anchors[interactions[, 1]]))/2
    w[w < 0] <- 0
    rowData <- as.data.frame(as.integer(w))
    colnames(rowData) <- c("loopWidth")
    
    # Remove 'chr' from anchors
    seqlevels(anchors) <- gsub("^chr(.*)$", "\\1", seqlevels(anchors))
    
    # Remove rownames from matrices
    row.names(interactions) <- NULL
    row.names(counts) <- NULL
    
    # Initialize colData slot
    groups <- rep("group1", length(samples))
    if(length(samples) == 1){
        sizeFactor <- 1
    } else {
        lc <- log2(counts)
        keep <- rowSums(counts > 0) == ncol(lc)
        lc <- lc[keep, ]
        target <- 2^rowMeans(lc)
        sizeFactor <- colMedians(sweep(2^lc, 1, target, FUN = "/"), na.rm = TRUE)
    }
    dfcd <- data.frame(sizeFactor, groups)
    rownames(dfcd) <- samples
    
    # Create loops object
    dlo <- loops()
    slot(dlo, "anchors", check = TRUE) <- anchors
    slot(dlo, "interactions", check = TRUE) <- interactions
    slot(dlo, "counts", check = TRUE) <- counts
    slot(dlo, "colData", check = TRUE) <- dfcd
    slot(dlo, "rowData", check = TRUE) <- rowData
    
    return(dlo)
}

#' @rdname loopsMake.mango
setMethod(f = "loopsMake.mango", def = function(beddir, samples=NA, mergegap = 500, ext = "all") {
    if (sum(is.na(samples)) > 0) {
        file.path(beddir, paste("", ext, "mango", sep = "."))
        samples <- dir(beddir, pattern = paste("", ext, "mango", sep = "."))
        samples <- sub(paste(".interactions", ext, "mango", sep = "."), "", samples)
    }
    .loopsMake.mango(beddir, samples, mergegap, ext)
})

