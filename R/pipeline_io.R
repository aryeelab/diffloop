#' @include core.R
NULL

#' Read preprocessed ChiA-PET data
#'
#' \code{loopdataMake} reads in a data directory created by the
#' \code{dnaloop} preprocessing pipeline and returns a loopdata object
#'
#' This function reads in preprocessed ChIA-PET data produced by the
#' \code{dnaloop} preprocessing pipeline. The preprocessed directory 
#' contains one subdirectory per sample. The \code{samples} argument specifies
#' which samples are read. if \code{samples} is not specified all samples will
#' be read. \code{type} restricts loops whether they are on the same "inter" or
#' different "intra" chormosome. Default is "all"
#'
#' @param beddir A string. The preprocessed data directory
#' @param samples A character vector. Optional list of samples to read in
#' @param mergegap An integer value of the radius to merge anchors; default 0
#' @param type Specificies "intra", "inter", or "all" looping. Default "all"
#'
#' @return A loopdata object
#'
#' @examples
#' # Reading in all samples, no mergegap, all loops
#' bd<- system.file('extdata', 'esc_jurkat', package='diffloopdata')
#' loops <- loopdataMake(bd)
#'
#' # Reading in a subset of samples, 1kb mergegap, only intrachromosomal
#' # looping
#' samples <- c('naive_esc_1', 'naive_esc_2')
#' naive.intra <- loopdataMake(bd, samples, 1000, "intra")
#'

#' @import foreach
#' @importFrom utils read.table stack
#' @importFrom dplyr group_by full_join
#' @import readr

#' @export
setGeneric(name = "loopdataMake", def = function(beddir, samples = NA, 
    mergegap = 0, type = "all") standardGeneric("loopdataMake"))

#' @import GenomicRanges
.loopdataMake <- function(beddir, samples, mergegap, type) {
    
    ct <- list(col_character(), col_integer(), col_integer(), 
        col_character(), col_integer(), col_integer(), col_character(), 
        col_integer())
    
    restrictPets <- function(bt){
        if (type == "intra"){ return(bt[bt$X1 == bt$X4, ])
        } else if (type == "inter"){ return(bt[bt$X1 != bt$X4, ])
        } else { return(bt)}
    }
    
    # Iterate through files to set up anchors
    anchorsraw <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "loop_counts.bedpe", 
            sep = "."))
        bt <- read_delim(fullfile, " ", col_types = ct, col_names = FALSE)
        bt <- restrictPets(bt)
        plyr::rbind.fill(bt[, 1:3], setNames(bt[, 4:6], names(bt[, 
            1:3])))
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
        g <- dplyr::group_by(df, left, right)
        d <- cbind(g, counts)
        colnames(d) <- c("left", "right", sample)
        d
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
        dplyr::full_join(a, b, by = c("left", "right"))
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
    loops <- iraw[order(iraw[, 1], iraw[, 2]), ]
    colnames(loops) <- c("left", "right")
    counts <- as.matrix(pets[, -c(1:2)])[order(iraw[, 1], iraw[, 
        2]), ]
    counts[is.na(counts)] <- 0
    
    # Initialize colData slot
    sizeFactor <- rep(1, length(samples))
    groups <- rep("group1", length(samples))
    dfcd <- data.frame(sizeFactor, groups)
    rownames(dfcd) <- samples
    
    # Create loopdata object
    dlo <- loopdata()
    slot(dlo, "anchors", check = TRUE) <- anchors
    slot(dlo, "loops", check = TRUE) <- as.matrix(loops)
    slot(dlo, "counts", check = TRUE) <- counts
    slot(dlo, "colData", check = TRUE) <- dfcd
    
    return(dlo)
}

#' @rdname loopdataMake
setMethod(f="loopdataMake", def=function(beddir, samples, mergegap = 0, type="all") {
        if(sum(is.na(samples)) > 0){
            samples <- dir(beddir, pattern = ".loop_counts.bedpe")
            samples <- sub(".loop_counts.bedpe", "", samples)
        }
        .loopdataMake(beddir, samples, mergegap, type)
    })
