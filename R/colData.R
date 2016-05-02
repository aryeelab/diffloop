#' @include assoc.R
NULL

#' Compute normalizing factors for each sample
#'
#' \code{calcLDSizeFactors} takes a loopdata object computes size 
#' factors based for each sample
#'
#' This function updates the \code{loopdata} object with new
#' \code{sizeFactor} values for
#' each sample in the \code{colData} slot using a method identical to that
#' employed in \code{DESeq2.} 
#'
#' @param dlo A loopdata object with unnormalized size factors
#'
#' @return A loopdata object with new size factors in \code{colData}
#'
#' @examples
#' # Computing normalizing factors from the full ChIA-PET Data
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_chr1reg <- calcLDSizeFactors(jpn_chr1reg)

#' @importMethodsFrom matrixStats colMedians
#' @export
setGeneric(name = "calcLDSizeFactors", def = function(dlo) standardGeneric("calcLDSizeFactors"))

#' @rdname calcLDSizeFactors
setMethod("calcLDSizeFactors", c("loopdata"), function(dlo) {
    lc <- log2(dlo@counts)
    keep <- rowSums(dlo@counts > 0) == ncol(lc)
    lc <- lc[keep, ]
    target <- 2^rowMeans(lc)
    sizeFactor <- colMedians(sweep(2^lc, 1, target, FUN = "/"), 
        na.rm = TRUE)
    dlo@colData$sizeFactor <- sizeFactor
    return(dlo)
})

#' Update groups in colData for loopdata object
#'
#' \code{updateLDGroups} changes the \code{groups} column in \code{colData} for
#' a \code{loopdata} object
#'
#' This function updates the \code{groups} column in \code{colData} for
#' a \code{loopdata} object. Make sure that the length of \code{groups}
#' the number of samples in \code{colData}!
#'
#' @param dlo A loopdata object 
#' @param groups A character vector. Lists the groups each sample belongs in
#'
#' @return A loopdata object with new groups in \code{colData}
#'
#' @examples
#' # Updating groups from all 'group1' to meaningful designations
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' celltypes <- c('naive1','naive1','primed2','primed2','jurkat3','jurkat3')
#' jpn_chr1reg <- updateLDGroups(jpn_chr1reg, celltypes)
#' 
#' @export
setGeneric(name = "updateLDGroups", def = function(dlo, groups) standardGeneric("updateLDGroups"))

#' @rdname updateLDGroups
setMethod("updateLDGroups", c("loopdata"), function(dlo, groups) {
    dlo@colData$groups <- groups
    return(dlo)
})

#' @importMethodsFrom Biobase sampleNames
#' @export
setMethod("sampleNames", "loopdata", function(object) {
    rownames(object@colData)
})

#' @importMethodsFrom Biobase sampleNames<-
#' @export
setReplaceMethod("sampleNames", c("loopdata", "ANY"), function(object, value) {

    dfcd <- object@colData
    rownames(dfcd) <- value
    ncounts <- object@counts
    colnames(ncounts) <- value

    slot(object, "counts", check = TRUE) <- ncounts
    slot(object, "colData", check = TRUE) <- dfcd
    
    return(object)
})
