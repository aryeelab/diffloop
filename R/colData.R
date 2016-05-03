#' @include assoc.R
NULL

#' Compute normalizing factors for each sample
#'
#' \code{calcLDSizeFactors} takes a loops object computes size 
#' factors based for each sample
#'
#' This function updates the \code{loops} object with new
#' \code{sizeFactor} values for
#' each sample in the \code{colData} slot using a method identical to that
#' employed in \code{DESeq2.} 
#'
#' @param dlo A loops object with unnormalized size factors
#'
#' @return A loops object with new size factors in \code{colData}
#'
#' @examples
#' # Computing normalizing factors from the full ChIA-PET Data
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' loops.small <- calcLDSizeFactors(loops.small)

#' @importMethodsFrom matrixStats colMedians
#' @export
setGeneric(name = "calcLDSizeFactors", def = function(dlo) standardGeneric("calcLDSizeFactors"))

#' @rdname calcLDSizeFactors
setMethod("calcLDSizeFactors", c("loops"), function(dlo) {
    lc <- log2(dlo@counts)
    keep <- rowSums(dlo@counts > 0) == ncol(lc)
    lc <- lc[keep, ]
    target <- 2^rowMeans(lc)
    sizeFactor <- colMedians(sweep(2^lc, 1, target, FUN = "/"), 
        na.rm = TRUE)
    dlo@colData$sizeFactor <- sizeFactor
    return(dlo)
})

#' Update groups in colData for loops object
#'
#' \code{updateLDGroups} changes the \code{groups} column in \code{colData} for
#' a \code{loops} object
#'
#' This function updates the \code{groups} column in \code{colData} for
#' a \code{loops} object. Make sure that the length of \code{groups}
#' the number of samples in \code{colData}!
#'
#' @param dlo A loops object 
#' @param groups A character vector. Lists the groups each sample belongs in
#'
#' @return A loops object with new groups in \code{colData}
#'
#' @examples
#' # Updating groups from all 'group1' to meaningful designations
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' celltypes <- c('naive1','naive1','primed2','primed2','jurkat3','jurkat3')
#' loops.small <- updateLDGroups(loops.small, celltypes)
#' 
#' @export
setGeneric(name = "updateLDGroups", def = function(dlo, groups) standardGeneric("updateLDGroups"))

#' @rdname updateLDGroups
setMethod("updateLDGroups", c("loops"), function(dlo, groups) {
    dlo@colData$groups <- groups
    return(dlo)
})

#' Grab/Update Sample Names
#'
#' \code{sampleNames} takes a loops object returns the names of 
#' the samples in the structure. One can also update the names using
#' set replace.
#'
#' The examples show both accession and updating sample names. 
#'
#' @param object A loops object 
#' @param value New names when using set replace
#'
#' @return Vector of sample names
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' sampleNames(loops.small)
#' nnames <- c('one', 'two', 'three', 'four', 'five', 'six')
#' sampleNames(loops.small) <- nnames

#' @importMethodsFrom Biobase sampleNames
#' @export
setMethod("sampleNames", "loops", function(object) {
    rownames(object@colData)
})


#' @importMethodsFrom Biobase sampleNames<-
#' @rdname sampleNames-loops-method
#' @export
setReplaceMethod("sampleNames", c("loops", "ANY"), function(object, 
    value) {
    
    dfcd <- object@colData
    rownames(dfcd) <- value
    ncounts <- object@counts
    colnames(ncounts) <- value
    
    slot(object, "counts", check = TRUE) <- ncounts
    slot(object, "colData", check = TRUE) <- dfcd
    
    return(object)
})


#' Split samples into their own loops object
#'
#' \code{splitSamples} takes a loops object and returns a list of loops
#' objects where each sample populates its own loops object
#'
#' This function splits the colData and counts slots for each sample
#' but makes copies of the anchors, interactions, and rowdata
#'
#' @param dlo A loops object 
#'
#' @return A list of loops objects w
#'
#' @examples
#' # Updating groups from all 'group1' to meaningful designations
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' split <- splitSamples(loops.small)
#' 
#' @export
setGeneric(name = "splitSamples", def = function(dlo) standardGeneric("splitSamples"))

#' @rdname splitSamples
setMethod("splitSamples", c("loops"), function(dlo) {
    llo <- lapply(1:as.numeric(dim(dlo)[3]), function(i) { assign(sampleNames(dlo[,i]) , dlo[,i]) })
    names(llo) <- sampleNames(dlo)
    return(llo)
})