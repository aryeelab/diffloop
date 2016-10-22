#' @include sugar.R
NULL

#' Generalized differential Loop Calling
#'
#' \code{loopAssoc} takes a \code{loops} object and prepares it for the
#' returns another \code{loops} object with summary statistics
#' per-loop in the \code{rowData} 
#'
#' By the default, we generate is to generate a design matrix from
#' \code{loops@colData$groups}. Currently, 'edgeR' and 'Voom' are the
#' two supported association  
#' methods, but new association tests may be added in later developments.
#'
#' @param y A loops object for association
#' @param method Specifies association; either "Voom" or "edgeR"
#' @param design A design matrix of the samples; required for "Voom"
#' @param coef A vector for the coefficient of GLM. See edgeR manual
#' @param contrast A vector for the contrast. See edgeR manual
#'
#' @return A loops object
#'
#' @examples
#' # Differential loop fit
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' # assoc <- loopAssoc(loops.small, coef = 2)

#' @import edgeR
#' @import locfit
#' @import statmod
#' @import limma
#' @export
setGeneric(name = "loopAssoc", def = function(y, method = "edgeR", design = NULL, coef = NULL, contrast = NULL)
    standardGeneric("loopAssoc"))

#' @rdname loopAssoc
setMethod(f = "loopAssoc", signature = c("loops", "ANY", "ANY", "ANY", "ANY"), 
    definition = function(y, method = "edgeR", design = NULL,  coef = NULL, contrast = NULL) {
        stopifnot(method == "edgeR" | method == "Voom")
        if(is.null(coef) & is.null(contrast) & method == "edgeR"){
            stop("specify either the coefficient or the contrast when using the edgeR method")
        }
        
        groups <- y@colData$groups
        z <- DGEList(counts = y@counts, group = groups)
        z <- calcNormFactors(z)
        if(is.null(design)) design <- model.matrix(~groups)
        cat("The coefficients of the fitted GLM object are:\n")
        cat(colnames(model.matrix(~groups)))
        
        if(method == "edgeR"){
            yy <- estimateDisp(z, design)
            fit <- glmQLFit(yy, design, robust = TRUE)
            if(!is.null(coef)){
                qlf <- glmQLFTest(fit, coef = coef)
            } else {
                qlf <- glmQLFTest(y@fit, contrast = contrast)
            }
            results <- as.data.frame(topTags(qlf, n = nrow(y@counts), sort.by = "none"))
        } else {
            v <- voom(z,design,plot=FALSE)
            fit <- lmFit(v, design)
            fit <- eBayes(fit, robust = TRUE)
            results <- as.data.frame(topTable(fit, number = nrow(y@counts), sort.by = "none"))
        }
        newRowData <- as.data.frame(cbind(y@rowData, results))
        row.names(newRowData) <- NULL
        y@rowData <- newRowData
        return(y)
    })

#' Combined association test for all loops in a defined region 
#'
#' \code{slidingWindowTest} takes a \code{loops} object and 
#' integer values of the association window and the distance between 
#' consecutive windows. 
#'
#' This function returns a data.frame sorted by FDR of each region. The engine
#' loops over each chromosome and defines the first window at the left-most
#' loop and slides the window right until no more loops are present in \code{x}
#' Each region is determined from a sliding window of fixed length.
#' The combined significance measure per feature is computed via the Simes 
#' method for intrachromosomal loops where at least one anchor from the loop
#' overlaps with the region. Requires PValue column in the rowData slot.
#' 
#' @param x A loops object with PValue column (from association testing)
#' @param window The length a window will be for combined association
#' @param step The size that the window will shift for each association
#'
#' @return A data.frame sorted by FDR
#'
#' @examples
#' # Sliding window test 100kb at a time between naive and jurkat
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' # assoc_jn <- loopAssoc(loops.small, coef = 2)
#' # sw_jn <- slidingWindowTest(assoc_jn, 100000, 100000)
#' 
#' @import GenomicRanges
#' @importFrom stats complete.cases model.matrix p.adjust prcomp setNames
#' @export
setGeneric(name = "slidingWindowTest", def = function(x, window, 
    step) standardGeneric("slidingWindowTest"))

#' @rdname slidingWindowTest
setMethod(f = "slidingWindowTest", signature = c("loops", "numeric", 
    "numeric"), definition = function(x, window, step) {
    dlo <- x
    pvals <- x@rowData$PValue
    
    # Generate data.frame of anchor locations and pvalues for
    # loops
    keepcols <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", 
        "end_2")
    bigdf <- cbind(summary(dlo)[keepcols], pvals)
    df.use <- bigdf[bigdf$chr_1 == bigdf$chr_2, ]  #only intra
    
    all.chromosomes <- unique(c(levels(droplevels(df.use$chr_1)), 
        levels(droplevels(df.use$chr_2))))
    
    # Count how many regions will be present
    nRegions.perchr <- sapply(all.chromosomes, function(t) {
        if (any(df.use$chr_1 == t)) {
            endBP <- max(df.use[df.use$chr_1 == t, c("end_1", 
                "end_2")])
            start <- min(df.use[df.use$chr_1 == t, c("start_1", 
                "start_2")])
            nRegions <- round((endBP - start)/step) + 1
            nRegions
        } else {
            0
        }
    })
    
    # Define results data.frame
    resdf <- data.frame(matrix(NA, nrow = sum(nRegions.perchr), 
        ncol = 6))
    
    # i row of resdf; j element of all.chromosomes; kth step
    for (j in 1:length(all.chromosomes)) {
        chrom <- all.chromosomes[j]
        last <- sum(nRegions.perchr[1:j])
        i <- last - nRegions.perchr[j]
        k <- 0
        chrdf <- df.use[df.use$chr_1 == all.chromosomes[j], ]
        startchr <- min(chrdf[, c("start_1", "start_2")])
        while (i < last) {
            start <- window * k + startchr
            end <- start + window
            
            # Determine if either anchor overlaps with window Keep
            # pvalues where they do
            lAnchorIn <- ((chrdf$end_1 >= start & chrdf$end_1 <= 
                end) | (chrdf$start_1 >= start & chrdf$start_1 <= 
                end))
            rAnchorIn <- ((chrdf$end_2 >= start & chrdf$end_2 <= 
                end) | (chrdf$start_2 >= start & chrdf$start_2 <= 
                end))
            wPvals <- chrdf[(lAnchorIn | rAnchorIn), "pvals"]
            
            if (length(wPvals) != 0) {
                # Perform Simes Method
                r <- rank(wPvals)
                combinedP <- min(length(wPvals) * wPvals/r)
                FDR <- combinedP
                combinedResult <- cbind(combinedP, FDR)
                row <- data.frame(chrom, start, end, length(wPvals), 
                  combinedResult, stringsAsFactors = FALSE)
                resdf[i, ] <- row
            }
            i <- i + 1
            k <- k + 1
        }
    }
    # Adjust for multiple testing using FDR keep instances where
    # we actually had a hit.
    resdf <- resdf[complete.cases(resdf), ]
    resdf$X5 <- p.adjust(resdf$X5, method = "fdr")
    colnames(resdf) <- c("chr", "start", "stop", "n", "FDR", 
        "PValue")
    return(resdf[with(resdf, order(FDR)), ])
})

#' Combined association test for all loops in a defined region 
#'
#' \code{featureTest} takes a \code{loops} and 
#' genomic coordinates of regions and computes combined significance 
#' metrics for each region using the Simes procedure
#'
#' This function returns a data.frame sorted by FDR of each region. Assumes
#' the region name is specified in the GRanges object with \code{id} column. 
#' Each feature is a one row in the GRanges object. The combined significance
#' measure per feature is computed via the Simes method for intrachromosomal 
#' loops where at least one anchor from the loop overlaps with the region of 
#' interest.
#' 
#' @param x A loops object 
#' @param features A GRanges object defining regions for a combined test
#'
#' @return A data.frame sorted by FDR
#'
#' @examples
#' # Human genes chromosome 1 regional association
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' # assoc <- loopAssoc(loops.small, coef = 2)
#' # Gene based association
#' # sw_jn <- featureTest(assoc, getHumanGenes(c('1')))

#' @import GenomicRanges
#' @importFrom stats complete.cases model.matrix p.adjust prcomp setNames
#' @export
setGeneric(name = "featureTest", def = function(x, features) standardGeneric("featureTest"))

#' @rdname featureTest
setMethod(f = "featureTest", signature = c("loops", "GRanges"), 
    definition = function(x, features) {
        dlo <- x
        pvals <- x@rowData$PValue
        
        # Generate data.frame of anchor locations and pvalues for
        # loops
        keepcols <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", 
            "end_2")
        bigdf <- cbind(summary(dlo)[keepcols], pvals)
        df.use <- bigdf[bigdf$chr_1 == bigdf$chr_2, ]  #intra
        featuredf <- as.data.frame(features)
        
        # Define results data.frame
        resdf <- data.frame(matrix(NA, nrow = length(features), 
            ncol = 7))
        
        for (i in 1:length(features)) {
            # Feature coordinates
            start <- featuredf[i, ]$start
            end <- featuredf[i, ]$end
            chr <- as.character(featuredf[i, ]$seqnames)
            feat <- as.character(featuredf[i, ]$id)
            
            # Grab pvalues of anchors in individual feature locus
            lAnchorIn <- ((df.use$end_1 >= start & df.use$end_1 <= 
                end) | (df.use$start_1 >= start & df.use$start_1 <= 
                end))
            rAnchorIn <- ((df.use$end_2 >= start & df.use$end_2 <= 
                end) | (df.use$start_2 >= start & df.use$start_2 <= 
                end))
            samechromo <- as.character(df.use$chr_1) == chr
            wPvals <- df.use[(lAnchorIn | rAnchorIn) & samechromo, 
                "pvals"]
            
            # Perform Simes method; update dataframe
            if (length(wPvals) != 0) {
                r <- rank(wPvals)
                combinedP <- min(length(wPvals) * wPvals/r)
                FDR <- combinedP
                combinedResult <- cbind(combinedP, FDR)
                row <- data.frame(chr, start, end, length(wPvals), 
                  feat, combinedResult, stringsAsFactors = FALSE)
                resdf[i, ] <- row
            }
            i <- i + 1
        }
        # Adjust for multiple testing using FDR keep instances where
        # we actually had a hit.
        resdf <- resdf[complete.cases(resdf), ]
        resdf$X6 <- p.adjust(resdf$X6, method = "fdr")
        colnames(resdf) <- c("chr", "start", "stop", "n", "feature", 
            "FDR", "PValue")
        return(resdf[with(resdf, order(FDR)), ])
    })

#' Perform quick differential loop calling
#'
#' \code{quickAssoc} takes a loops object and performs a basic
#' \code{edgeR} association on the counts matrix and groups from \code{colData}
#'
#' This function returns the output of fitting an edgeR model using
#' the groups defined in \code{colData} for the specific loops
#' object. The factor normalization is based on the \code{edgeR} model.
#' For quick association, the number of groups is restricted to two. If
#' a more complex group structure exists, consider using the \code{loopAssoc}
#' function. 
#' 
#' @param y A loops object for association
#'
#' @return A loops object
#'
#' @examples
#' # Differential loop calling between naive and primed
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' np <- loops.small[,1:4]
#' assoc_np <- quickAssoc(np)

#' @import edgeR
#' @export
setGeneric(name = "quickAssoc", def = function(y) standardGeneric("quickAssoc"))

#' @rdname quickAssoc
setMethod(f = "quickAssoc", signature = c("loops"), definition = function(y) {
    # Check that there's only two groups, if not, escape
    if (length(unique(y@colData$groups)) != 2) {
        stop("Must be two groups for quickAssoc; use loopAssoc instead!")
    }
    groups <- y@colData$groups
    
    z <- DGEList(counts = y@counts, group = groups)
    design <- model.matrix(~groups)
     z <- calcNormFactors(z)
    yy <- estimateDisp(z, design)
    fit <- glmQLFit(yy, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = 2)
    results <- as.data.frame(topTags(qlf, n = nrow(y@counts), 
        sort.by = "none"))
    newRowData <- as.data.frame(cbind(y@rowData, results))
    row.names(newRowData) <- NULL
    y@rowData <- newRowData
    return(y)
})

#' Perform quick differential loop calling
#'
#' \code{quickAssocVoom} takes a loops object and performs a basic
#' \code{voom} association on the counts matrix and groups from \code{colData}
#'
#' This function returns the output of fitting an \code{voom} model using
#' the groups defined in \code{colData} for the specific loops
#' object. The factor normalization is based on the \code{voom} model.
#' For quick association, the number of groups is restricted to two. If
#' a more complex group structure exists, consider using the \code{loopAssoc}
#' function. 
#'
#' @param y A loops object for association
#'
#' @return A loops object
#'
#' @examples
#' # Differential loop calling between naive and primed
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' np <- loops.small[,1:4]
#' assoc_np_voom <- quickAssocVoom(np)

#' @import limma
#' @export
setGeneric(name = "quickAssocVoom", def = function(y) standardGeneric("quickAssocVoom"))

#' @rdname quickAssocVoom
setMethod(f = "quickAssocVoom", signature = c("loops"), definition = function(y) {
    # Check that there's only two groups, if not, escape
    if (length(unique(y@colData$groups)) != 2) {
        stop("Must be two groups for quickAssoc; use loopAssoc instead!")
    }
    groups <- y@colData$groups
    
    z <- DGEList(counts = y@counts, group = groups)
    design <- model.matrix(~groups)

    z <- calcNormFactors(z)
    v <- voom(z,design,plot=FALSE)
    fit <- lmFit(v,design)
    fit <- eBayes(fit, robust = TRUE)
    results <- as.data.frame(topTable(fit, number = nrow(y@counts), 
        sort.by = "none"))
    newRowData <- as.data.frame(cbind(y@rowData, results))
    row.names(newRowData) <- NULL
    y@rowData <- newRowData
    return(y)
})


#' Grab top loops
#'
#' \code{topLoops} takes a loops object and performs basic filtering
#' for \code{FDR} or \code{PValue}
#'
#' This function returns a subsetted \code{loops} object where all
#' loops meet the significance threshold specificed by the parameters
#' in the function call.
#'
#' @param dlo A loops object 
#' @param FDR Maximum threshold for False Discovery Rate; default = 1
#' @param PValue Maximum threshold for P-value; default = 1 
#'
#' @return A loops object subsetted by specified parameters
#'
#' @examples
#' # Differential loop calling between naive and primed
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' np <- loops.small[,1:4]
#' assoc_np <- quickAssoc(np)
#' top_np <- topLoops(assoc_np, FDR = 0.3)

#' @export
setGeneric(name = "topLoops", def = function(dlo, FDR, PValue) standardGeneric("topLoops"))

.topLoops <- function(dlo, FDR, PValue) {
    idxF <- dlo@rowData$FDR <= FDR
    idxP <- dlo@rowData$PValue <= PValue
    idxA <- idxF & idxP
    
    return(subsetLoops(dlo, idxA))
}

#' @rdname topLoops
setMethod(f = "topLoops", signature = c("loops", "numeric", "numeric"), 
    definition = function(dlo, FDR, PValue) {
        .topLoops(dlo, FDR, PValue)
    })

#' @rdname topLoops
setMethod(f = "topLoops", signature = c("loops", "numeric", "missing"), 
    definition = function(dlo, FDR, PValue) {
        .topLoops(dlo, FDR, 1)
    })

#' @rdname topLoops
setMethod(f = "topLoops", signature = c("loops", "missing", "numeric"), 
    definition = function(dlo, FDR, PValue) {
        .topLoops(dlo, 1, PValue)
    })


#' Retain loops spanning some genomic feature
#'
#' \code{filterSpanningLoops} returns a loops object where the ends of the 
#' anchors completely span one or more genomic feature (e.g. boundary)
#'
#' Rather than a simple overlap, the function by default requires a
#' genomic locus to be completely contacints
#'
#' @param dlo A loops object
#' @param gf A GRanges object of features
#'
#' @return An loops object
#'
#' @examples
#' # Return the width for loops 
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' w <- loopWidth(loops.small)
#'
#' @export
setGeneric(name = "filterSpanningLoops", def = function(dlo, gf) standardGeneric("filterSpanningLoops"))

#' @rdname filterSpanningLoops
setMethod(f = "filterSpanningLoops", signature = c("loops", "GRanges"), definition = function(dlo, gf) {
    sdf <- summary(dlo)
    span <- makeGRangesFromDataFrame(sdf, seqnames.field = "chr_1", start.field = "start_1", end.field = "end_2")
    co <- findOverlaps(gf, span, type = "within")
    return(subsetLoops(dlo, subjectHits(co)))
})

