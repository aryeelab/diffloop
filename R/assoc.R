#' @include sugar.R
NULL

#' Fit model for association testing
#'
#' \code{loopFit} takes a \code{loopdata} object and prepares it for the
#' \code{loopTest} function. 
#'
#' This function returns a \code{loopfit} object, which combines
#' the \code{loopdata} object in the input with a \code{DGEGLM} object
#' that is the normal output of an \code{edgeR glmQLFit}. To set up a 
#' different design matrix, pass that parameter through the function. 
#' Otherwise, the default is to generate a new matrix from
#' \code{loopdata@colData$groups}. Currently, 'QLF' is the only supported
#' method, but new features may be added in later developments
#'
#' @param y A loopdata object for association
#' @param design A design matrix (optional)
#' @param method Specifies association; currently only 'QLF' is supported
#'
#' @return A loopfit object
#'
#' @examples
#' # Differential loop fit
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_loopfit <- loopFit(jpn_chr1reg)
#' # Differential loop calling between naive and jurkat
#' assoc_jn <- loopTest(jpn_loopfit, coef = 2)

#' @import edgeR
#' @import locfit
#' @import statmod
#' @export
setGeneric(name = "loopFit", def = function(y, design, method = "QLF") standardGeneric("loopFit"))

#' @rdname loopFit
setMethod(f = "loopFit", signature = c("loopdata", "missing", 
    "missing"), definition = function(y, design, method) {
    groups <- y@colData$groups
    z <- DGEList(counts = y@counts, group = groups)
    z <- calcNormFactors(z)
    design <- model.matrix(~groups)
    cat("The coefficients of the fitted GLM object are:\n")
    cat(colnames(model.matrix(~groups)))
    yy <- estimateDisp(z, design)
    fit <- glmQLFit(yy, design, robust = TRUE)
    return(loopfit(loopdata = y, fit = fit))
})

#' @rdname loopFit
setMethod(f = "loopFit", signature = c("loopdata", "matrix", 
    "missing"), definition = function(y, design, method) {
    groups <- y@colData$groups
    z <- DGEList(counts = y@counts, group = groups)
    z <- calcNormFactors(z)
    yy <- estimateDisp(z, design)
    fit <- glmQLFit(yy, design, robust = TRUE)
    return(loopfit(loopdata = y, fit = fit))
})

#' Compute looptest object
#'
#' \code{loopTest} takes a \code{loopfit} object from the 
#' \code{loopFit} function and creates a \code{looptest} object.
#'
#' This function returns a \code{looptest} object, which combines the results
#' from an association with the \code{loopdata}. The default association is
#' using coefficient 2 from the model matrix (e.g. good for pair
#' comparisons) but the user may specify a different coefficient. Currently, 
#' 'QLF' is the only supported method, but new features may be added in
#' later developments. Users may also specify the contrast between the columns
#' in the design matrix as used in \code{edgeR}.
#'
#' @param y A loopfit object for association
#' @param coef Specifies coefficient of design matrix
#' @param contrast Specifies comparison of groups from design matrix
#' @param method Specifies association method; only QLF is currently supported
#'
#' @return A looptest object
#'
#' @examples
#' # Differential loop fit
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_loopfit <- loopFit(jpn_chr1reg)
#' # Differential loop calling between naive and jurkat
#' assoc_jn <- loopTest(jpn_loopfit, coef = 2)

#' @import edgeR
#' @export 
setGeneric(name = "loopTest", function(y, coef = 2, contrast, 
    method = "QLF") standardGeneric("loopTest"))

#' @rdname loopTest
setMethod(f = "loopTest", signature = c("loopfit", "missing", 
    "missing", "missing"), definition = function(y, coef, contrast, 
    method) {
    qlf <- glmQLFTest(y@fit, coef = 2)
    return(looptest(loopdata = y@loopdata, results = as.data.frame(topTags(qlf, 
        n = nrow(y@loopdata@counts), sort.by = "none"))))
})

#' @rdname loopTest
setMethod(f = "loopTest", signature = c("loopfit", "numeric", 
    "missing", "missing"), definition = function(y, coef, contrast, 
    method) {
    qlf <- glmQLFTest(y@fit, coef)
    return(looptest(loopdata = y@loopdata, results = as.data.frame(topTags(qlf, 
        n = nrow(y@loopdata@counts), sort.by = "none"))))
})

#' @rdname loopTest
setMethod(f = "loopTest", signature = c("loopfit", "missing", 
    "numeric", "missing"), definition = function(y, coef, contrast, 
    method) {
    qlf <- glmQLFTest(y@fit, contrast = contrast)
    return(looptest(loopdata = y@loopdata, results = as.data.frame(topTags(qlf, 
        n = nrow(y@loopdata@counts), sort.by = "none"))))
})

#' Combined association test for all loops in a defined region 
#'
#' \code{slidingWindowTest} takes a \code{looptest} object and 
#' integer values of the association window and the distance between 
#' consecutive windows. 
#'
#' This function returns a data.frame sorted by FDR of each region. The engine
#' loops over each chromosome and defines the first window at the left-most
#' loop and slides the window right until no more loops are present in \code{x}
#' Each region is determined from a sliding window of fixed length.
#' The combined significance measure per feature is computed via the Simes 
#' method for intrachromosomal loops where at least one anchor from the loop
#' overlaps with the region.
#' 
#' @param x A looptest object 
#' @param window The length a window will be for combined association
#' @param step The size that the window will shift for each association
#'
#' @return A data.frame sorted by FDR
#'
#' @examples
#' # Sliding window test 100kb at a time between naive and jurkat
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_loopfit <- loopFit(jpn_chr1reg)
#' # Differential loop calling between naive and jurkat
#' assoc_jn <- loopTest(jpn_loopfit, coef = 2)
#' sw_jn <- slidingWindowTest(assoc_jn, 100000, 100000)
#' 
#' @import GenomicRanges
#' @importFrom stats complete.cases model.matrix p.adjust prcomp setNames
#' @export
setGeneric(name = "slidingWindowTest", def = function(x, window, 
    step) standardGeneric("slidingWindowTest"))

#' @rdname slidingWindowTest
setMethod(f = "slidingWindowTest", signature = c("looptest", 
    "numeric", "numeric"), definition = function(x, window, step) {
    dlo <- x@loopdata
    pvals <- x@results$PValue
    
    # Generate data.frame of anchor locations and pvalues for
    # loops
    keepcols <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", 
        "end_2")
    bigdf <- cbind(summarize(dlo)[keepcols], pvals)
    df.use <- bigdf[bigdf$chr_1 == bigdf$chr_2, ]  #only intra
    
    all.chromosomes <- unique(c(levels(droplevels(df.use$chr_1)), levels(droplevels(df.use$chr_2))))
    
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
                row <- data.frame(chrom, start, end, length(wPvals), combinedResult, 
                  stringsAsFactors = FALSE)
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
    colnames(resdf) <- c("chr", "start", "stop", "n", "FDR", "PValue")
    return(resdf[with(resdf, order(FDR)), ])
})

#' Combined association test for all loops in a defined region 
#'
#' \code{featureTest} takes a \code{looptest} and 
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
#' @param x A looptest object 
#' @param features A GRanges object defining regions for a combined test
#'
#' @return A data.frame sorted by FDR
#'
#' @examples
#' # Human genes chromosome 1 regional association
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' jpn_loopfit <- loopFit(jpn_chr1reg)
#' # Differential loop calling between naive and jurkat
#' assoc_jn <- loopTest(jpn_loopfit, coef = 2)
#' # Gene based association
#' sw_jn <- featureTest(assoc_jn, getHumanGenes(c('1')))

#' @import GenomicRanges
#' @importFrom stats complete.cases model.matrix p.adjust prcomp setNames
#' @export
setGeneric(name = "featureTest", def = function(x, features) standardGeneric("featureTest"))

#' @rdname featureTest
setMethod(f = "featureTest", signature = c("looptest", "GRanges"), 
    definition = function(x, features) {
        dlo <- x@loopdata
        pvals <- x@results$PValue
        
        # Generate data.frame of anchor locations and pvalues for
        # loops
        keepcols <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", 
            "end_2")
        bigdf <- cbind(summarize(dlo)[keepcols], pvals)
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
                row <- data.frame(chr, start, end, length(wPvals), feat, combinedResult, 
                  stringsAsFactors = FALSE)
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
#' \code{quickAssoc} takes a loopdata object and performs a basic
#' \code{edgeR} association on the counts matrix and groups from \code{colData}
#'
#' This function returns the output of fitting an edgeR model using
#' the groups defined in \code{colData} for the specific loopdata
#' object. The factor normalization is based on the \code{edgeR} model.
#' For quick association, the number of groups is restricted to two. If
#' a more complex group structure exists, consider using the \code{loopFit}
#' and \code{loopTest} functions
#'
#' @param y A loopdata object for association
#'
#' @return A looptest object
#'
#' @examples
#' # Differential loop calling between naive and primed
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' np <- jpn_chr1reg[,1:4]
#' assoc_np <- quickAssoc(np)

#' @import edgeR
#' @export
setGeneric(name = "quickAssoc", def = function(y) standardGeneric("quickAssoc"))

#' @rdname quickAssoc
setMethod(f = "quickAssoc", signature = c("loopdata"), definition = function(y) {
    # Check that there's only two groups, if not, escape
    if (length(unique(y@colData$groups)) != 2) {
        stop("Must be two groups for quickAssoc; use loopFit instead!")
    }
    groups <- y@colData$groups
    z <- DGEList(counts = y@counts, group = groups)
    design <- model.matrix(~groups)
    z <- calcNormFactors(z)
    yy <- estimateDisp(z, design)
    fit <- glmQLFit(yy, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = 2)
    return(looptest(loopdata = y, results = as.data.frame(topTags(qlf, 
        n = nrow(y@counts), sort.by = "none"))))
})

#' Grab top loops
#'
#' \code{topLoops} takes a looptest object and performs basic filtering
#' for \code{FDR} or \code{PValue}
#'
#' This function returns a subsetted \code{looptest} object where all
#' loops meet the significance threshold specificed by the parameters
#' in the function call
#'
#' @param dlo A looptest object 
#' @param FDR Maximum threshold for False Discovery Rate; default = 1
#' @param PValue Maximum threshold for P-value; default = 1 
#'
#' @return A looptest object subsetted by specified parameters
#'
#' @examples
#' # Differential loop calling between naive and primed
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' np <- jpn_chr1reg[,1:4]
#' assoc_np <- quickAssoc(np)
#' top_np <- topLoops(assoc_np, FDR = 0.3)

#' @export
setGeneric(name = "topLoops", def = function(dlo, FDR, PValue) standardGeneric("topLoops"))

.topLoops <- function(dlo, FDR, PValue) {
    idxF <- dlo@results$FDR <= FDR
    idxP <- dlo@results$PValue <= PValue
    idxA <- idxF & idxP
    
    slot(dlo, "loopdata", check = TRUE) <- subsetLoops(dlo@loopdata, 
        idxA)
    slot(dlo, "results", check = TRUE) <- dlo@results[idxA, ]
    return(dlo)
}

#' @rdname topLoops
setMethod(f = "topLoops", signature = c("looptest", "numeric", 
    "numeric"), definition = function(dlo, FDR, PValue) {
    .topLoops(dlo, FDR, PValue)
})

#' @rdname topLoops
setMethod(f = "topLoops", signature = c("looptest", "numeric", 
    "missing"), definition = function(dlo, FDR, PValue) {
    .topLoops(dlo, FDR, 1)
})

#' @rdname topLoops
setMethod(f = "topLoops", signature = c("looptest", "missing", 
    "numeric"), definition = function(dlo, FDR, PValue) {
    .topLoops(dlo, 1, PValue)
})
