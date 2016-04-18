#' @include core.R
NULL

#' Visualize looping
#'
#' \code{loopPlot} takes a \code{loopdata/looptest} object 
#' and a \code{GRanges} objectand shows all loops in region
#' (where both anchors are present)
#'
#' Basic plot function shows the looping in each sample. The 
#' intensity of the color is proportional to the number of counts
#' observed for the particular loop relative to the other loops
#' in the entire plot.
#'
#' @param x A loopdata or looptest object 
#' @param y A GRanges object containing region of interest
#' @param organism 'h' for human or 'm' for mouse supported
#' @param geneinfo A data.frame manually specifying annotation (see Examples)
#'
#' @return A plot object
#'
#' @examples
#' # Print loops in region chr1:36000000-36300000
#' library(GenomicRanges)
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' regA <- GRanges(c('1'),IRanges(start=c(36000000),end=c(36300000)))
#' plot1 <- subsetRegion(jpn_chr1reg, regA)
#' #Example of \code{geneinfo} table
#' geneinfo <- data.frame(1,359345,359681,'RP5-8572K21.15','.',-1)
#' names(geneinfo) <- c('chrom','start','stop','gene','strand')
#' @import GenomicRanges
#' @import Sushi
#' @import biomaRt
#' @importFrom grDevices recordPlot
#' @importFrom graphics mtext par
#' 
#' @export
setGeneric(name = "loopPlot", def = function(x, y, organism, 
    geneinfo) standardGeneric("loopPlot"))

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "loopdata", y = "GRanges", 
    organism = "missing", geneinfo = "missing"), definition = function(x, 
    y, organism, geneinfo) {
    return(.loopPlot(x, y, "h", "NA"))
})

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "loopdata", y = "GRanges", 
    organism = "character", geneinfo = "missing"), definition = function(x, 
    y, organism, geneinfo) {
    return(.loopPlot(x, y, organism, "NA"))
})

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "loopdata", y = "GRanges", 
    organism = "missing", geneinfo = "data.frame"), definition = function(x, 
    y, organism, geneinfo) {
    return(.loopPlot(x, y, "", geneinfo))
})

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "looptest", y = "GRanges", 
    organism = "missing", geneinfo = "missing"), definition = function(x, 
    y, organism, geneinfo) {
    return(.loopPlot(x@loopdata, y, "h", "NA"))
})

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "looptest", y = "GRanges", 
    organism = "character", geneinfo = "missing"), definition = function(x, 
    y, organism, geneinfo) {
    return(.loopPlot(x@loopdata, y, organism, "NA"))
})

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "looptest", y = "GRanges", 
    organism = "missing", geneinfo = "data.frame"), definition = function(x, 
    y, organism, geneinfo) {
    return(.loopPlot(x@loopdata, y, "", geneinfo))
})

.loopPlot <- function(x, y, organism = "h", geneinfo = "NA") {
    
    # Immediately restrict the loopdata object to the region of
    # interest
    objReg <- removeSelfLoops(subsetRegion(x, y))
    
    # Grab Regional Coordinates
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    if (geneinfo == "NA") {
        # Get gene annotation from bioMart
        if (organism == "h") {
            mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
            chrom_biomart = gsub("chr", "", chrom)
            geneinfo = getBM(attributes = c("chromosome_name", 
                "exon_chrom_start", "exon_chrom_end", "external_gene_name", 
                "strand"), filters = c("chromosome_name", "start", 
                "end"), values = list(chrom_biomart, start, end), 
                mart = mart)
        } else if (organism == "m") {
            mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
            chrom_biomart = gsub("chr", "", chrom)
            geneinfo = getBM(attributes = c("chromosome_name", 
                "start_position", "end_position", "external_gene_name", 
                "strand"), filters = c("chromosome_name", "start", 
                "end"), values = list(chrom_biomart, start, end), 
                mart = mart)
        }
        # make names the same
        names(geneinfo) = c("chrom", "start", "stop", "gene", 
            "strand")
        
        # reorder and make proper bed format
        geneinfo$score = "."
        geneinfo = geneinfo[, c(1, 2, 3, 4, 6, 5)]
    }
    
    # Dimensions of dataframe
    n <- dim(objReg@loops)[1]  #number of loops
    m <- dim(objReg@counts)[2]  #number of samples
    
    # Setup Dataframe for Plot
    leftAnchor <- as.data.frame(objReg@anchors[objReg@loops[, 
        1]])[c(1, 2, 3)]
    LA <- do.call("rbind", replicate(m, leftAnchor, simplify = FALSE))
    rightAnchor <- as.data.frame(objReg@anchors[objReg@loops[, 
        2]])[c(1, 2, 3)]
    RA <- do.call("rbind", replicate(m, rightAnchor, simplify = FALSE))
    colnames(LA) <- c("chr_1", "start_1", "end_1")
    colnames(RA) <- c("chr_2", "start_2", "end_2")
    name <- rep(NA, n)
    strand_1 <- rep(".", n * m)
    strand_2 <- rep(".", n * m)
    score <- matrix(objReg@counts, ncol = 1)
    sample_id <- matrix(sapply(colnames(objReg@counts), function(x) rep(x, 
        n)), ncol = 1)
    bedPE <- data.frame(LA, RA, name, score, strand_1, strand_2, 
        sample_id)
    
    # Plot
    w <- loopWidth(objReg)
    h <- sqrt(w/max(w))
    
    samples <- colnames(objReg@counts)
    lwd <- 5 * (bedPE$score/max(bedPE$score))
    
    loplot <- recordPlot()
    par(mfrow = c(m + 1, 1), mar = c(3, 1, 1, 1), oma = c(0, 
        0, 3, 0))
    for (sample in samples[-m]) {
        idx <- which(bedPE$sample_id == sample)
        plotBedpe(bedPE[idx, ], chrom, start, end, color = rep("red", 
            n), lwd = lwd[idx], plottype = "loops", heights = h, 
            lwdrange = c(0, 5), main = sample)
        labelgenome(chromchr, start, end, side = 1, scipen = 20, 
            n = 3, scale = "Mb", line = 0.18, chromline = 0.5, 
            scaleline = 0.5)
    }
    sample = samples[m]
    idx <- which(bedPE$sample_id == sample)
    plotBedpe(bedPE[idx, ], chrom, start, end, color = rep("red", 
        n), lwd = lwd[idx], plottype = "loops", heights = h, 
        lwdrange = c(0, 5), main = sample)
    labelgenome(chromchr, start, end, side = 1, scipen = 20, 
        n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    
    pg = plotGenes(geneinfo = geneinfo, chrom = chromchr, chromstart = start, 
        chromend = end, bheight = 0.1, plotgenetype = "box", 
        bentline = FALSE, labeloffset = 0.4, fontsize = 1, arrowlength = 0.025,
        labeltext = TRUE)
    mtext(paste0("Region: ", chrom, ":", start, "-", end), outer = TRUE, 
        line = 1)
    return(loplot)
}


#' Visualize sample relationships
#'
#' \code{pcaPlot} takes a loopdata object plots the individual samples
#' based on the principal components of the loop counts matrix
#'
#' Groups for the principal component plots are derived from \code{colData}
#' and the normalizing factors are also taken from \code{colData}. While some
#' loopdata objects may have non-informative groups or size factors, they 
#' should always be present. 
#'
#' @param dlo A loopdata/looptest object 
#'
#' @return A ggplot2 plot
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' p1 <- pcaPlot(jpn_chr1reg)
#' @import ggplot2
#' 
#' @export
setGeneric(name = "pcaPlot", def = function(dlo) standardGeneric("pcaPlot"))

#' @rdname pcaPlot
setMethod(f = "pcaPlot", signature = c("loopdata"), definition = function(dlo) {
    pca <- prcomp(t(log2(1 + sweep(dlo@counts, 2, dlo@colData$sizeFactor, 
        FUN = "/"))))
    groups <- dlo@colData$groups
    p1 <- qplot(pca$x[, 1], pca$x[, 2], xlab = "Principal Component 1", 
        ylab = "Principal Component 2", col = groups) + theme_bw()
})

#' @rdname pcaPlot
setMethod(f = "pcaPlot", signature = c("looptest"), definition = function(dlo) {
    pcaPlot(dlo@loopdata)
})