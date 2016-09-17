#' @include core.R
NULL

#' Visualize looping
#'
#' \code{loopPlot} takes a \code{loops} object 
#' and a \code{GRanges} objectand shows all loops in region
#' (where both anchors are present)
#'
#' Basic plot function shows the looping in each sample. The 
#' intensity of the color is proportional to the number of counts
#' observed for the particular loop relative to the other loops
#' in the entire plot. If colorLoops is specified at TRUE, then 
#' the x object must be loops and it must have a loop.type
#' column which can be generated from the \code{annotateLoops}
#' function. Blue loops are CTCF loops; black are none; red are 
#' enhancer-promoter loops; orange are promoter-promoter loops;
#' and purple are enhancer-enhancer loops. 
#' Plots use hg19 and mm9 annotation by default.
#'
#' @param x A loops object 
#' @param y A GRanges object containing region of interest
#' @param organism 'h' for human or 'm' for mouse supported
#' @param geneinfo A data.frame manually specifying annotation (see Examples)
#' @param colorLoops Differentiates loops based on loop.type in loops object
#' @param cache logic variable (default = TRUE) to use gene annotation from July.2015 freeze
#' @param maxCounts Number of counts associated with thickest loop. Default is largest count
#' in region displayed
#' @param showAnchorWidths Display the width of the anchor on the plot? Default = FALSE
#'
#' @return A plot object
#'
#' @examples
#' # Print loops in region chr1:36000000-36300000
#' library(GenomicRanges)
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' regA <- GRanges(c('1'),IRanges(start=c(36000000),end=c(36300000)))
#' plot1 <- loopPlot(loops.small, regA)
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
setGeneric(name = "loopPlot", def = function(x, y, organism = "h", 
    geneinfo = "NA", colorLoops = TRUE, cache = TRUE, showAnchorWidths = FALSE,  maxCounts = -1) standardGeneric("loopPlot"))

#' @rdname loopPlot
setMethod("loopPlot", signature(x = "loops", y = "GRanges", organism = "ANY", 
    geneinfo = "ANY", colorLoops = "ANY", cache = "ANY", maxCounts = "ANY"),
    definition = function(x, y, organism = "h", geneinfo = "NA", 
                          colorLoops = TRUE, cache = TRUE, showAnchorWidths = FALSE, maxCounts = -1) {
    return(.loopPlot(x, y, organism, geneinfo, colorLoops, cache, showAnchorWidths, maxCounts))

})

.loopPlot <- function(x, y, organism = "h", geneinfo = "NA", colorLoops = TRUE,
                      cache = TRUE, showAnchorWidths = FALSE, maxCounts = -1) {
    
    # Immediately restrict the loops object to the region
    objReg <- removeSelfLoops(subsetRegion(x, y))
    
    # Grab Regional Coordinates
    chrom <- as.character(seqnames(y))
    chromchr <- paste(c("chr", as.character(chrom)), collapse = "")
    start <- as.integer(start(ranges(range(y))))
    end <- as.integer(end(ranges(range(y))))
    
    if (geneinfo == "NA" && !cache) {
        # Get gene annotation from bioMart
        if (organism == "h") {
            mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                   path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
            chrom_biomart = gsub("chr", "", chrom)
            geneinfo = getBM(attributes = c("chromosome_name", 
                "exon_chrom_start", "exon_chrom_end", "external_gene_name", 
                "strand"), 
                mart = mart)
        } else if (organism == "m" && geneinfo == "NA" && !cache) {
            mart=useMart(host='may2012.archive.ensembl.org', 
                biomart='ENSEMBL_MART_ENSEMBL', dataset = "mmusculus_gene_ensembl")
            chrom_biomart = gsub("chr", "", chrom)
            geneinfo = getBM(attributes = c("chromosome_name", 
                "start_position", "end_position", "external_gene_id", 
                "strand"), 
                mart = mart)
        }
        # make names the same
        names(geneinfo) = c("chrom", "start", "stop", "gene", "strand")
        
        # reorder and make proper bed format
        geneinfo$score = "."
        geneinfo = geneinfo[, c(1, 2, 3, 4, 6, 5)]
    } else if (cache){
        # load and subset geneinfo; condition on mouse/human
        if(organism == "h") {
            rda <- paste(system.file("rda", package = "diffloop"), 
                "geneinfo.h.rda", sep = "/")
        } else if (organism == "m") {
            rda <- paste(system.file("rda", package = "diffloop"), 
                "geneinfo.m.rda", sep = "/")
        }
        load(rda)
        geneinfo <- geneinfo[geneinfo$chrom == chrom & geneinfo$start > start - 10000 & geneinfo$stop < end + 10000,]
    } 
    
    # Dimensions of dataframe
    n <- dim(objReg@interactions)[1]  #number of interactions
    m <- dim(objReg@counts)[2]  #number of samples
    
    res <- objReg@rowData
    
    # Setup colors for plotting
    cs <- 0
    if(!is.null(res$loop.type) & colorLoops){
        cs <- res$loop.type
        cs <- gsub("e-p", "red", cs)
        cs <- gsub("p-p", "orange", cs)
        cs <- gsub("e-e", "mediumpurple1", cs)
        cs <- gsub("ctcf", "blue", cs)
        cs <- gsub("none", "black", cs)
    } else {
        cs <- rep("black", n)
    }
    
    # Setup Dataframe for Plot
    leftAnchor <- as.data.frame(objReg@anchors[objReg@interactions[, 1]])[c(1, 2, 3)]
    LA <- do.call("rbind", replicate(m, leftAnchor, simplify = FALSE))
    rightAnchor <- as.data.frame(objReg@anchors[objReg@interactions[,2]])[c(1, 2, 3)]
    RA <- do.call("rbind", replicate(m, rightAnchor, simplify = FALSE))
    colnames(LA) <- c("chr_1", "start_1", "end_1")
    colnames(RA) <- c("chr_2", "start_2", "end_2")
    name <- rep(NA, n)
    strand_1 <- rep(".", n * m)
    strand_2 <- rep(".", n * m)
    score <- matrix(objReg@counts, ncol = 1)
    sample_id <- matrix(sapply(colnames(objReg@counts), function(x) rep(x, n)), ncol = 1)
    if(length(score != 0)){
        bedPE <- data.frame(LA, RA, name, score, strand_1, strand_2, sample_id)
    } else {
         bedPE <- data.frame(score = integer(0), sample_id = character(0), stringsAsFactors = FALSE)
    }
    # Plot
    w <- loopWidth(objReg)
    h <- sqrt(w/max(suppressWarnings(max(w)),1))
    
    malw <- 1
    if(maxCounts == -1){
        malw <- max(suppressWarnings(max(bedPE$score)), 1)
    } else {
        malw <- maxCounts
    }
    
    samples <- colnames(objReg@counts)
    lwd <- 5 * (bedPE$score/malw)
    
    loplot <- recordPlot()
    par(mfrow = c(m + 1, 1), mar = c(3, 1, 1, 1), oma = c(0, 0, 3, 0))
    for (sample in samples[-m]) {
        idx <- which(bedPE$sample_id == sample)
        bedPE_sample <- bedPE[idx , ]
        lwd_sample <- lwd[idx]
        
        ## Show anchor widths
        if(showAnchorWidths){
            loos <- bedPE_sample[bedPE_sample$score !=0,]
            
            #Make new data frame
            tdf <- unique(rbind(loos[,c(1,2,3)],setNames(loos[,c(4,5,6)], c("chr_1", "start_1", "end_1"))))
            a1df <- data.frame(
                chr_1 = tdf$chr_1, 
                start_1 = tdf$start_1,
                end_1 = tdf$start_1,
                chr_2 = tdf$chr_1,
                start_2 = tdf$end_1,
                end_2 = tdf$end_1,
                name = NA,
                score = 1,
                strand_1 = ".",
                strand_2 = ".",
                sample_id = sample
            )
            bedPE_sample <- rbind(bedPE_sample, a1df)
            
            #Update vectors
            cssamp <- c(cs, rep("black", dim(a1df)[1]))
            hsamp <- c(h, rep(0.01, dim(a1df)[1]))
            lwd_sample <- c(lwd_sample, rep(4, dim(a1df)[1]))
        } else {
            cssamp <- cs
            hsamp <- h
        }
        
        plotBedpe(bedPE_sample, chrom, start, end, color = cssamp, lwd = lwd_sample, plottype = "loops", heights = hsamp, 
                  lwdrange = c(0, 5), main = sample, adj=0)
        labelgenome(chromchr, start, end, side = 1, scipen = 20, n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    }
    sample = samples[m]
    idx <- which(bedPE$sample_id == sample)
    bedPE_sample <- bedPE[idx , ]
    lwd_sample <- lwd[idx]
    
    ## Show anchor widths
    if(showAnchorWidths & dim(bedPE_sample)[1] != 0){
        loos <- bedPE_sample[bedPE_sample$score !=0,]
        
        #Make new data frame
        tdf <- unique(rbind(loos[,c(1,2,3)],setNames(loos[,c(4,5,6)], c("chr_1", "start_1", "end_1"))))
        a1df <- data.frame(
            chr_1 = tdf$chr_1, 
            start_1 = pmax(start, tdf$start_1),
            end_1 = pmax(start, tdf$start_1), # to get the loop to fully fit on the window
            chr_2 = tdf$chr_1,
            start_2 = pmin(end, tdf$end_1),
            end_2 = pmin(end, tdf$end_1),
            name = NA,
            score = 1,
            strand_1 = ".",
            strand_2 = ".",
            sample_id = sample
        )
        bedPE_sample <- rbind(bedPE_sample, a1df)
        
        #Update vectors
        cssamp <- c(cs, rep("black", dim(a1df)[1]))
        hsamp <- c(h, rep(0.01, dim(a1df)[1]))
        lwd_sample <- c(lwd_sample, rep(4, dim(a1df)[1]))
    } else {
        cssamp <- cs
        hsamp <- h
    }
    
    plotBedpe(bedPE_sample, chrom, start, end, color = cssamp, lwd = lwd_sample, plottype = "loops", heights = hsamp, 
        lwdrange = c(0, 5), main = sample, adj=0)
    labelgenome(chromchr, start, end, side = 1, scipen = 20,  n = 3, scale = "Mb", line = 0.18, chromline = 0.5, scaleline = 0.5)
    
    if(dim(geneinfo)[1] == 0){ #Dummy plot
        plotBedpe(data.frame(), chrom, start, end, color = c("blue"), lwd = 0, 
                  plottype = "loops", heights = 0, lwdrange = c(0, 0), 
                  main = "", adj=0)
    } else {
        pg <- plotGenes(geneinfo = geneinfo, chrom = chromchr, chromstart = start, 
            chromend = end, bheight = 0.1, plotgenetype = "box", 
            bentline = FALSE, labeloffset = 0.4, fontsize = 1, arrowlength = 0.025, 
            labeltext = TRUE)
    }
    
    mtext(paste0("Region: ", chrom, ":", start, "-", end), outer = TRUE, line = 1)
    return(loplot)
}


#' Visualize sample relationships
#'
#' \code{pcaPlot} takes a loops object plots the individual samples
#' based on the principal components of the loop counts matrix
#'
#' Groups for the principal component plots are derived from \code{colData}
#' and the normalizing factors are also taken from \code{colData}. While some
#' loops objects may have non-informative groups or size factors, they 
#' should always be present. 
#'
#' @param dlo A loops object 
#'
#' @return A ggplot2 plot
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' p1 <- pcaPlot(loops.small)
#' @import ggplot2
#' 
#' @export
setGeneric(name = "pcaPlot", def = function(dlo) standardGeneric("pcaPlot"))

#' @rdname pcaPlot
setMethod(f = "pcaPlot", signature = c("loops"), definition = function(dlo) {
    pca <- prcomp(t(log2(1 + sweep(dlo@counts, 2, dlo@colData$sizeFactor, 
        FUN = "/"))))
    groups <- dlo@colData$groups
    p1 <- qplot(pca$x[, 1], pca$x[, 2], xlab = "Principal Component 1", 
        ylab = "Principal Component 2", col = groups) + theme_bw()
})


#' Plot the most significant loops
#'
#' \code{plotTopLoops} takes a loops object and creates a time-stamped .pdf
#' file with loop plots (one per page) of the top loops. 
#'
#' Each plot will show the region +/- 1 loopwidth of the loop with annotation specified
#' for either human or mouse. Assumes columns Pvalue and FDR are specified in the loops
#' object. We recommend removing self loops before using this function (and in reality, 
#' before any association testing was called.)
#'
#' @param lto loops object 
#' @param n number of loops to print (can remain 0 to specify from other parameters)
#' determined by PValue
#' @param PValue Maximum pvalue threshold for loop inclusion when printing loop plot
#' @param FDR False discovery rate threshold for inclusion
#' @param organism Either 'm' for mouse or 'h' for human. 
#' @param colorLoops Default FALSE; specify true if rowData slot contains
#' loop.type from annotateLoops to visualize plots with varying colors for
#' CTCF looping and enhancer-promoter looping
#'
#' @return Prints a time stamped .pdf file of top loops
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' jpn.u <- removeSelfLoops(loops.small)
#' jpn_loopfit <- loopFit(jpn.u)
#' assoc_jn <- loopTest(jpn_loopfit, coef = 2)
#' plotTopLoops(assoc_jn, n=2)
#' 
#' @import plyr
#' @import GenomicRanges
#' @import grDevices  
#' @import utils 
#' 
#' @export
setGeneric(name = "plotTopLoops", def = function(lto, n = 0, 
    PValue = 1, FDR = 1, organism = "h", colorLoops = FALSE) standardGeneric("plotTopLoops"))

#' @rdname plotTopLoops
setMethod(f = "plotTopLoops", signature = c("loops", "ANY", "ANY", 
    "ANY", "ANY", "ANY"), definition = function(lto, n = 0, PValue = 1, 
    FDR = 1, organism = "h", colorLoops = FALSE) {
    if (n > 0) {
        if (n > dim(lto)[2]) {
            stop("Too many loops to print; there aren't that many in the data!")
        }
        d <- lto@rowData
        tl <- subsetLoops(lto, as.integer(rownames(head(d[order(d$PValue), , drop = FALSE], n = n))))
    } else if (FDR < 1 | PValue < 1) {
        tl <- topLoops(lto, FDR = FDR, PValue = PValue)
        n <- dim(tl)[2]
    } else {
        stop("Specify valid parameters (either n > 0 or PValue/FDR < 1")
    }
    fname <- gsub(":", ".", gsub(" ", "at", paste("top-", as.character(n), 
        "loops-", Sys.time(), ".pdf", sep = ""), fixed = TRUE))
    pdf(file = fname)
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    
    for (i in 1:n) {
        one <- subsetLoops(tl, i)
        lw <- loopWidth(one)
        regPlot <- GRanges(c(one@anchors[1]@seqnames), IRanges(c(start(one@anchors[1]@ranges)), 
            c(end(one@anchors[1]@ranges))))
        regPlot <- padGRanges(regPlot, pad = lw * 1.5)
        loopPlot(lto, regPlot, organism = organism, colorLoops = colorLoops)
        setTxtProgressBar(pb, i)
    }
    dev.off()
    close(pb)
})

#' Plot several loop regions
#'
#' \code{manyLoopPlots} takes a loops object and creates a time-stamped .pdf
#' file with loop plots (one per page) of all regions specified in the GRanges object. 
#'
#' Each plot will show one region sequentially that is supplied in the GRanges object.
#'
#' @param x loops object 
#' @param y GRanges object with many regions to be visualized
#' @param organism 'h' for human or 'm' for mouse supported
#' @param geneinfo A data.frame manually specifying annotation (see Examples)
#' @param colorLoops Differentiates loops based on loop.type in loops object
#' @param cache logic variable (default = TRUE) to use gene annotation from July.2015 freeze
#' @param maxCounts Number of counts associated with thickest loop. Default is largest count
#' in region displayed
#' @param showAnchorWidths Display the width of the anchor on the plot? Default = FALSE
#'
#' @return Prints a time stamped .pdf file of top loops
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' jpn.u <- removeSelfLoops(loops.small)
#' jpn_loopfit <- loopFit(jpn.u)
#' assoc_jn <- loopTest(jpn_loopfit, coef = 2)
#' #manyLoopPlots(assoc_jn, regs) #define regs as multiple GRanges
#' 
#' @import GenomicRanges
#' 
#' @export
setGeneric(name = "manyLoopPlots", def = function(x, y, organism = "h", 
    geneinfo = "NA", colorLoops = FALSE, cache = TRUE, maxCounts = -1, showAnchorWidths = FALSE) standardGeneric("manyLoopPlots"))

#' @rdname manyLoopPlots
setMethod("manyLoopPlots", signature(x = "loops", y = "GRanges", organism = "ANY", 
    geneinfo = "ANY", colorLoops = "ANY", cache = "ANY", maxCounts = "ANY", showAnchorWidths = "ANY"),
    definition = function(x, y, organism = "h", geneinfo = "NA", colorLoops = FALSE, cache = TRUE,
                          maxCounts = -1, showAnchorWidths = FALSE) {
   
    fname <- gsub(":", ".", gsub(" ", "at", paste0("manyLoopsPlotted-", Sys.time(), ".pdf"), fixed = TRUE))
    pdf(file = fname)
    pb <- txtProgressBar(min = 0, max = length(y), style = 3)
    
    for (i in 1:length(y)) {
        yy <- y[i]
        loopPlot(x, yy, organism = organism, colorLoops = colorLoops, maxCounts = maxCounts, showAnchorWidths = showAnchorWidths)
        setTxtProgressBar(pb, i)
    }
    dev.off()
    close(pb)
})