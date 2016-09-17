#' @include union.R
NULL

#' Add meta data column to anchors based on GRanges
#'
#' \code{annotateAnchors} adds a logical variable to meta data columns in the
#' anchors based on a GRanges object of features' genomic coordinates
#'
#' This function adds column of TRUE/FALSE values on the loops object
#' anchors whether a feature is observed nearby in \code{features}. The name
#' of this column that will be in the anchors GRanges object is specified by
#' a user defined string \code{featureName}. Gap tolerance between a feature
#' and an anchor is specified by \code{maxgap}, where the default is 1,000bp.
#'
#' @param dlo A loops object whose anchors will be annotated
#' @param features A Granges object corresponding to locations of interest
#' @param featureName A string that will be the mcol name in anchors
#' @param maxgap A value of max permissible gap between a feature and anchor
#'
#' @return A loops object with new meta data column in anchors
#'
#' @examples
#' # Annotate whether anchors are near a gene body; within 1kb
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' gb <-getHumanGenes()
#' loops.small <- annotateAnchors(loops.small,gb,'nearGeneBody')
#'
#' # Adding close to gene bodies with no gap tolerance
#' loops.small <- annotateAnchors(loops.small,gb,'inGeneBody',0)
#'
#' @import GenomicRanges
#' 
#' @export
setGeneric(name = "annotateAnchors", def = function(dlo, features, 
    featureName, maxgap) standardGeneric("annotateAnchors"))

.annotateAnchors <- function(dlo, features, featureName, maxgap) {
    hits <- suppressWarnings(findOverlaps(features, dlo@anchors, 
        maxgap = maxgap))
    idx <- unique(subjectHits(hits))
    values <- data.frame(matrix(FALSE, ncol = 1, nrow = length(ranges(dlo@anchors))))
    values[idx, ] <- TRUE
    colnames(values) <- featureName
    mcols(dlo@anchors) <- c(mcols(dlo@anchors), values)
    return(dlo)
}

#' @rdname annotateAnchors
setMethod(f = "annotateAnchors", signature = c("loops", "GRanges", 
    "character", "missing"), definition = function(dlo, features, 
    featureName, maxgap = 1000) {
    maxgap <- 1000
    .annotateAnchors(dlo, features, featureName, maxgap)
})

#' @rdname annotateAnchors
setMethod(f = "annotateAnchors", signature = c("loops", "GRanges", 
    "character", "numeric"), definition = function(dlo, features, 
    featureName, maxgap) {
    .annotateAnchors(dlo, features, featureName, maxgap)
})

#' Add meta data column to anchors based on .bigwig
#'
#' \code{annotateAnchors.bigwig} adds a numeric variable to meta data 
#' columns in the anchors slot based on a user-specified .bigwig
#' file. 
#'
#' This function adds a meta data column to anchors of the specified
#' loops object. All values from the .bigwig file that overlap with the 
#' each anchor are handled by the FUN (default is to average them) to 
#' produce a single value added to the mcols of the anchors. 
#'
#' @param dlo A loops object whose anchors will be annotated
#' @param file A Granges object corresponding to locations of interest
#' @param FUN A function used to combine multiple values observed in a single anchor; default is mean
#' @param pad An integer value of to pad the anchors of the loops object; default is 0
#'
#' @return A loops object with new numeric meta data column in anchors
#'
#' @examples
#' # Annotate whether anchors are near a gene body; within 1kb
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' gb <-getHumanGenes()
#' loops.small <- annotateAnchors(loops.small,gb,'nearGeneBody')
#'
#' @import GenomicRanges
#' @importFrom rtracklayer import.bw
#' @importFrom tools file_path_sans_ext
#' 
#' @export
setGeneric(name = "annotateAnchors.bigwig", function(dlo, file, FUN = mean, pad = 0)
    standardGeneric("annotateAnchors.bigwig"))

#' @rdname annotateAnchors.bigwig
setMethod(f = "annotateAnchors.bigwig", definition = function(dlo, file, FUN = mean, pad = 0) {
    
    sample <- basename(file_path_sans_ext(file))
    bw.vals <- import.bw(file, which = addchr(dlo@anchors))
    ovl.k <- findOverlaps(addchr(dlo@anchors), bw.vals)
    qh.k <- queryHits(ovl.k) # anchors
    sh.k <- subjectHits(ovl.k) 
    values.t <- as.data.frame(tapply(mcols(bw.vals[sh.k])$score, qh.k, FUN))
    
    #A lot of extra effort to handle anchor regions with no values
    colnames(values.t) <- "bwvalues"
    vNA <- data.frame(matrix(NA, ncol = 1, nrow = length(ranges(dlo@anchors))))
    colnames(vNA) <- "NAss"
    ugly <- merge(vNA, values.t, by=0, all = TRUE, sort = F)
    ugly <- ugly[order(as.numeric(ugly$Row.names)), ]
    values <- data.frame(ugly$bwvalues)
    colnames(values) <- sample
    mcols(dlo@anchors) <- as.data.frame(c(mcols(dlo@anchors), values))
    return(dlo)
})

#' Get protein coding gene regions
#'
#' \code{getHumanGenes} returns a \code{GRanges} object of all protein
#' coding genes genome-wide or within specified chromosomes. Annotation
#' is from regions from hg19/Gr37 and protein coding genes.
#'
#' This function returns a \code{GRanges} object with the coordinates and
#' gene IDs of all protein coding genes either genome-wide 
#' (by default) orspecified within a particular chromosome. 
#'
#' @param chr A vector of chromosomes 
#' @param cache logic variable (default = TRUE) to use genes from July.2015 freeze
#'
#' @return A GRanges object
#'
#' @examples
#' # Grab all protein coding gene locations genome-wide
#' pc.genes <- getHumanGenes()
#' # Grab all protein coding gene loctions on chromosome 1
#' chr1 <- getHumanGenes(c('1'))
#' @import GenomicRanges
#' @import biomaRt
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels seqlevels<- 
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @export
setGeneric(name = "getHumanGenes", def = function(chr, cache = TRUE) standardGeneric("getHumanGenes"))

#' @import GenomicRanges
.getHumanGenesNoCache <- function(chr) {
    vals = list(chr, "protein_coding")
    mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                   path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
    geneinfo = getBM(attributes = c("chromosome_name", "start_position", 
        "end_position", "external_gene_name"), filters = c("chromosome_name", 
        "biotype"), values = vals, mart = mart)
    colnames(geneinfo) <- c("chr", "start", "end", "id")
    raw <- makeGRangesFromDataFrame(geneinfo, keep.extra.columns = TRUE, 
        ignore.strand = TRUE, seqnames.field = c("chr"), start.field = "start", 
        end.field = c("end"), starts.in.df.are.0based = FALSE)
    gr <- sortSeqlevels(raw)
    gr <- sort(gr)
    return(gr)
}

#' @rdname getHumanGenes
setMethod(f = "getHumanGenes", signature = c("missing", "ANY"), 
    definition = function(chr, cache = TRUE) {
        all <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
            "10", "11", "12", "13", "14", "15", "16", "17", "18", 
            "19", "20", "21", "22", "X", "Y")
        if (cache) {
            human.genes = NULL
            rda <- paste(system.file("rda", package = "diffloop"), 
                "human.genes.rda", sep = "/")
            load(rda)
            return(human.genes[as.vector(as.data.frame(human.genes)$seqnames) %in%
                as.vector(all)])
        } else {
            return(.getHumanGenesNoCache(all))
        }
    })

#' @rdname getHumanGenes
setMethod(f = "getHumanGenes", signature = c("character", "ANY"), 
    definition = function(chr, cache = TRUE) {
        if (cache) {
            human.genes = NULL
            rda <- paste(system.file("rda", package = "diffloop"), 
                "human.genes.rda", sep = "/")
            load(rda)
            return(human.genes[is.element(as.vector(as.data.frame(human.genes)$seqnames), 
                as.vector(chr))])
        } else {
            return(.getHumanGenesNoCache(chr))
        }
    })

#' Get Human Transcription Start Sites
#'
#' \code{getHumanTSS} returns a \code{GRanges} object of all 
#' transcription start sites for humans. Regions from hg19/Gr37 for 
#' protein coding regions. 
#'
#' This function returns a \code{GRanges} object with the coordinates and
#' gene TSS. The start and end of the IRanges slot will be the same number,
#' so consider using the \code{padGRanges} function after calling this function.
#'
#' @param chr Specifies what chromosomes are desired for the TSS#'  
#' @return A GRanges object
#'
#' @examples
#' # Grab all transition start sites genome-wide
#' human.TSS <- getHumanTSS()
#' @import GenomicRanges
#' @import biomaRt
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels seqlevels<- 
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @export
setGeneric(name = "getHumanTSS", def = function(chr) standardGeneric("getHumanTSS"))

#' @rdname getHumanTSS
setMethod(f = "getHumanTSS", signature = c("missing"), 
    definition = function(chr) {
        all <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                 "10", "11", "12", "13", "14", "15", "16", "17", "18", 
                 "19", "20", "21", "22", "X", "Y")
        geneinfo = NULL
        rda <- paste(system.file("rda", package = "diffloop"), 
                     "geneinfo.h.rda", sep = "/")
        load(rda)
        human.TSS <- GRanges(geneinfo[,c(1,2,3,4)])
        end(human.TSS[geneinfo$strand == 1]) <- start(human.TSS[geneinfo$strand == 1])
        start(human.TSS[geneinfo$strand == -1]) <- end(human.TSS[geneinfo$strand == -1])
        return(human.TSS[as.vector(as.data.frame(human.TSS)$seqnames) %in% as.vector(all)])
            
    })

#' @rdname getHumanTSS
setMethod(f = "getHumanTSS", signature = c("character"), 
    definition = function(chr) {
        geneinfo = NULL
        rda <- paste(system.file("rda", package = "diffloop"), 
                     "geneinfo.h.rda", sep = "/")
        load(rda)
        human.TSS <- GRanges(geneinfo[,c(1,2,3,4)])
        end(human.TSS[geneinfo$strand == 1]) <- start(human.TSS[geneinfo$strand == 1])
        start(human.TSS[geneinfo$strand == -1]) <- end(human.TSS[geneinfo$strand == -1])
        return(human.TSS[is.element(as.vector(as.data.frame(human.TSS)$seqnames),  as.vector(chr))])
    })


#' Annotate loops as Enhancer-Promoter or CTCF-CTCF
#'
#' \code{annotateLoops} adds a column to the rowData slot of a loops
#' object categorizing loops as either e-p (enhancer-promoter), p-p
#' (promoter-promoter), e-e (enhancer-enhancer), ctcf (CTCF-CTCF)
#' or none (no biological annotation). 
#'
#' Function annotates loops where both anchors are near CTCF peaks or where 
#' one anchor is near an enhancer and the other near a promoter. Consider using
#' functions \code{addchr}, \code{rmchr}, \code{bedToGRanges}, and  \code{padGRanges}
#' when setting up the 3 GRanges inputs. Provide a blank GRanges objects to ignore
#' classification for one set. 
#'
#' @param lto A loops object whose loops will be annotated
#' @param ctcf GRanges object corresponding to locations of CTCF peaks
#' @param enhancer GRanges object corresponding to locations of enhancer peaks
#' @param promoter GRanges object corresponding to locations of promoter regions
#'
#' @return A loops object with an additional row 'loop.type' in the rowData slot
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' ctcf_j <- system.file('extdata','Jurkat_CTCF_chr1.narrowPeak',package='diffloop')
#' ctcf <- rmchr(padGRanges(bedToGRanges(ctcf_j), pad = 1000))
#' h3k27ac_j <- system.file('extdata','Jurkat_H3K27ac_chr1.narrowPeak',package='diffloop')
#' h3k27ac <- rmchr(padGRanges(bedToGRanges(h3k27ac_j), pad = 1000))
#' promoter <- padGRanges(getHumanTSS(c('1')), pad = 1000)
#' annotated_small <- annotateLoops(loops.small, ctcf, h3k27ac, promoter)
#' 
#' @import GenomicRanges
#' @export
setGeneric(name = "annotateLoops", def = function(lto, ctcf, 
    enhancer, promoter) standardGeneric("annotateLoops"))

#' @rdname annotateLoops
setMethod(f = "annotateLoops", signature = c("loops", "GRanges", 
    "GRanges", "GRanges"), definition = function(lto, ctcf, enhancer, 
    promoter) {
    
    lto.df <- summary(lto)
    Ranchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
    Lanchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))
    
    # Determine if right anchor is near CTCF peak
    Rhits.c <- suppressWarnings(findOverlaps(ctcf, Ranchors, 
        maxgap = 0))
    Rvalues.c <- rep(FALSE, dim(lto.df)[1])
    Rvalues.c[unique(subjectHits(Rhits.c))] <- TRUE
    
    # Determine if left anchor is near CTCF peak
    Lhits.c <- suppressWarnings(findOverlaps(ctcf, Lanchors, 
        maxgap = 0))
    Lvalues.c <- rep(FALSE, dim(lto.df)[1])
    Lvalues.c[unique(subjectHits(Lhits.c))] <- TRUE
    
    ####### 
    
    # Determine if right anchor is near promoter region
    Rhits.p <- suppressWarnings(findOverlaps(promoter, Ranchors, 
        maxgap = 0))
    Rvalues.p <- rep(FALSE, dim(lto.df)[1])
    Rvalues.p[unique(subjectHits(Rhits.p))] <- TRUE
    
    # Determine if left anchor is near promoter region
    Lhits.p <- suppressWarnings(findOverlaps(promoter, Lanchors, 
        maxgap = 0))
    Lvalues.p <- rep(FALSE, dim(lto.df)[1])
    Lvalues.p[unique(subjectHits(Lhits.p))] <- TRUE
    
    ####### 
    
    # Determine if right anchor is near enhancer peak
    Rhits.e <- suppressWarnings(findOverlaps(enhancer, Ranchors, 
        maxgap = 0))
    Rvalues.e <- rep(FALSE, dim(lto.df)[1])
    Rvalues.e[unique(subjectHits(Rhits.e))] <- TRUE
    
    # Determine if left anchor is near enhancer peak
    Lhits.e <- suppressWarnings(findOverlaps(enhancer, Lanchors, 
        maxgap = 0))
    Lvalues.e <- rep(FALSE, dim(lto.df)[1])
    Lvalues.e[unique(subjectHits(Lhits.e))] <- TRUE
    
    ####### 
    
    ctcf.loops <- Lvalues.c & Rvalues.c
    ee.loops <- as.integer(Lvalues.e & Rvalues.e) * 10
    pp.loops <- as.integer(Lvalues.p & Rvalues.p) * 100
    ep.loops <- as.integer((Lvalues.e & Rvalues.p) | (Lvalues.p & Rvalues.e)) * 1000
    loop.types <- ctcf.loops + ee.loops + pp.loops + ep.loops
    description <- rep("none", length(loop.types))

    description[loop.types >= 1] <- "ctcf"
    description[loop.types >= 10] <- "e-e"
    description[loop.types >= 100] <- "p-p"
    description[loop.types >= 1000] <- "e-p"

    lto@rowData$loop.type <- description
    return(lto)
})

#' @rdname annotateLoops
setMethod(f = "annotateLoops", signature = c("loops", "missing", 
    "GRanges", "GRanges"), definition = function(lto, ctcf, enhancer, 
    promoter) {
    
    lto.df <- summary(lto)
    Ranchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
    Lanchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))

    ####### 
    
    # Determine if right anchor is near promoter region
    Rhits.p <- suppressWarnings(findOverlaps(promoter, Ranchors, 
        maxgap = 0))
    Rvalues.p <- rep(FALSE, dim(lto.df)[1])
    Rvalues.p[unique(subjectHits(Rhits.p))] <- TRUE
    
    # Determine if left anchor is near promoter region
    Lhits.p <- suppressWarnings(findOverlaps(promoter, Lanchors, 
        maxgap = 0))
    Lvalues.p <- rep(FALSE, dim(lto.df)[1])
    Lvalues.p[unique(subjectHits(Lhits.p))] <- TRUE
    
    ####### 
    
    # Determine if right anchor is near enhancer peak
    Rhits.e <- suppressWarnings(findOverlaps(enhancer, Ranchors, 
        maxgap = 0))
    Rvalues.e <- rep(FALSE, dim(lto.df)[1])
    Rvalues.e[unique(subjectHits(Rhits.e))] <- TRUE
    
    # Determine if left anchor is near enhancer peak
    Lhits.e <- suppressWarnings(findOverlaps(enhancer, Lanchors, 
        maxgap = 0))
    Lvalues.e <- rep(FALSE, dim(lto.df)[1])
    Lvalues.e[unique(subjectHits(Lhits.e))] <- TRUE
    
    ####### 
    
    ee.loops <- as.integer(Lvalues.e & Rvalues.e) * 10
    pp.loops <- as.integer(Lvalues.p & Rvalues.p) * 100
    ep.loops <- as.integer((Lvalues.e & Rvalues.p) | (Lvalues.p & Rvalues.e)) * 1000
    loop.types <- ee.loops + pp.loops + ep.loops
    description <- rep("none", length(loop.types))

    description[loop.types >= 10] <- "e-e"
    description[loop.types >= 100] <- "p-p"
    description[loop.types >= 1000] <- "e-p"

    lto@rowData$loop.type <- description
    return(lto)
})

#' Keep enhancer-promoter loops
#'
#' \code{keepEPloops} adds a column to the rowData slot of a loops
#' object that shows the corresponding TSS of a gene name based on
#' the promoter GRanges. The loops object is then subsetted and returns
#' only loops that are enhancer-promoter. 
#'
#' This function works similar to the \code{annotateLoops} function but
#' returns only enhancer-promoter loops that are defined in this function.
#' Additionally, this function returns the gene name(s) of the nearby 
#' transcription start sites in a comma-separted list if there are multiple.
#' These gene names are defined by the promoter GRanges mcol slot.
#'
#' @param lto A loops object whose loops will be annotated
#' @param enhancer GRanges object corresponding to locations of enhancer peaks
#' @param promoter GRanges object corresponding to locations of promoter regions
#'
#' @return A loops object with an additional row 'loop.type' in the rowData slot
#' in addition to the gene.tss (which has the gene name) and the
#' anchor.tss which shows the anchor(s) near the promoter region for the gene.
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' h3k27ac_j <- system.file('extdata','Jurkat_H3K27ac_chr1.narrowPeak',package='diffloop')
#' h3k27ac <- rmchr(padGRanges(bedToGRanges(h3k27ac_j), pad = 1000))
#' promoter <- padGRanges(getHumanTSS(c('1')), pad = 1000)
#' small.ep <- keepEPloops(loops.small, h3k27ac, promoter)
#' 
#' @import GenomicRanges
#' @importFrom stats aggregate
#' @export
setGeneric(name = "keepEPloops", def = function(lto, 
    enhancer, promoter) standardGeneric("keepEPloops"))

#' @rdname keepEPloops
setMethod(f = "keepEPloops", signature = c("loops", 
    "GRanges", "GRanges"), definition = function(lto, enhancer, 
    promoter) {
    
    lto.df <- summary(lto)
    Ranchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
    Lanchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))
    
    # Determine if right anchor is near promoter region
    Rhits.p <- suppressWarnings(findOverlaps(promoter, Ranchors, 
        maxgap = 0))
    Rvalues.p <- rep(FALSE, dim(lto.df)[1])
    Rvalues.p[unique(subjectHits(Rhits.p))] <- TRUE

    # Determine if left anchor is near promoter region
    Lhits.p <- suppressWarnings(findOverlaps(promoter, Lanchors, 
        maxgap = 0))
    Lvalues.p <- rep(FALSE, dim(lto.df)[1])
    Lvalues.p[unique(subjectHits(Lhits.p))] <- TRUE
    
    # Aggregate TSS
    Rtss <- data.frame(Rhits.p)
    Rtss <- cbind(Rtss, mcols(promoter[Rtss$queryHits]))
    Ltss <- data.frame(Lhits.p)
    Ltss <- cbind(Ltss, mcols(promoter[Ltss$queryHits]))
    ttss <- rbind(Rtss, Ltss)[,c(2,3)]
    tss <- aggregate(gene~subjectHits,paste,collapse=",",data=ttss)
    
    ####### 
    
    # Determine if right anchor is near enhancer peak
    Rhits.e <- suppressWarnings(findOverlaps(enhancer, Ranchors, maxgap = 0))
    Rvalues.e <- rep(FALSE, dim(lto.df)[1])
    Rvalues.e[unique(subjectHits(Rhits.e))] <- TRUE
    
    # Determine if left anchor is near enhancer peak
    Lhits.e <- suppressWarnings(findOverlaps(enhancer, Lanchors, 
        maxgap = 0))
    Lvalues.e <- rep(FALSE, dim(lto.df)[1])
    Lvalues.e[unique(subjectHits(Lhits.e))] <- TRUE
    
    ####### 
    
    #Add annotation and subset
    ep.loops <- (Lvalues.e & Rvalues.p) | (Lvalues.p & Rvalues.e)
    gene.tss <- rep("none", dim(lto)[2])
    
    anchor.tss <- rep("none", dim(lto)[2])
    anchor.tss[which(Rvalues.p & Lvalues.p)] <- "1,2"
    anchor.tss[which(Lvalues.p)] <- "1"
    anchor.tss[which(Lvalues.p)] <- "2"
    
    gene.tss[tss$subjectHits] <- tss$gene
    lto@rowData$loop.type <- "e-p"
    lto@rowData$gene.tss <- gene.tss
    lto@rowData$anchor.tss <- anchor.tss
    new.loops <- subsetLoops(lto, ep.loops)
    
    return(new.loops)
})

#' Keep CTCF loops
#'
#' \code{keepCTCFloops} returns loops that are nearby CTCF peaks
#' as determined by some external data in a GRanges object 
#'
#' This function works similar to the \code{annotateLoops} function but
#' returns only CTCF loops that are defined in this function.
#' However, loops in \code{annotateLoops} may have a different annotation
#' based on their priority scheme. For example, an e-p loop from 
#' \code{annotateLoops} may be returned as a CTCF loop by this function
#' if the loop had both annotations. These peaks don't necessarily 
#' need to be CTCF peaks, so using a GRanges object with enhancers
#' or promoters to determine e-e loops and p-p loops could also
#' be used in this function
#'
#' @param lto A loops object whose loops will be annotated
#' @param ctcf GRanges object corresponding to locations of CTCF peaks
#'
#' @return A loops object with all loops having both anchors in the GRanges region
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' ctcf_j <- system.file('extdata','Jurkat_CTCF_chr1.narrowPeak',package='diffloop')
#' ctcf <- rmchr(padGRanges(bedToGRanges(ctcf_j), pad = 1000))
#' small.ctcf <- keepCTCFloops(loops.small, ctcf)
#' 
#' @import GenomicRanges
#' @importFrom stats aggregate
#' @export
setGeneric(name = "keepCTCFloops", def = function(lto, 
    ctcf) standardGeneric("keepCTCFloops"))

#' @rdname keepCTCFloops
setMethod(f = "keepCTCFloops", signature = c("loops", 
    "GRanges"), definition = function(lto, ctcf) {
    
    lto.df <- summary(lto)
    Ranchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
    Lanchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))
    
    # Determine if right anchor is near CTCF peak
    Rhits.c <- suppressWarnings(findOverlaps(ctcf, Ranchors, 
        maxgap = 0))
    Rvalues.c <- rep(FALSE, dim(lto.df)[1])
    Rvalues.c[unique(subjectHits(Rhits.c))] <- TRUE
    
    # Determine if left anchor is near CTCF peak
    Lhits.c <- suppressWarnings(findOverlaps(ctcf, Lanchors, 
        maxgap = 0))
    Lvalues.c <- rep(FALSE, dim(lto.df)[1])
    Lvalues.c[unique(subjectHits(Lhits.c))] <- TRUE
   
    ctcf.loops <- Lvalues.c & Rvalues.c
    new.loops <- subsetLoops(lto, ctcf.loops)
    
    return(new.loops)
})


#' Annotate enhancer-promoter loops with differential gene expression
#'
#' \code{annotateLoops.dge} adds columns to the rowData slot of a loops
#' object that shows summary statistics corresponding TSS of a gene name 
#' based on the genes.tss rowData column. This function should be used
#' following the \code{keepEPloops} function. 
#'
#' This function links enhancer-promoter loops and differential gene
#' expression results. The rownames of the deseq_res slots should correspond
#' to the gene names in the gene.tss column of the rowData slot of the loops
#' object. The function returns a loops object if multiple is specified as FALSE
#' which is the case by default. Otherwise, if multiple is TRUE, then this
#' function returns a data frame since each loop may have more than moe TSS. 
#' One can reproduce this dataframe quickly when multiple = FALSE using the
#' summary() function on the returned loops object. 
#'
#' @param lto A loops object whose loops will be annotated
#' @param deseq_res A data.frame 
#' @param multiple Annotate loops with multiple TSS? Default = FALSE
#'
#' @return A loops object if multiple = FALSE or data frame if multiple = TRUE
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
#' load(rda)
#' h3k27ac_j <- system.file('extdata','Jurkat_H3K27ac_chr1.narrowPeak',package='diffloop')
#' h3k27ac <- rmchr(padGRanges(bedToGRanges(h3k27ac_j), pad = 1000))
#' promoter <- padGRanges(getHumanTSS(c('1')), pad = 1000)
#' small.ep <- keepEPloops(loops.small, h3k27ac, promoter)
#' #ADD SOMETHING HERE.
#' 
#' @import GenomicRanges
#' @importFrom data.table as.data.table
#' @export
setGeneric(name = "annotateLoops.dge", def = function(lto, 
    deseq_res, multiple = FALSE) standardGeneric("annotateLoops.dge"))

#' @rdname annotateLoops.dge
setMethod(f = "annotateLoops.dge", definition = function(lto, deseq_res, multiple = FALSE) {
    
     if(!"gene.tss" %in% colnames(lto@rowData)){
        stop("Must have gene.tss column in rowData slot!")     
    }
    
    res.dt <- as.data.table(as.data.frame(deseq_res), keep.rownames=TRUE)

    diffloop_results <- lto@rowData
    diffloop_results <- dplyr::add_rownames(diffloop_results, "idx")

    #E-P loops with one gene
    dup <- grep(",", diffloop_results$gene.tss)
    onehit <- diffloop_results[-dup,]
    
    if(!multiple){
        idx <- "" # hack to avoid note
        onehit.dt <- as.data.table(onehit)
        c.onehit <- merge(onehit.dt, res.dt, by.x = c("gene.tss"),by.y=c("rn"))
        c.onehit <- dplyr::arrange(c.onehit, as.numeric(idx))
        ida <- c.onehit$idx
        gene.tss <- c.onehit$gene.tss
        c.onehit <- within(c.onehit, rm("idx", "gene.tss"))
        c.onehit$gene.tss <- gene.tss
        new.loops <- subsetLoops(lto, as.numeric(ida))
        new.loops@rowData <- c.onehit
        return(new.loops)
        
    } else {
        diffloop_results <- summary(lto)

        #E-P loops with one gene
        dup <- grep(",", diffloop_results$gene.tss)
        onehit <- diffloop_results[-dup,]
    
        #Handle E-P loops with multiple genes
        multhits <- diffloop_results[dup,]
        multgenes <- strsplit(multhits$gene.tss, ",")
        multhits.long <- data.frame()
        for(i in 1:length(multgenes)){
            t <- cbind(data.frame(multgenes[[i]]), rep(multhits[i,], length(multgenes[[i]])))
            t$gene.tss <- t[,1]
            t <- t[,-1]
            multhits.long <- rbind(multhits.long, t)
        }

        #Combine; make data.tables
        allhits <- rbind(onehit, multhits.long)
        allhits.dt <- as.data.table(allhits)

        #Link DESeq results
        c.allhits <-  merge(allhits.dt, res.dt, by.x = c("gene.tss"),by.y=c("rn"))
        gene.tss <- c.allhits$gene.tss
        c.allhits <- within(c.allhits, rm("gene.tss"))
        c.allhits$gene.tss <- gene.tss
        return(data.frame(c.allhits))
    }
})
