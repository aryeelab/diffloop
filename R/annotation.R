#' @include union.R
NULL

#' Add meta data column to anchors based on .bed file
#'
#' \code{annotateAnchors} adds a logical variable to meta data columns in the
#' anchors based on a GRanges object of features' genomic coordinates
#'
#' This function adds column of TRUE/FALSE values on the loopdata object
#' anchors whether a feature is observed nearby in \code{features}. The name
#' of this column that will be in the anchors GRanges object is specified by
#' a user defined string \code{featureName}. Gap tolerance between a feature
#' and an anchor is specified by \code{maxgap}, where the default is 1,000bp.
#'
#' @param dlo A loopdata object whose anchors will be annotated
#' @param features A Granges object corresponding to locations of interest
#' @param featureName A string that will be the mcol name in anchors
#' @param maxgap A value of max permissible gap between a feature and anchor
#'
#' @return A loopdata object with new meta data column in anchors
#'
#' @examples
#' # Annotate whether anchors are near a gene body; within 1kb
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' gb <-getHumanGenes()
#' jpn_chr1reg <- annotateAnchors(jpn_chr1reg,gb,'nearGeneBody')
#'
#' # Adding close to gene bodies with no gap tolerance
#' jpn_chr1reg <- annotateAnchors(jpn_chr1reg,gb,'inGeneBody',0)
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
setMethod(f = "annotateAnchors", signature = c("loopdata", "GRanges", 
    "character", "missing"), definition = function(dlo, features, 
    featureName, maxgap = 1000) {
    maxgap <- 1000
    .annotateAnchors(dlo, features, featureName, maxgap)
})

#' @rdname annotateAnchors
setMethod(f = "annotateAnchors", signature = c("loopdata", "GRanges", 
    "character", "numeric"), definition = function(dlo, features, 
    featureName, maxgap) {
    .annotateAnchors(dlo, features, featureName, maxgap)
})


#' Get protein coding gene regions
#'
#' \code{getHumanGenes} returns a \code{GRanges} object of all protein
#' coding genes genome-wide or within specified chromosomes
#'
#' This function returns a \code{GRanges} object with the coordinates and
#' gene IDs of all protein coding genes either genome-wide 
#' (by default) orspecified within a particular chromosome. 
#'
#' @param chr A vector of chromosomes 
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
setGeneric(name = "getHumanGenes", def = function(chr) standardGeneric("getHumanGenes"))

#' @import GenomicRanges
.getHumanGenes <- function(chr) {
    vals = list(chr, "protein_coding")
    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
        host = "jul2015.archive.ensembl.org")
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
setMethod(f = "getHumanGenes", signature = c("missing"), definition = function(chr) {
    all <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
        "11", "12", "13", "14", "15", "16", "17", "18", "19", 
        "20", "21", "22", "X", "Y")
    return(.getHumanGenes(all))
})

#' @rdname getHumanGenes
setMethod(f = "getHumanGenes", signature = c("character"), definition = function(chr) {
    return(.getHumanGenes(chr))
})

#' Get Human Transcription Start Sites
#'
#' \code{getHumanTSS} returns a \code{GRanges} object of all 
#' transcription start sites for humans
#'
#' This function returns a \code{GRanges} object with the coordinates and
#' gene TSS. The start and end of the IRanges slot will be the same number,
#' so consider using the \code{padGRanges} function after calling this function.
#'
#' @param chr Specifies what chromosomes are desired for the TSS
#'  
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
setMethod(f = "getHumanTSS", signature = c("missing"), definition = function(chr) {
    chr <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
        "11", "12", "13", "14", "15", "16", "17", "18", "19", 
        "20", "21", "22", "X", "Y")
    return(.getHumanTSS(chr))
})

#' @rdname getHumanTSS
setMethod(f = "getHumanTSS", signature = c("character"), definition = function(chr) {
    return(.getHumanTSS(chr))
})

#' @import GenomicRanges
.getHumanTSS <- function(chr) {
    vals = list(chr)
    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
        host = "jul2015.archive.ensembl.org")
    geneinfo = getBM(attributes = c("chromosome_name", "start_position", 
        "external_gene_name"), filters = c("chromosome_name"), 
        values = vals, mart = mart)
    colnames(geneinfo) <- c("chr", "start", "id")
    geneinfo$end <- geneinfo$start
    raw <- makeGRangesFromDataFrame(geneinfo, keep.extra.columns = TRUE, 
        ignore.strand = TRUE, seqnames.field = c("chr"), start.field = "start", 
        end.field = c("end"), starts.in.df.are.0based = FALSE)
    gr <- sortSeqlevels(raw)
    gr <- sort(gr)
    return(gr)
}

#' Annotate loops as Enhancer-Promoter or CTCF-CTCF
#'
#' \code{annotateLoops} adds a column to the results slot of a looptest
#' object categorizing loops as either e-p (enhancer-promoter), ctcf 
#' (CTCF-CTCF) or none (no biological annotation). If both ctcf and e-p,
#' then categorized as e-p. 
#'
#' Function annotates loops where both anchors are near CTCF peaks or where 
#' one anchor is near an enhancer and the other near a promoter. Consider using
#' functions \code{addchr}, \code{rmchr}, \code{bedToGRanges}, and  \code{padGRanges}
#' when setting up the 3 GRanges inputs. Provide a blank GRanges objects to ignore
#' classification for one set. 
#'
#' @param lto A looptest object whose loops will be annotated
#' @param ctcf GRanges object corresponding to locations of CTCF peaks
#' @param enhancer GRanges object corresponding to locations of enhancer peaks
#' @param promoter GRanges object corresponding to locations of promoter regions
#'
#' @return A looptest object with an additional row 'loop.type' in the results slot
#'
#' @examples
#' rda<-paste(system.file('rda',package='diffloop'),'jpn_chr1reg.rda',sep='/')
#' load(rda)
#' ctcf_j <- system.file('extdata','Jurkat_CTCF_chr1.narrowPeak',package = 'diffloop')
#' ctcf <- rmchr(padGRanges(bedToGRanges(ctcf_j), pad = 1000))
#' h3k27ac_j <- system.file('extdata','Jurkat_H3K27ac_chr1.narrowPeak',package = 'diffloop')
#' h3k27ac <- rmchr(padGRanges(bedToGRanges(h3k27ac_j), pad = 1000))
#' promoter <- padGRanges(getHumanTSS(c('1')), pad = 1000)
#' jn <- jpn_chr1reg[,c(1,2,5,6)]
#' assoc_jn <- quickAssoc(jn)
#' assoc_jn <- removeSelfLoops(assoc_jn)
#' annotated_jn <- annotateLoops(assoc_jn, ctcf, h3k27ac, promoter)
#' 
#' @import GenomicRanges
#' @export
setGeneric(name = "annotateLoops", def = function(lto, ctcf, 
    enhancer, promoter) standardGeneric("annotateLoops"))

#' @rdname annotateLoops
setMethod(f = "annotateLoops", signature = c("looptest", "GRanges", 
    "GRanges", "GRanges"), definition = function(lto, ctcf, enhancer, 
    promoter) {
    
    lto.df <- summary(lto)
    Ranchors <- GRanges(lto.df[, 1], IRanges(lto.df[, 2], lto.df[, 
        3]))
    Lanchors <- GRanges(lto.df[, 4], IRanges(lto.df[, 5], lto.df[, 
        6]))
    
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
    Rhits.e <- suppressWarnings(findOverlaps(ctcf, Ranchors, 
        maxgap = 0))
    Rvalues.e <- rep(FALSE, dim(lto.df)[1])
    Rvalues.e[unique(subjectHits(Rhits.e))] <- TRUE
    
    # Determine if left anchor is near enhancer peak
    Lhits.e <- suppressWarnings(findOverlaps(ctcf, Lanchors, 
        maxgap = 0))
    Lvalues.e <- rep(FALSE, dim(lto.df)[1])
    Lvalues.e[unique(subjectHits(Lhits.e))] <- TRUE
    
    ####### 
    
    ctcf.loops <- Lvalues.c & Rvalues.c
    ep.loops <- (Lvalues.e & Rvalues.p) | (Lvalues.p & Rvalues.e)
    loop.types <- as.integer(ctcf.loops) + as.integer(ep.loops) * 
        2
    loop.types <- gsub("3", "e-p", loop.types)
    loop.types <- gsub("2", "e-p", loop.types)
    loop.types <- gsub("1", "ctcf", loop.types)
    loop.types <- gsub("0", "none", loop.types)
    
    lto@results$loop.type <- loop.types
    return(lto)
})
