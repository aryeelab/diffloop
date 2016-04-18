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
    featureName, maxgap=1000) {
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
    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
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