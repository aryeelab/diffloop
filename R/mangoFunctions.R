#' @include union.R
NULL

# First three functions taken from mango; occassionally gently modified;
# all three are hidden and used in the mangoCorrection function
.binmaker <- function(vectortobin, binmethod = "equalocc", numberbins = 30) {
    sortedvec = sort(vectortobin)
    
    # Set up bins based on method
    if (binmethod == "equalocc") {
        itemsperbin = length(vectortobin)/numberbins
        splitdata = split(sortedvec, ceiling(seq_along(sortedvec)/itemsperbin))
        mins = c()
        maxes = c()
        for (i in (1:length(splitdata))) {
            mins = c(mins, min(splitdata[[i]]))
            maxes = c(maxes, max(splitdata[[i]]))
        }
        borders = (mins[-1] + maxes[-length(maxes)])/2
        return(borders)
    }
    
    if (binmethod == "log2") {
        log10sortedvec = log2(sortedvec)
        borders = seq(min(log10sortedvec), max(log10sortedvec), length.out = numberbins + 1)[2:numberbins]
        return(borders)
    }
    
    if (binmethod == "log10") {
        log10sortedvec = log10(sortedvec)
        borders = seq(min(log10sortedvec), max(log10sortedvec), length.out = numberbins + 1)[2:numberbins]
        return(borders)
    }
    
    if (binmethod == "equalsize") {
        borders = seq(min(sortedvec), max(sortedvec), length.out = numberbins + 1)[2:numberbins]
        return(borders)
    }
}

.model_chia <- function(x, y = NA, borders, yvals = TRUE) {
    sumofy = c()
    meanofx = c()
    sumofx = c()
    countofx = c()
    bin = findInterval(x, borders)
    for (b in 0:length(borders)) {
        binindexes = which(bin == b)
        if (yvals == FALSE) {
            sumofy = c(sumofy, length(binindexes))
        }
        if (yvals == TRUE) {
            sumofy = c(sumofy, sum(y[binindexes]))
        }
        meanofx = c(meanofx, mean(as.numeric(x[binindexes])))
        sumofx = c(sumofx, sum(as.numeric(x[binindexes])))
        countofx = c(countofx, length(binindexes))
    }
    pvals = sumofy/sum(sumofy)
    
    return(cbind(meanofx, sumofy, pvals, sumofx, countofx))
}

.makecombos <- function(chrom, chrpeaks, mindist = 0, maxdist = 1e+08) {
    chrpeaks = chrpeaks[which(chrpeaks[, 1] == chrom), ]
    
    # sort by distance
    chrpeaks = chrpeaks[order(chrpeaks[, 2]), ]
    
    npeaks = nrow(chrpeaks)
    n = npeaks^2
    
    # intialize dataframe
    start1 = numeric(n)
    end1 = numeric(n)
    start2 = numeric(n)
    end2 = numeric(n)
    score1 = numeric(n)
    score2 = numeric(n)
    
    curlength = 0
    i = 1
    
    for (i in (1:npeaks)) {
        j = i + 1
        temp_score2 = chrpeaks[j:npeaks, 5]
        temp_score1 = rep(chrpeaks[i, 5], length(temp_score2))
        temp_end2 = chrpeaks[j:npeaks, 3]
        temp_end1 = rep(chrpeaks[i, 3], length(temp_end2))
        temp_start2 = chrpeaks[j:npeaks, 2]
        temp_start1 = rep(chrpeaks[i, 2], length(temp_start2))
        
        dist = abs((temp_start1 + temp_end1)/2 - (temp_start2 + temp_end2)/2)
        keepers = which(dist < maxdist & dist > mindist)
        
        if (length(keepers) > 0) {
            first = curlength + 1
            curlength = curlength + length(keepers)
            start1[first:curlength] = temp_start1[keepers]
            end1[first:curlength] = temp_end1[keepers]
            start2[first:curlength] = temp_start2[keepers]
            end2[first:curlength] = temp_end2[keepers]
            score1[first:curlength] = temp_score1[keepers]
            score2[first:curlength] = temp_score2[keepers]
        }
    }
    
    chrpairs = data.frame(chr1 = rep(chrom, curlength),
                          start1 = start1[1:curlength],
                          end1 = end1[1:curlength],
                          chr2 = rep(chrom, curlength),
                          start2 = start2[1:curlength],
                          end2 = end2[1:curlength],
                          score1 = score1[1:curlength],
                          score2 = score2[1:curlength])
    
    return(chrpairs)
}

#' Perform mango bias correction
#'
#' \code{mangoCorrection} takes a loops object and filters loops based
#' on the binomial model used in the mango ChIA-PET pipeline.
#'
#' This function processes ChIA-PET data in a loops object and filters loops
#' that may be biased due to proximity or low PET counts as previously 
#' described by the \code{mango} pipeline. PET and anchor counts are aggregated
#' across all samples to compute statistical significance. Consider using a larger
#' number of bins (e.g. 30) for a larger data object when possible. 
#'
#' @param lo A loops object. 
#' @param FDR Minimum FDR value for loop to be included; default 1
#' @param PValue Minimum p0value for loop to be included; default 1
#' @param nbins Number of bins for mango computation
#'
#' @return A loops object where loops are filtered using mango bias correction
#'
#' @examples
#' rda <- paste(system.file("rda", package = "diffloop"), "loops.small.rda", sep = "/")
#' load(rda)
#' loops.small <- removeSelfLoops(loops.small)
#' loops.small.mango <- mangoCorrection(loops.small, PValue = 0.05)
#' 
#' @import GenomicRanges
#' @importFrom utils stack
#' @importFrom stats pbinom predict smooth.spline

#' @export
setGeneric(name = "mangoCorrection", def = function(lo, FDR = 1, PValue = 1, 
    nbins = 10) standardGeneric("mangoCorrection"))

#' @rdname mangoCorrection
setMethod(f = "mangoCorrection", def = function(lo, FDR = 1, PValue = 1, nbins = 10) {
    totalPetCounts <- rowSums(lo@counts)
    lAnchor <- cbind(as.data.frame(lo@anchors[lo@interactions[, 1]]), totalPetCounts)
    rAnchor <- cbind(as.data.frame(lo@anchors[lo@interactions[, 2]]), totalPetCounts)
    sAnchors <- subset(rbind(lAnchor, rAnchor), select = -c(width, strand))
    tc.GR <- GRanges(aggregate(totalPetCounts ~ seqnames + start + end, data = sAnchors, sum, na.rm = TRUE))
    
    # Set up chrpeaks for pairwise computation
    chrpeaks <- as.data.frame(tc.GR)
    chrpeaks <- chrpeaks[, c(1, 2, 3, 4, 6, 5)]
    names(chrpeaks) <- c("chr", "start", "end", "name", "score", "strand")
    
    # Assume ordering in tc.GR is same as the anchors; seems OK
    mcols(lo@anchors) <- mcols(tc.GR)
    
    df <- summary(lo)
    df$PETS <- totalPetCounts
    df$depths <- df$totalPetCounts_1 * df$totalPetCounts_2
    
    totalcombos <- 0
    chromosomes <- unique(as.character(seqnames(tc.GR)@values))
    
    # Distance normalization
    distanceborders <- .binmaker(df$loopWidth, binmethod = "equalocc", numberbins = nbins)
    distance_IAB_model <- .model_chia(df$loopWidth, df$PETS, borders = distanceborders, yvals = TRUE)
    distance_IAB_spline <- smooth.spline(log10(distance_IAB_model[, 1]), distance_IAB_model[, 3], spar = 0.75)
    
    # Depth normalization
    depthborders <- .binmaker(df$depths, binmethod = "equalocc", numberbins = nbins)
    depth_IAB_model <- .model_chia(x = df$depths, y = df$PETS, borders = depthborders, yvals = TRUE)
    depth_IAB_model.complete <- depth_IAB_model[complete.cases(depth_IAB_model),]
    depth_IAB_spline <- smooth.spline(log10(depth_IAB_model.complete[, 1]), depth_IAB_model.complete[, 3], spar = 0.75)
    
    # model Combos vs distance
    meanofx_dist = rep(0, nbins)
    sumofy_dist = rep(0, nbins)
    pvals_dist = rep(0, nbins)
    sumofx_dist = rep(0, nbins)
    countofx_dist = rep(0, nbins)
    meanofx_depth = rep(0, nbins)
    sumofy_depth = rep(0, nbins)
    pvals_depth = rep(0, nbins)
    sumofx_depth = rep(0, nbins)
    countofx_depth = rep(0, nbins)
    
    for (chrom in chromosomes) {
        
        # Look at all combinations of anchors
        combos = .makecombos(chrom, chrpeaks)
        combos$distance = abs((combos[, 2] + combos[, 3])/2 - (combos[, 5] + combos[, 6])/2)
        combos$depths = combos[, 7] * combos[, 8]
        distance_combo_model_chrom = .model_chia(x = combos$distance, y = NA, borders = distanceborders, yvals = FALSE)
        sumofy_dist = sumofy_dist + distance_combo_model_chrom[, 2]
        sumofx_dist = sumofx_dist + distance_combo_model_chrom[, 4]
        countofx_dist = countofx_dist + distance_combo_model_chrom[, 5]
        
        # model combinations of anchors vs depth
        depth_combo_model_chrom = .model_chia(x = combos$depths, y = NA, borders = depthborders, yvals = FALSE)
        sumofy_depth = sumofy_depth + depth_combo_model_chrom[, 2]
        sumofx_depth = sumofx_depth + depth_combo_model_chrom[, 4]
        countofx_depth = countofx_depth + depth_combo_model_chrom[, 5]
    }
    
    # combine data from all chromosomes
    depth_combo_model = cbind(sumofx_depth/countofx_depth, sumofy_depth, sumofy_depth/sum(sumofy_depth))
    depth_combo_model.complete <- depth_combo_model[complete.cases(depth_combo_model),]

    distance_combo_model = cbind(sumofx_dist/countofx_dist, sumofy_dist, sumofy_dist/sum(sumofy_dist))
    
    depth_combo_spline = smooth.spline(log10(depth_combo_model.complete[, 1]), depth_combo_model.complete[, 3], spar = 0.75)
    mOI <- is.infinite(distance_combo_model[, 1]) | is.na(distance_combo_model[, 1])
    distance_combo_spline = smooth.spline(log10(distance_combo_model[!mOI, 1]), distance_combo_model[!mOI, 3], spar = 0.75)
    
    df$P_IAB_distance = predict(distance_IAB_spline, log10(df$loopWidth))$y
    df$P_combos_distance = predict(distance_combo_spline, log10(df$loopWidth))$y
    df$P_IAB_depth = predict(depth_IAB_spline, log10(df$depths))$y
    df$P_combos_depth = predict(depth_combo_spline, log10(df$depths))$y
    
    # cap values to min and max
    df$P_IAB_distance[which(df$P_IAB_distance <= min(distance_IAB_model[, 3]))] = min(distance_IAB_model[, 3])
    df$P_IAB_distance[which(df$P_IAB_distance >= max(distance_IAB_model[, 3]))] = max(distance_IAB_model[, 3])
    df$P_combos_distance[which(df$P_combos_distance <= min(distance_combo_model[,3]))] = min(distance_combo_model[,3])
    df$P_combos_distance[which(df$P_combos_distance >= max(distance_combo_model[,3]))] = max(distance_combo_model[,3])
    df$P_IAB_depth[which(df$P_IAB_depth <= min(depth_IAB_model[, 3]))] = min(depth_IAB_model[, 3])
    df$P_IAB_depth[which(df$P_IAB_depth >= max(depth_IAB_model[, 3]))] = max(depth_IAB_model[, 3])
    df$P_combos_depth[which(df$P_combos_depth <= min(depth_combo_model[, 3]))] = min(depth_combo_model[, 3])
    df$P_combos_depth[which(df$P_combos_depth >= max(depth_combo_model[, 3]))] = max(depth_combo_model[, 3])
    
    # calculate the binomial probability
    totalcombos = sum(sumofy_dist)
    df$p_binom = (df$P_IAB_distance * df$P_IAB_depth)/
        (df$P_combos_distance * df$P_combos_depth * totalcombos)
    
    # calculate the total IABs
    totalIAB = sum(distance_IAB_model[, 2])
    
    # calculate the final interaction P values
    df$P <- apply(cbind(df$PETS, rep(totalIAB, nrow(df)), df$p_binom), 1, function(v) {
        return(1 - pbinom(q = v[1] - 1, size = v[2], prob = v[3]))
    })
    
    df$Q <- p.adjust(df$P, method = "BH")
    
    # Add to row data
    lo@rowData$mango.FDR <- df$Q
    lo@rowData$mango.P <- df$P
    
    # Scrub 
    mcols(lo@anchors) <- NULL
    
    # Filter as needed
    idxF <- df$Q <= FDR
    idxP <- df$P <= PValue
    idxA <- idxF & idxP
    
    return(subsetLoops(lo, idxA))
})