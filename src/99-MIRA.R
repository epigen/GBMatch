#NOTE: Slightly adapted from N. Sheffields analysis for the DNA methylation Ewing sarcoma paper
#' Methylation-based Inference of Regulatory Activity (MIRA)
#' MIRA is a score that measures the degree of dip in DNA methylation
#' level surrounding a regulatory site of interest, such as a
#' transcription factor binding sites.
#' This script (an incipient R package) provides functions
#' for aggregating methylation data across region sets, in bins


#' A function to score the degree of "dip" in DNA methylation data
#' Methylation-based Inference of Regulatory Activity (MIRA)
#' calculates the ratio between the flank (shoulders) and midpoint.
#'
#' @param flankShift The number of bins away from the center to use as the
#' flank.
scoreDip = function(values, binCount, flankShift = 5) {
  centerSpot = ceiling(binCount/2)
  leftSide = centerSpot - flankShift
  rightSide = centerSpot + flankShift
  midpoint = (values[centerSpot] + values[centerSpot+1] + values[centerSpot-1] ) /3
  # log ratio...
  return ( log ( ((values[leftSide] + values[rightSide])/2) / midpoint ) )
}


#' Convenience function to aggregate bisulfite scores across bins, across regions.
#' You need a list of ranges (a bed file as a data.table, in a list)
#' and a list of BSDTs.

#' @param rangeDTList This is a list of data.tables, with each just having chr,
#' start and end.
#' @param anno An annotation table (data.table), keyed with the same key as the
#' BSDTs; also with a `cell_type` column, used to append cell_type to the binned
#' result tables.
binProcess = function(rangeDTList, BSDTNames, anno, binCount=11) {
  # aggregate in bins
  if (! "list" %in% class(BSDTNames)) {
    BSDTNames = list(BSDTNames)
  }
  if (! "list" %in% class(rangeDTList)) {
    rangeDTList = list(rangeDTList)
  }
  
  binResults = list()
  # Aggregate each Range into bins.
  
  for (DTName in names(rangeDTList)) {
    if (! "data.table" %in% class(get(DTName)) ) { message("ranges must be a DT") }
    message(DTName)
    
    binnedBSDTList = list()
    for (i in 1:length(BSDTNames)) {
      BSDTName = BSDTNames[[i]]
      message("BSDT ", i, ": ", BSDTName)
      
      cacheName = paste0("binnedBSDT_", DTName, "_", BSDTName, "_", binCount)
      rangeDT = rangeDTList[[DTName]]
      binnedBSDTList[[i]] = simpleCache(cacheName , {
        loadCaches(BSDTName)
        BSDT = get(BSDTName)
        BSBinAggregate(BSDT, rangeDT, binCount)
      },buildEnvir=nlist(BSDTName, rangeDT, binCount), cacheSubDir="binnedBSDT")
    } # end for binned BSDT list
    binnedBSDT = rbindlist(binnedBSDTList)
    setkey(binnedBSDT, id)
    binnedBSDT=binnedBSDT[anno]
    
    # clear any NA cell_types
    binnedBSDT = normalizeEwingsToNonEwings(binnedBSDT)
    
    # Calculate Dipscores (MIRA)
    dipScores = binnedBSDT[,list(cell_type=unique(cell_type), dipScore=unique(scoreDip(methyl, binCount, flankShift = 5))), by=id]
    setkey(dipScores,id)
    dipScores=dipScores[anno]
    
    binResults[[DTName]] = nlist(binnedBSDT, dipScores)
    
  } # end for ranges
  return(binResults)
}

#' Aggregating signals in bins across a set of regions
#'
#' given a start, end, and number of bins, to divide,
#' this function will split the regions into bins.
#' Bins will be only approximately the same size, due to rounding.
#' (they should not be more than 1 different).
#'
#' Use case: take a set of regions, like CG islands, and bin them; now you can
#' aggregate signal scores across the bins, giving you an aggregate signal
#' in bins across many regions of the same type.
#'
#' In theory, this just runs on 3 values, but you can run it inside a
#' data.table j expression to divide a bunch of regions in the same way.
#' @param start
#'
#' @return
#' A data.table, expanded to nrow= number of bins, with these id columns:
#' 		id: region ID
#' 		binID: repeating ID (this is the value to aggregate across)
#' 		ubinID: unique bin IDs
#' @examples
#' loadCGData("hg19")
#' cgIslandsDT = data.table(...)
#' binnedCGI = cgIslandsDT[, binRegion(start,end, 50)]
binRegion = function(start, end, bins, idDF=NULL) {
  binSize = (end-start)/(bins)
  breaks = round(rep(start, each=(bins+1)) + (0:(bins)) * rep(binSize, each=(bins+1)))
  
  endpoints = (bins+1) * (1:(length(start)))
  startpoints = 1 + (bins+1)  * (0:(length(start)-1))
  #TODO: remove this code split
  if (is.null(idDF)) {
    dt = data.table(start=breaks[-endpoints], end=breaks[-startpoints], id=rep((1:length(start)), each=bins), binID= 1:bins, ubinID=1:length(breaks[-startpoints]), key="id")
  } else {
    chr = rep(idDF, each=bins)
    dt = data.table(chr, start=breaks[-endpoints], end=breaks[-startpoints], id=rep((1:length(start)), each=bins), binID= 1:bins, ubinID=1:length(breaks[-startpoints]), key="id")
    
  }
  return(dt)
}

#' A wrapper of BSAggregate that first bins regions and then aggregates
#' each bin across a set of regions, individually.
#'
#' Produced originally for binning Ewing RRBS data across various region sets
#'
#' @param rangeDT A data table with the sets of regions to be binned,
#' with columns named start, end
#' @param binCount Number of bins across the region
#' @param byRegionGroup Pass along to binCount (see ?binCount)
#' @param minReads Filter out bins with fewer than X reads before returning.
BSBinAggregate = function(BSDT, rangeDT, binCount, minReads = 500, byRegionGroup=TRUE) {
  if (! "data.table" %in% class(rangeDT)) {
    stop("rangeDT must be a data.table")
  }
  seqnamesColName = "seqnames"  # default column name
  if (! "seqnames" %in% colnames(rangeDT)) {
    if ("chr" %in% colnames(rangeDT)) {
      message("seqnames column name set to: chr")
      seqnamesColName = "chr"
    } else {
      # Got neither.
      stop("rangeDT must have a seqnames column")
    }
  }
  
  message("Binning...")
  binnedDT = rangeDT[, binRegion(start, end, binCount, get(seqnamesColName))]
  binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
  message("Aggregating...")
  binnedBSDT = BSAggregate(BSDT, regionsGRL=GRangesList(binnedGR), jCommand=buildJ(c("methyl", "readCount"), c("mean", "sum")), byRegionGroup=byRegionGroup, splitFactor="id")
  # If we aren't aggregating by bin, then don't restrict to min reads!
  if (byRegionGroup) {
    binnedBSDT = binnedBSDT[readCount > minReads,]
  }
  return(binnedBSDT)
}

#' BSaggregate -- Aggregate a BSDT across regions or region groups,
#' for multiple samples at a time.
#' This function is as BScombineByRegion, but can handle not only multiple
#' samples in BSDT, but also simultaneously multiple region sets by passing
#' a regionsGRL (GRangesList object).
#' you can use jCommand to do other functions.

#' Given a bisulfite data table as input, with an identifier column for
#' different samples; plus a GRanges objects with regions to aggregate.
#'
#' @param BSDT The bisulfite data.table (output from one of the parsing
#' functions for DNA methylation calls) that you wish to aggregate. It can
#' be a combined table, with individual samples identified by column passed
#' to splitFactor.
#' @param regionsGRL Regions across which you want to aggregate.
#' @param excludeGR A GenomicRanges object with regions you want to
#' exclude from the aggregation function. These regions will be eliminated
#' from the input table and not counted.
#' @param jCommand You can pass a custom command in the j slot to data.table
#' specifying which columns to aggregate, and which functions to use. You
#' can use buildJ() to build a jCommand argument easily.
#' @param byRegionGroup You can aggregate by regionID or by regionGroupID;
#' this reflects the regionsGRL that you pass; by default, BSAggregate will
#' aggregate each region individually -- scores will then be contiguous, and
#' the output is 1 row per region.
#' Turn on this flag to aggregate across all region groups, making the result
#' uncontiguous, and resulting in 1 row per *region group*.
#'
#' @export
BSAggregate = function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, sumCols=NULL, jCommand=NULL, byRegionGroup=FALSE, keep.na=FALSE) {
  
  # Assert that regionsGRL is a GRL.
  # If regionsGRL is given as a GRanges, we convert to GRL
  if( "GRanges" %in% class(regionsGRL)) {
    regionsGRL = GRangesList(regionsGRL);
  } else if (! "GRangesList" %in% class(regionsGRL)) {
    stop("regionsGRL is not a GRanges or GRangesList object");
  }
  
  if(! is.null(excludeGR)) {
    BSDT = BSFilter(BSDT, minReads=0, excludeGR)
  }
  
  bsgr = BSdtToGRanges(list(BSDT));
  
  additionalColNames = setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor));
  
  colModes = sapply(BSDT,mode);
  if (is.null(sumCols)) {
    sumCols = setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
    # Restrict to numeric columns.
    sumCols = intersect(sumCols, names(colModes[which(colModes == "numeric")]))
    
  }
  # It's required to do a findoverlaps on each region individually,
  # Not on a GRL, because of the way overlaps with GRLs work. So,
  # we must convert the GRL to a GR, but we must keep track of which
  # regions came from which group.
  regionsGR = unlist(regionsGRL)
  
  if(is.null(regionsGRL.length)) {
    if (length(regionsGRL) > 100) {
      message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
    }
    regionsGRL.length = sapply(regionsGRL, length)
    message("Done counting regionsGRL lengths.");
  }
  
  # Build a table to keep track of which regions belong to which group
  region2group = data.table(
    regionID=1:length(regionsGR),
    chr=as.vector(seqnames(regionsGR)),
    start=as.vector(start(regionsGR)),
    end=as.vector(end(regionsGR)),
    withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
    regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
  setkey(region2group, regionID)
  
  
  message("Finding overlaps...");
  fo = findOverlaps(bsgr[[1]], regionsGR)
  
  setkey(BSDT, chr, start)
  # Gut check:
  # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))
  
  message("Setting regionIDs...");
  BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
  
  if (NROW(BSDT) < 1) {
    warning("No BSDT sites in the given region list; please expand your regionsGRL")
    return(NULL)
  }
  
  BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
  
  if (is.null(jCommand)) {
    cols=c(sumCols, keepCols)
    funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
    jCommand = buildJ(cols, funcs)
  }
  message("jCommand: ", jCommand)
  
  # Define aggregation column. aggregate by region or by region group?
  if (byRegionGroup) {
    agCol = "regionGroupID";
  } else {
    agCol = "regionID"; # Default
  }
  
  # Build the by string
  if (is.null(splitFactor)) {
    byString = paste0("list(regionID)");
  } else {
    byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
  }
  
  # Now actually do the aggregation:
  message("Combining...");
  bsCombined = BSDT[,eval(parse(text=jCommand)), by=eval(parse(text=byString))]
  setkey(bsCombined, regionID)
  # Now aggregate across groups.
  # I do this in 2 steps to avoid assigning regions to groups,
  # which takes awhile. I think this preserve memory and is faster.
  
  # Define aggregation column. aggregate by region or by region group?
  if (byRegionGroup) {
    # must set allow=TRUE here in case there are multiple IDs (splitCol)
    bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
    if (! is.null(splitFactor) ) {
      byStringGroup = paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
    } else {
      byStringGroup = "list(regionGroupID)"
    }
    bsCombined=bsCombined[,eval(parse(text=jCommand)), by=eval(parse(text=byStringGroup))]
    return(bsCombined);
  } else {
    e = region2group[bsCombined,]
    setkey(e, regionID);
    return(e);
  }
  # WARNING: There are now 2^2 ways to aggregate, sum vs mean
  # at each level: across regions, then across region sets. This
  # doesn't give you a choice at this point.
}

#' Given a binned BSDT (which is produced from BSBinAggregate), this plot will
#' produce and return a few ggplots (which could then be printed).
#' It produces a combined, and a faceted by cell-type plot, showing how signal
#' varies across bins.
BSBinPlots = function(binnedBSDT, binCount, regionName="regions") {
  byclass = binnedBSDT[, list(methyl= mean(methyl), methylNormRel=mean(methylNormRel)), by=list(regionGroupID, grepl("EWS", cell_type))]
  
  combinedPlot = ggplot(binnedBSDT, aes(y=log2(methylNormRel), x=factor(regionGroupID), color=grepl("EWS", cell_type))) + ylim(c(-1,1.25)) + theme_ns() + geom_hline(aes(yintercept=0), alpha=.2) + geom_point(position="jitter", alpha=.3)  + geom_line(data=byclass, aes(y=log2(methylNormRel), x=regionGroupID, group=grepl, color=grepl), size=1) + ylab("Methylation log_2 ratio vs normal") + xlab(paste0("Genome bins surrounding ", regionName, " sites")) + theme(legend.position=c(1,0), legend.justification=c(1,0))+ labs(colour = "EWS") +scale_colour_brewer(palette="Set1") + scale_x_discrete(labels=labelBinCenter(binCount))
  
  inclassPlot = ggplot(binnedBSDT[grepl("EWS", cell_type),], aes(y=log2(methylNormRel), x=factor(regionGroupID), color=grepl("EWS", cell_type))) + ylim(c(-1, 1)) + theme_ns() + geom_hline(aes(yintercept=0), alpha=checkTransparency(.2)) + geom_point(position="jitter", alpha=checkTransparency(.3))  + geom_smooth(data=byclass[grepl==TRUE,], aes(y=log2(methylNormRel), x=regionGroupID, group=grepl, color=grepl), size=1, fill="red") + ylab("Methylation log_2 ratio vs normal") + xlab(paste0("Genome bins surrounding ", regionName, " sites")) + theme(legend.position=c(1,0), legend.justification=c(1,0))+ labs(colour = "EWS") +scale_colour_brewer(palette="Set1") + scale_x_discrete(labels=labelBinCenter(binCount))
  
  bytype = binnedBSDT[, list(methyl= mean(methyl), methylNormRel=median(methylNormRel)), by=list(regionGroupID, cell_type)]
  
  facetPlot = ggplot(bytype, aes(y=log(methylNormRel), x=factor(regionGroupID), group=cell_type))  + theme_classic() + geom_hline(aes(yintercept=0), colour="#aaaaaa", linetype="dashed")  + geom_line(alpha=checkTransparency(.8), size=1)  + ylab("Methylation log ratio vs normal") + xlab(paste0("Genome bins surrounding ", regionName, " sites")) + theme(legend.position=c(1,0), legend.justification=c(1,0))+ labs(colour = "EWS") +scale_colour_brewer(palette="Set1")  + facet_wrap(~ cell_type) + geom_hline(aes(yintercept=0), colour="#999999", linetype="dashed") + scale_x_discrete(breaks=seq_len(binCount),labels=labelBinCenter(binCount))  + theme_blank_facet_label()
  
  return(nlist(inclassPlot, combinedPlot, facetPlot))
}


#' Calculate the median DNA methylation in each bin, for all non-EWS samples, and
#' then divide DNA methylation by this median value to get a relative increase/decrease
#' in each bin.
normalizeEwingsToNonEwings = function(binnedBSDT) {
  #medMethyl = binnedBSDT[, list(medMethyl= median(methyl)), by=list(regionGroupID)]
  
  # Normalize to MSCs:

  
  # Normalize to control:
  medMethyl = binnedBSDT[grepl("control", category), list(medMethyl= mean(methyl)), by=list(regionGroupID)]
  
  binnedBSDT = merge(binnedBSDT, medMethyl, by.x="regionGroupID", by.y="regionGroupID")
  
  binnedBSDT[, methylNormRel:=methyl/medMethyl]
  #	binnedBSDT[, methylNormRel:=methyl/binnedBSDT[,median(methyl)]]
  print(binnedBSDT) # data.table quirk requires print to print on 2nd call.
  return(binnedBSDT)
}
