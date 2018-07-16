#NOTE: helper script to summarize differential methylation results across promoter regions
#Adapted from N. Sheffields analysis for the Ewing sarcoma DNA methylation study
BSAggregate_prom=function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, sumCols=NULL, jCommand=NULL, byRegionGroup=FALSE, keep.na=FALSE) {
  
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
    withinGroupID= unlist(sapply(regionsGRL.length, seq)),
    regionGroupID=rep(1:length(regionsGRL), regionsGRL.length),
    gene_name=regionsGR$gene_name,
    ensG=regionsGR$ensG)
  setkey(region2group, regionID)
  
  message("Finding overlaps...");
  fo = findOverlaps(bsgr[[1]], regionsGR)
  
  setkey(BSDT, chr, start)
  # Gut check:
  # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))
  
  message("Setting regionIDs...");
  BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
  BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
  #BSDT[queryHits(fo),regionID:=subjectHits(fo)]
  #if (!keep.na) {
  #	BSDT = BSDT[queryHits(fo),]
  #}
  
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
  
  # Now actually do the aggregate:
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
  # at each level: across regions, then across region sets. THis
  # doesn't give you a choice at this point. 
}



sampleSummaryByProm=function(id, regions, excludeGR, annotationTable,
                             idColumn, fileColumn, cachePrepend, cacheSubDir, 
                             jCommand, readFunction, byRegionGroup= FALSE,
                             regionsGRL.length=NULL, ...) {
  # Check 
  if (! "GRanges" %in% class(regions) & !"GRangesList" %in% class(regions)) {
    stop("regions is not a GRanges");
  }
  
  # identify the sample:
  i = which (id == annotationTable[,get(idColumn)]); 
  cacheName = paste0(cachePrepend, id);
  infile = annotationTable[i, get(fileColumn)];
  message(i, ": ", cacheName);
  # Produce a cache; read in the file and summarize by regions
  simpleCache(cacheName, {
    inObject = readFunction(infile);
    inObject = BSAggregate_prom(inObject, regions, excludeGR=excludeGR, regionsGRL.length=regionsGRL.length, jCommand=jCommand, byRegionGroup=byRegionGroup)
    inObject[, id:=id]
    inObject
  },
  cacheSubDir=cacheSubDir, 
  buildEnvir=nlist(infile, regions, jCommand, id, readFunction, excludeGR, regionsGRL.length, byRegionGroup), ...) # end simpleCache
  return(get(cacheName))
}