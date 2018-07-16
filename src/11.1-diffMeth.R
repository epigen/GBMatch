#NOTE:This script performs differential DNA methyaltion analysis between two samples (primary and recurring tumor) 
library(project.init)
project.init2("GBMatch")
library(LOLA)
library(RGenomeUtils)

nenv()

eload(loadAnnotation())
BSSamples = SV$psa[file.exists(BSFile) & library=="RRBS", sample_name]; BSSamples

annot=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
setnames(annot,"N_number_seq","sample_name")
SV$psa=merge(SV$psa[file.exists(BSFile) & library=="RRBS",],annot,by="sample_name")

#perform CpG wise differential methylation analysis and combine into meaningful regions afterwards
#use fishers exact test (no replicates)

##functions:
# speed-improved Fisher test (a version of the R function that is reduced to the bare essentials)
fast.fisher = function (x, y = NULL, workspace = 2e+05, hybrid = FALSE, control = list(), 
                        or = 1, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95, 
                        simulate.p.value = FALSE, B = 2000, cache=F) 
{
  if (nrow(x)!=2 | ncol(x)!=2) stop("Incorrect input format for fast.fisher")
  if (max(is.na(x))!=0) {
    
    RVAL <- list(p.value=NA, estimate=NA, null.value=NA, alternative=NA, method="fast.fisher", data.name=NA)
    attr(RVAL, "class") <- "htest"
    return(RVAL)
  }   
  if (cache) {
    key = paste(x,collapse="_")
    cachedResult = hashTable[[key]]
    if (!is.null(cachedResult)) {
      return(cachedResult)
    }
  }
  # ---- START: cut version of fisher.test ----
  DNAME <- deparse(substitute(x))
  METHOD <- "fast.fisher"
  nr <- nrow(x)
  nc <- ncol(x)
  PVAL <- NULL
  if ((nr == 2) && (nc == 2)) {
    m <- sum(x[, 1])
    n <- sum(x[, 2])
    k <- sum(x[1, ])
    x <- x[1, 1]
    lo <- max(0, k - n)
    hi <- min(k, m)
    NVAL <- or
    names(NVAL) <- "odds ratio"
    support <- lo:hi
    logdc <- dhyper(support, m, n, k, log = TRUE)
    dnhyper <- function(ncp) {
      d <- logdc + log(ncp) * support
      d <- exp(d - max(d))
      d/sum(d)
    }
    mnhyper <- function(ncp) {
      if (ncp == 0) 
        return(lo)
      if (ncp == Inf) 
        return(hi)
      sum(support * dnhyper(ncp))
    }
    pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
      if (ncp == 1) {
        if (upper.tail) 
          return(phyper(x - 1, m, n, k, lower.tail = FALSE))
        else return(phyper(x, m, n, k))
      }
      if (ncp == 0) {
        if (upper.tail) 
          return(as.numeric(q <= lo))
        else return(as.numeric(q >= lo))
      }
      if (ncp == Inf) {
        if (upper.tail) 
          return(as.numeric(q <= hi))
        else return(as.numeric(q >= hi))
      }
      d <- dnhyper(ncp)
      if (upper.tail) 
        sum(d[support >= q])
      else sum(d[support <= q])
    }
    if (is.null(PVAL)) {
      PVAL <- switch(alternative, less = pnhyper(x, or), 
                     greater = pnhyper(x, or, upper.tail = TRUE), 
                     two.sided = {
                       if (or == 0) 
                         as.numeric(x == lo)
                       else if (or == Inf) 
                         as.numeric(x == hi)
                       else {
                         relErr <- 1 + 10^(-7)
                         d <- dnhyper(or)
                         sum(d[d <= d[x - lo + 1] * relErr])
                       }
                     })
      RVAL <- list(p.value = PVAL)
    }
    mle <- function(x) {
      if (x == lo) 
        return(0)
      if (x == hi) 
        return(Inf)
      mu <- mnhyper(1)
      if (mu > x) 
        uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
      else if (mu < x) 
        1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                                                  1))$root
      else 1
    }
    ESTIMATE <- mle(x)
    #names(ESTIMATE) <- "odds ratio"
    RVAL <- c(RVAL, estimate = ESTIMATE, null.value = NVAL)
  }
  RVAL <- c(RVAL, alternative = alternative, method = METHOD, data.name = DNAME)
  attr(RVAL, "class") <- "htest"
  # ---- END: cut version of fisher.test ----    
  if (cache) hashTable[[key]] <<- RVAL # write to global variable
  return(RVAL)                                                                         
}

# this function is a wrapper around the normal t-test that returns NA instead of raising an error if assumptions are not met
fast.fisherFailsafe = function(...) { 
  obj = try(fast.fisher(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

# function for pairwise differential DNA methylation analysis
pairDiffmeth_fast.fisher=function(id,annotationTable,idColumn,sample_col1,sample_col2,file_col1,file_col2,readFunction,cacheSubDir,...){
  # identify the sample:
  i = which (id == annotationTable[,get(idColumn)]); 
  cacheName = paste0("diffMeth_",id,"__",annotationTable[i, get(sample_col1)],"vs",annotationTable[i, get(sample_col2)]);
  infile_s1 = annotationTable[i, get(file_col1)];
  infile_s2 = annotationTable[i, get(file_col2)];
  message(i, ": ", cacheName);
  # Produce a cache; read in the file and summarize by regions
  simpleCache(cacheName, {
    inObject_s1 = readFunction(infile_s1);
    inObject_s2 = readFunction(infile_s2); 
    inObject = merge(inObject_s1,inObject_s2,by=c("chr","start"),all=TRUE,suffixes=c(".1",".2"));
    if (NROW(inObject) == 0) { return(NULL) };
    inObject[, id:=id];
    message(i, ": now calculating p-values.");
    inObject[,p.value:=as.double(fast.fisherFailsafe(matrix(c(hitCount.1,readCount.1-hitCount.1,hitCount.2,readCount.2-hitCount.2),nrow = 2),alternative = "two.sided")),by=1:nrow(inObject)];
    inObject[,p.value.adj:=p.adjust(p.value,"fdr")]
    inObject
  },
  cacheSubDir=cacheSubDir, 
  buildEnvir=nlist(infile_s1,infile_s2, id, readFunction,i),...) # end simpleCache
  
  tryCatch( { 
  }, error = function(e) { 
    return(NULL)
  })
}

# combine p-values using a generalization of Fisher's method
rrbsExpectedTestCorrelation = 0.8
combineTests = function(pvalues,testWeights=NULL,correlated=F) {
  if (!is.null(testWeights)) {
    # check if weights are valid
    if (length(pvalues) != length(testWeights)) stop("Number of items in <pvalues> and in <testWeights> must be identical if weights are to be used")
    if (sum(is.na(testWeights))>0) stop("NA values are not permitted for the test weights")
    if (sum(testWeights<0)>0) stop("Weights must be positive")
    # standardize weights
    testWeights = testWeights/sum(testWeights)
  } else {
    # use equal weighting
    testWeights = rep(1/length(pvalues),length(pvalues))
  }
  pvalues[is.na(pvalues)] = 1
  if (correlated==F & is.null(testWeights)) {
    # use Fisher's classical method for combining p-values in the case of independence
    tcombined = sum(-2*log(pvalues))
    return(pchisq(tcombined,2*length(pvalues),lower.tail=FALSE))
  } else {
    # use the method proposed in Makambi (2003) Journal of Applied Statistics
    r = ifelse(correlated,rrbsExpectedTestCorrelation,0)
    m = length(pvalues)
    M.Fm = sum(-2*log(pvalues)*testWeights)
    var.M.Fm = 4*sum(testWeights^2)
    for (i in 1:m) {
      for (j in setdiff(1:m,i)) {
        var.M.Fm = var.M.Fm + testWeights[i]*testWeights[j]*(3.25*r+0.75*r^2)
      }
    }
    nu = 8/var.M.Fm
    tcombined = nu*M.Fm/2 # chi square test generalization of Makambi	   	   
    return(pchisq(tcombined,nu,lower.tail=FALSE))
  }
}


buildJadd=function(cols,funcs,special){
  r = paste("list(", paste(c(paste0(cols, "=", funcs, "(", cols, ")"),special), collapse=","), ")")
  return(r);
}

##---------------------------------------FUNCTIONS END-----------------------------------------##

#create pairs of samples surgery1 vs. surgery2
test_combinations=SV$psa[category=="GBMatch",list(s1=sample_name[surgery.x==1],s2=sample_name[surgery.x==2],f1=BSFile[surgery.x==1],f2=BSFile[surgery.x==2]),by=c("patientID","patID")]
test_combinations=na.omit(test_combinations)
test_combinations[,comparisonID:=paste0(patID,"_",c(1:length(s1))),by="patID"]
test_combinations[,combi:=paste0(comparisonID,"__",s1,"vs",s2),]

#create testcombinations with fresh frozen
test_combinations_frozen=SV$psa[category=="GLASS",list(s1=sample_name[surgery.x==1],s2=sample_name[surgery.x==2],f1=BSFile[surgery.x==1],f2=BSFile[surgery.x==2]),by=c("patientID","patID")]
test_combinations_frozen=na.omit(test_combinations_frozen)
test_combinations_frozen[,comparisonID:=paste0(patID,"_f",c(1:length(s1))),by="patID"]
test_combinations_frozen[,combi:=paste0(comparisonID,"__",s1,"vs",s2),]


# calculate diffmeth using simple cache
setLapplyAlias(12)
subDir="diffMeth/cg_fastFisher_1vs2"
test_combinations[,diffCache:=paste0(getOption("RCACHE.DIR"),subDir,"/diffMeth_",comparisonID,"__",s1,"vs",s2,".RData"),]

sel_comp=test_combinations$comparisonID

simpleCache("rrbsDiffmeth_cg_fastFisher", {
  sampleSummaryList = lapplyAlias(sel_comp, pairDiffmeth_fast.fisher,
                                  annotationTable=test_combinations, idColumn = "comparisonID", sample_col1="s1",
                                  sample_col2="s2",file_col1="f1",file_col2="f2", cacheSubDir=subDir, 
                                  readFunction=function(x) {
                                    tmp = BSreadBiSeq(x);
                                    tmp[,methyl:=round(hitCount/readCount,3)]
                                    return(tmp)
                                  }, recreate=FALSE, mc.preschedule=FALSE)
}, recreate=TRUE, noload=TRUE)


#diffmeth for frozen samples

subDir="diffMeth/cg_fastFisher_frozen"
test_combinations_frozen[,diffCache:=paste0(getOption("RCACHE.DIR"),subDir,"/diffMeth_",comparisonID,"__",s1,"vs",s2,".RData"),]

sel_comp=test_combinations_frozen$comparisonID

simpleCache("rrbsDiffmeth_cg_fastFisher_frozen", {
  sampleSummaryList = lapplyAlias(sel_comp, pairDiffmeth_fast.fisher,
                                  annotationTable=test_combinations_frozen, idColumn = "comparisonID", sample_col1="s1",
                                  sample_col2="s2",file_col1="f1",file_col2="f2", cacheSubDir=subDir, 
                                  readFunction=function(x) {
                                    tmp = BSreadBiSeq(x);
                                    tmp[,methyl:=round(hitCount/readCount,3)]
                                    return(tmp)
                                  }, recreate=FALSE, mc.preschedule=FALSE)
}, recreate=TRUE, noload=TRUE)


#summarize diffmeth for region sets
cols=c("hitCount.1", "readCount.1", "methyl.1","hitCount.2", "readCount.2", "methyl.2","sampleName.1","sampleName.2","id")
funcs = c("sum", "sum", "mean","sum", "sum", "mean","unique","unique","unique")
special=c("p.value.adj=combineTests(p.value.adj,testWeights=min_cov,correlated=TRUE)","CpGcount=length(na.omit(p.value.adj))")
jCommand = buildJadd(cols,funcs,special)


#summarize by promoters
source(file.path(getOption("PROJECT.DIR"),"src/99-summarize_promoters.R"))
eload(loadGencodeGenes("human",versNum=87))
SV$gencodeContainer$genesGR$gene_name=SV$gencodeContainer$genes$external_gene_name
SV$gencodeContainer$genesGR$ensG=SV$gencodeContainer$genes$ensembl_gene_id

prom1k=promoters(SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]], upstream=1000, downstream=500)
prom1k=prom1k[!duplicated(prom1k$gene_name)]

simpleCache("rrbsDiffmeth_Promoters1k", {
  sampleSummaryList = lapplyAlias(test_combinations$combi, sampleSummaryByProm,
                                  regions=prom1k, excludeGR = NULL,
                                  test_combinations, idColumn = "combi", fileColumn="diffCache", 
                                  cachePrepend="diffmeth.Prom1k.", cacheSubDir="diffMeth/Prom1k", 
                                  jCommand=jCommand, 
                                  readFunction=function(x) {
                                    load(x);
                                    ret[,min_cov:=pmin(readCount.1,readCount.2),];
                                    ret=na.omit(ret);
                                  }, mc.preschedule=FALSE)
  
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong}, # Cache this.
  recreate=TRUE, noload=TRUE)


simpleCache("rrbsDiffmeth_frozen_Promoters1k", {
  sampleSummaryList = lapplyAlias(test_combinations_frozen$combi, sampleSummaryByProm,
                                  regions=prom1k, excludeGR = NULL,
                                  test_combinations_frozen, idColumn = "combi", fileColumn="diffCache", 
                                  cachePrepend="diffmeth.Prom1k.frozen.", cacheSubDir="diffMeth/Prom1k_frozen", 
                                  jCommand=jCommand, 
                                  readFunction=function(x) {
                                    load(x);
                                    ret[,min_cov:=pmin(readCount.1,readCount.2),];
                                    ret=na.omit(ret);
                                  }, mc.preschedule=FALSE)
  
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong}, # Cache this.
  recreate=TRUE, noload=TRUE)



