library(project.init)
project.init2("GBMatch")
library(LOLA)
nenv()

eload(loadGencodeGenes("human",versNum=87))
eload(loadUCSCRepeats())
repeats = SV$repeats

eload(loadAnnotation())
BSSamples = SV$psa[file.exists(BSFile) & library=="RRBS", sample_name]; BSSamples


##functions##########
buildJadd=function(cols,funcs,special){
  r = paste("list(", paste(c(paste0(cols, "=", funcs, "(", cols, ")"),special), collapse=","), ")")
  return(r);
}

####################

setLapplyAlias(12)
simpleCache("rrbsCg", {
	sampleSummaryList = lapplyAlias(BSSamples, sampleFilter,
		minReads=10, excludeGR = repeats,
		SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
		cachePrepend="meth.cg.", cacheSubDir="meth/cg_sub", 
		readFunction=function(x) {
			tmp = BSreadBiSeq(x);
			tmp[,methyl:=round(hitCount/readCount,3)]
			return(tmp)
		}, recreate=FALSE, mc.preschedule=FALSE)
	sampleSummaryLong = rbindlist(sampleSummaryList)
	sampleSummaryLong # Cache this.
}, recreate=TRUE, noload=TRUE)

setLapplyAlias(12)
simpleCache("rrbsCgNoMin", {
  sampleSummaryList = lapplyAlias(BSSamples, sampleFilter,
  minReads=0, excludeGR = repeats,
  SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
  cachePrepend="meth.cg.", cacheSubDir="meth/cg_subNoMin", 
  readFunction=function(x) {
    tmp = BSreadBiSeq(x);
    tmp[,methyl:=round(hitCount/readCount,3)]
    return(tmp)
  }, recreate=FALSE, mc.preschedule=FALSE)
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
}, recreate=TRUE, noload=TRUE)


### TILES 
# Choose regions across which to summarize
eload(loadTiles(genomeBuild=genome, tileSize=5000))
eload(loadTiles(genomeBuild=genome, tileSize=1000))
agRegions = get(paste0("tiles5000", genome), env=SV)
agRegions1k = get(paste0("tiles1000", genome), env=SV)


cols=c("hitCount", "readCount", "methyl")
funcs = c("sum", "sum", "mean")
special=c("CpGcount=length(na.omit(methyl))")
jCommand = buildJadd(cols,funcs,special)

setLapplyAlias(12)


# 5kb tiles after subtracting repeats
simpleCache("rrbsTiled5ksub", {
  sampleSummaryList = lapplyAlias(BSSamples, sampleSummaryByRegion,
                                  regions=agRegions, excludeGR = repeats,
                                  SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
                                  cachePrepend="meth.5k.", cacheSubDir="meth/tile5k_sub", 
                                  jCommand=jCommand, 
                                  readFunction=function(x) {
                                    ino = BSreadBiSeq(x);
                                    ino[,methyl:=hitCount/readCount]
                                  }, mc.preschedule=FALSE)
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
}, recreate=TRUE, noload=TRUE)


# 1kb tiles after subtracting repeats
simpleCache("rrbsTiled1ksub", {
  sampleSummaryList = lapplyAlias(BSSamples, sampleSummaryByRegion,
                                  regions=agRegions1k, excludeGR = repeats,
                                  SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
                                  cachePrepend="meth.1k.", cacheSubDir="meth/tile1k_sub", jCommand=jCommand, 
                                  readFunction=function(x) {
                                    ino = BSreadBiSeq(x);
                                    ino[,methyl:=hitCount/readCount]
                                  }, mc.preschedule=FALSE)
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
}, recreate=TRUE, noload=TRUE)


#promoters
prom1k=promoters(SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]], upstream=1000, downstream=500)
prom1k=prom1k[!duplicated(prom1k)]

simpleCache("rrbsProm1kb", {
  sampleSummaryList = lapplyAlias(BSSamples, sampleSummaryByRegion,
                                  regions=prom1k, excludeGR = NULL,
                                  SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
                                  cachePrepend="prom.1k.", cacheSubDir="meth/prom1k", 
                                  jCommand=jCommand, 
                                  readFunction=function(x) {
                                    ino = BSreadBiSeq(x);
                                    ino[,methyl:=hitCount/readCount]
                                  }, mc.preschedule=FALSE)
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
}, recreate=TRUE, noload=TRUE)


#CGI
cpgIslands <- LOLA::getRegionSet("/data/groups/lab_bock/shared/resources/regions/LOLACore/hg38", collections="ucsc_features", "cpgIslandExt.bed")
cpgIslands=cpgIslands[!duplicated(cpgIslands)]

simpleCache("rrbsCGI", {
  sampleSummaryList = lapplyAlias(BSSamples, sampleSummaryByRegion,
                                  regions=cpgIslands, excludeGR = NULL,
                                  SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
                                  cachePrepend="CGI.", cacheSubDir="meth/CGI", 
                                  jCommand=jCommand, 
                                  readFunction=function(x) {
                                    ino = BSreadBiSeq(x);
                                    ino[,methyl:=hitCount/readCount]
                                  }, mc.preschedule=FALSE)
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
}, recreate=TRUE, noload=TRUE)

