library(project.init)
project.init2("GBMatch")
library(LOLA)
nenv()
eload(loadAnnotation())
eload(loadGencodeGenes("human",versNum=87))
eload(loadUCSCRepeats())
repeats = SV$repeats

# Choose regions across which to summarize
eload(loadTiles(genomeBuild=genome, tileSize=5000))
eload(loadTiles(genomeBuild=genome, tileSize=1000))

agRegions = get(paste0("tiles5000", genome), env=SV)
agRegions1k = get(paste0("tiles1000", genome), env=SV)
prom1k=promoters(SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]], upstream=1000, downstream=500)
prom1k=prom1k[!duplicated(prom1k)]

##functions##########
buildJadd=function(cols,funcs,special){
  r = paste("list(", paste(c(paste0(cols, "=", funcs, "(", cols, ")"),special), collapse=","), ")")
  return(r);
}

####################


# Now, aggregate PDR across the same region sets.

# Build summary J command
cols=c("ConcordantMethylatedReadCount", "ConcordantUnmethylatedReadCount", "DiscordantReadCount", "PDRa")
funcs = c("sum", "sum", "sum", "mean")
special=c("CpGcount=length(na.omit(PDRa))")
jCommand = buildJadd(cols,funcs,special)

# Choose samples on which to run this summary:
pdrSamples = SV$psa[file.exists(PDRFile), sample_name]; pdrSamples



# This lapply will run on every sample in pdrSamples (defined above)
# A summary by region; It just reads in the file, and then summarizes the
# columns (defined by the jCommand) across the regions defined in 'regions'
setLapplyAlias(12)

simpleCache("pdrTiled5ksub", {
  sampleSummaryList = lapplyAlias(pdrSamples, sampleSummaryByRegion,
                                  regions=agRegions, excludeGR = repeats,
                                  SV$psa, idColumn = "sample_name", fileColumn="PDRFile", 
                                  cachePrepend="pdr.5k.", cacheSubDir="pdr/tile5k_sub", jCommand, 
                                  readFunction=function(x) { 
                                    ino = fread(x);
                                    ino[, PDRa := DiscordantReadCount/(ConcordantMethylatedReadCount+ConcordantUnmethylatedReadCount+DiscordantReadCount)];
                                  }, mc.preschedule=FALSE)
  
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
},recreate=TRUE)


#cache 1kb tiles
simpleCache("pdrTiled1ksub", {
  sampleSummaryList = lapplyAlias(pdrSamples, sampleSummaryByRegion,
                                  regions=agRegions1k, excludeGR = repeats,
                                  SV$psa, idColumn = "sample_name", fileColumn="PDRFile", 
                                  cachePrepend="pdr.1k.", cacheSubDir="pdr/tile1k_sub", jCommand, 
                                  readFunction=function(x) { 
                                    ino = fread(x);
                                    ino[, PDRa := DiscordantReadCount/(ConcordantMethylatedReadCount+ConcordantUnmethylatedReadCount+DiscordantReadCount)];
                                  }, mc.preschedule=FALSE)
  
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
},recreate=TRUE)



#summarize promoters
simpleCache("pdrPromoters1ksub", {
  sampleSummaryList = lapplyAlias(pdrSamples, sampleSummaryByRegion,
                                  regions=prom1k, excludeGR = repeats,
                                  SV$psa, idColumn = "sample_name", fileColumn="PDRFile", 
                                  cachePrepend="pdr.prom1k.", cacheSubDir="pdr/prom1k_sub", jCommand, 
                                  readFunction=function(x) { 
                                    ino = fread(x);
                                    ino[, PDRa := DiscordantReadCount/(ConcordantMethylatedReadCount+ConcordantUnmethylatedReadCount+DiscordantReadCount)];
                                  }, mc.preschedule=FALSE)
  
  sampleSummaryLong = rbindlist(sampleSummaryList)
  sampleSummaryLong # Cache this.
},recreate=TRUE)



#per CGI
cpgIslands <- LOLA::getRegionSet("/data/groups/lab_bock/shared/resources/regions/LOLACore/hg38", collections="ucsc_features", "cpgIslandExt.bed")
cpgIslands=unlist(cpgIslands[!duplicated(cpgIslands)])


simpleCache("pdrCgi", {
	sampleSummaryList = lapplyAlias(pdrSamples, sampleSummaryByRegion,
		regions=cpgIslands,excludeGR=NULL, 
		SV$psa, idColumn = "sample_name", fileColumn="PDRFile", 
		cachePrepend="pdr.cgi", cacheSubDir="pdr/cgi", jCommand, 
		readFunction=function(x) { 
			ino = fread(x);
			ino[, PDRa := DiscordantReadCount/(ConcordantMethylatedReadCount+ConcordantUnmethylatedReadCount+DiscordantReadCount)];
		}, mc.preschedule=FALSE)

	sampleSummaryLong = rbindlist(sampleSummaryList)
	sampleSummaryLong # Cache this.
},recreate=TRUE)

combinePDRCache("pdrCgi")


