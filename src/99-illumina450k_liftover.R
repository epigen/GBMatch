#NOTE: This script was used for hg19 to hg38 liftover
create_liftover=function(){
  library(rtracklayer)
  library(data.table)
  source("http://bioconductor.org/biocLite.R")
  biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19, package = 'IlluminaHumanMethylation450kanno.ilmn12.hg19')
  
  Illumina450=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations
  hg19_probes_gr=with(Illumina450, GRanges(seqnames = Rle(chr), IRanges(start=pos, end=pos+1),strand=Rle(strand),ID=row.names(Illumina450)))
  
  chain = import.chain("/data/groups/lab_bock/shared/resources/ucsc/goldenpath/hg19/liftOver/hg19ToHg38.over.chain")
  
  hg38_probes_gr=unlist(liftOver(hg19_probes_gr,chain))
  
  hg38_probes_dt=data.table(chr=as.character(seqnames(hg38_probes_gr)),pos=as.numeric(start(hg38_probes_gr)),strand=as.character(strand(hg38_probes_gr)),ID=hg38_probes_gr$ID)
  return(hg38_probes_dt)
}

simpleCache(create_liftover(), cacheName="Illumina450k_to_hg38_forTurcanEtAl")
