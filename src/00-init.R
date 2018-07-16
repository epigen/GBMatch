#NOTE: This script contains generally used functions and settings. 
#It is used as an initiation file for every other script.

library(rtracklayer, quietly=TRUE)
library(data.table, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(reshape2, quietly=TRUE)
library(yaml, quietly=TRUE)

#load and apply configuration file
cfg = load.config("GBMatch")
cfg
nenv()
genome = cfg$genome$human
psa = fread(cfg$metadata$sample_annotation)

#set general plotting specifications
theme_set(theme_bw())
theme_update(axis.text=element_text(color="black"))
pdf.options(useDingbats=FALSE)

#path to external data
extData="/data/groups/lab_bock/jklughammer/projects/Glioblastoma_match/"

#functions to include sample numbers in ggplot plots
#use like this with ggplot: +stat_summary(fun.data=addN, geom="text", vjust=-0.5, col="blue")
addN <- function(y) 
  c(label=length(y), y=median(y))

addNmin <- function(y) 
  c(label=length(y), y=min(y))


#extracts legend from a ggplot (used to plot legend separately from plot)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


resolveLink=function(link){
  path=Sys.readlink(link)
  if (path ==""){path=link}
  return(path)
}

loadAnnotation = function() {
  psa[library%in%c("RRBS","rrbs"), BSFile := dirdata("results_pipeline/", sample_name, "/biseq_", genome, "/RRBS_cpgMethylation_", sample_name, ".bed")]
  psa[, PDRFile := dirdata("results_pipeline/", sample_name, "/pdr_", genome, "/", sample_name, ".pdr.bed")]
  psa[, NMMFile := dirdata("results_pipeline/", sample_name, "/nmm_", genome, "/", sample_name, ".nmm.bed")]
  
  psa[,file.exists(BSFile)]
  return(nlist(psa))
}

direpicomdata = function(...) {
  paste0(Sys.getenv("PROCESSED"), "epigenome_compendium/", ...)
}


