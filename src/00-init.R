library(rtracklayer, quietly=TRUE)
library(data.table, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(reshape2, quietly=TRUE)
library(yaml, quietly=TRUE)


cfg = load.config("GBMatch")
theme_set(theme_bw())
cfg
nenv()
genome = cfg$genome$human
psa = fread(cfg$metadata$sample_annotation)

extData="/data/groups/lab_bock/jklughammer/projects/Glioblastoma_match/"


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


loadEpigenomeCompendium = function() {
  # This merges 2 subprojects.
  
  epicomAnnotationFile = paste0(Sys.getenv("CODEBASE"), "epigenome_compendium/metadata/psa_rrbs_encode.csv")
  psaEC = fread(epicomAnnotationFile)
  esAnnotationFile = paste0(Sys.getenv("CODEBASE"), "epigenome_compendium/metadata/psa_rrbs_diverse.csv")
  psaES = fread(esAnnotationFile)
  psa=rbind(psaES, psaEC, fill=TRUE)
  psa[, BSFile := direpicomdata("results_pipeline/", sample_name, "/biseq_", genome, "/RRBS_cpgMethylation_", sample_name, ".bed")]
  psa[, NMMFile := direpicomdata("results_pipeline/", sample_name, "/nmm_", genome, "/", sample_name, ".nmm.bed")]
  psa[, PDRFile := direpicomdata("results_pipeline/", sample_name, "/prd_", genome, "/", sample_name, ".pdr.bed")]
  psa[, PDRFileK := direpicomdata("results_pipeline/", sample_name, "/pdr_", genome, "/", sample_name, ".pdr")]
  psa[cell_type == "EWS", cell_type:="EWS_line"]  # to match the ewing project
  
  
  psa[,]
  
  psa[, file.exists(BSFile)]
  if(any(duplicated(psa$sample_name))) {
    warning("duplicate sample names")
    which(duplicated(psa$sample_name))
  }
  
  epicom = psa
  setkey(epicom, sample_name)
  return(nlist(epicom))
}


loadCancerCompendium = function() {
  epicomAnnotationFile = paste0(Sys.getenv("CODEBASE"), "epigenome_compendium/metadata/psa_rrbs_intracancer.csv")
  psa = fread(epicomAnnotationFile)
  
  # For the final analysis, I remove the RMS samples because I am not confident in the data quality.
  psa = psa[!cell_type %in% c("RMS"), ]
  
  psa[, BSFile := direpicomdata("results_pipeline/", sample_name, "/biseq_", genome, "/RRBS_cpgMethylation_", sample_name, ".bed")]
  psa[, NMMFile := direpicomdata("results_pipeline/", sample_name, "/nmm_", genome, "/", sample_name, ".nmm.bed")]
  psa[, PDRFile := direpicomdata("results_pipeline/", sample_name, "/pdr_", genome, "/", sample_name, ".pdr.bed")]
  # Kendell's PDR file
  psa[, PDRFileK := direpicomdata("results_pipeline/", sample_name, "/pdr_", genome, "/", sample_name, ".pdr")]
  psa[, file.exists(BSFile)]
  if(any(duplicated(psa$sample_name))) {
    warning("duplicate sample names")
    which(duplicated(psa$sample_name))
  }
  
  cancerData = psa
  setkey(cancerData, sample_name)
  return(nlist(cancerData))
}



#' Loads the Relative Proportion of Intermediate Methyation; you must
#' first have created the cache for this.
loadRPIM = function() {
  loadCaches("rrbsGBMatchRPIM")
  loadCaches("rrbsGBMatchPIM")
  # The average RPIM is a score that is only relevant *in context* -- 
  # of which samples are included in the matrix.
  # I used to apply across columns (2) here, but that gives you the inverse
  # relationship (using the sample of interest on the denominator of the
  # log ratio test). This results in the negative of the true RPIM.
  # I switched back to 1 now.
  RPIM = apply(rrbsGBMatchRPIM, 1, mean)
  PIM = rrbsGBMatchPIM
  dfPIM = data.table(pim=as.numeric(PIM), id=names(PIM))
  setkey(dfPIM, id)
  
  # Convert to data frame
  dfRPIM = data.table(rpim=as.numeric(RPIM), id=names(RPIM))
  setkey(dfRPIM, id)
  # make sure the order matches
  stopifnot(all(names(SV$PIM) == names(SV$RPIM)))
  return(nlist(RPIM, dfRPIM, PIM, dfPIM))
}

