library(project.init)
project.init2("GBMatch")
library(Rsamtools)

setwd(getOption("PROCESSED.PROJECT"))
dir.create("merged_bams")


annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"/results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
annotation[,aligned_bam:=paste0(resolveLink("results_pipeline"),"/",N_number_seq,"/bsmap_hg38/",N_number_seq,".bam"),]


##create merged control bam
mergeBam(annotation[category=="control"]$aligned_bam,"merged_bams/controls.bam")

##create merged bams for controls m and f separately
controls[,try(mergeBam(aligned_bam,paste0("merged_bams/control_",sex[1],".bam"),overwrite=FALSE)),by="sex"]

