library(project.init)
project.init2("GBMatch")
library(CopywriteR)

setwd(file.path(getOption("PROCESSED.PROJECT")))
annotation=fread("results_analysis/01.1-combined_annotation/annotation_combined.tsv")
annotation[,aligned_bam:=paste0(resolveLink("results_pipeline"),"/",N_number_seq,"/bsmap_hg38/",N_number_seq,".bam"),]


out="results_analysis/03-CopywriteR/results"
ref="results_analysis/03-CopywriteR/reference"
dir.create(out, recursive=TRUE)
dir.create(ref, recursive=TRUE)

#make input for copywriteR
suffix="single_100kb"



#use this for single
control_bam=normalizePath(list.files("merged_bams","controls.bam$",full.names=TRUE))
#primary run
#sample.control=data.frame(samples=c(control_bam,annotation[material=="FFPE"]$aligned_bam),controls=control_bam)

#validation run
sample.control=data.frame(samples=c(control_bam,annotation$aligned_bam),controls=control_bam)



#only needs to be produced once for each bin size
preCopywriteR(output.folder = file.path(ref),
              bin.size = 100000,
              ref.genome = "hg38",
              prefix = "chr")
#parameter for parallelisation
bp.param <- SnowParam(workers = 12, type = "SOCK")

#actrually run CopywrteR (--> median normalisations)
CopywriteR(sample.control = sample.control,
               destination.folder = file.path(out),
               reference.folder = file.path(ref, "hg38_100kb_chr"),
               bp.param = bp.param, keep.intermediary.files=FALSE)

#use more conservative of control and median of all samples as control
#define samples that should not be included to calculate the set median
exclude_samples=c("controls")

log2_counts=fread(paste0(file.path(out),"/CNAprofiles/log2_read_counts.igv"))
log2_counts_samples=log2_counts[,-c("Chromosome","Start","End","Feature",grep(exclude_samples,names(log2_counts),value=TRUE)),with=FALSE]

log2_counts_summary=rbindlist(apply(log2_counts_samples,1,function(x){data.table(set_median=median(x[abs(x)<100]),set_mean=mean(x[abs(x)<100]),set_sd=sd(x[abs(x)<100]))})) #don't include tiles with rediculously small values
log2_counts_summary=cbind(log2_counts[,c("Chromosome","Start","End","Feature","log2.controls.bam"),with=FALSE],log2_counts_summary)
log2_counts_summary[,set_control:=pmin(abs(log2.controls.bam),abs(set_median)),]
log2_counts_summary[,set_control:=ifelse(set_control==abs(log2.controls.bam),log2.controls.bam,set_median),]

file.rename(file.path(out,"/CNAprofiles/log2_read_counts.igv"), file.path(out,"/CNAprofiles/log2_read_counts_orig.igv"))

write.table("#track viewLimits=-3:3 graphType=heatmap color=255,0,0", file = file.path(out,"/CNAprofiles/log2_read_counts.igv"), sep = "\t", col.names=FALSE,row.names = FALSE, quote = FALSE)
write.table(log2_counts[,log2.controls.bam:=log2_counts_summary$set_control], file =file.path(out,"/CNAprofiles/log2_read_counts.igv"), sep = "\t", row.names = FALSE, quote = FALSE,append=TRUE)


save(log2_counts_summary,file = file.path(out, "CNAprofiles/log2_counts_summary.Rdata"))

#normalisation to control and plotting
plotCNA(destination.folder = file.path(out))

file.rename(file.path(out,"CNAprofiles"), paste0(file.path(out,"CNAprofiles_"),suffix))
