library(project.init)
project.init2("GBMatch")
library(gsubfn)
library(gtools)
library(ggrepel)


CNA_dir="CNAprofiles_single_100kb_missingSamples"
wd=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/03-CopywriteR/results",CNA_dir,"validation")
dir.create(wd)
setwd(wd)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

##get cytoband
load(file="../summary/cytoband_annot.RData")

##get WGS CNV calls ("gold standard")
WGS_CNV_dir=file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/cnvkit_vs_pga")
WGS_CNV_files=list.files(WGS_CNV_dir,pattern="call.cns",full.names=TRUE)

CNV_all=data.table()
for (path in WGS_CNV_files){
  sample=gsub(".*/|_S[0-9]+.*","",path)
  CNV=fread(path,drop="gene")
  CNV[,sample:=sample]
  CNV_all=rbindlist(list(CNV_all,CNV))
}

CNV_all[,seg_width:=end-start,]

#annotate with cytoband (chromosome arms)
CNV_all_gr=with(CNV_all,GRanges(seqnames = Rle(chromosome), IRanges(start=start, end=end),strand=Rle("*"),sample=sample))
ov=as.data.frame(findOverlaps(CNV_all_gr,cytoband_annot$gr,type="any",select="all"))
CNV_all_cb=unique(cbind(CNV_all[ov$queryHits],cytoband_annot$dt[ov$subjectHits,c("chromArm","chromArm_length"),with=FALSE]))


CNV_all_cb[,significant:=ifelse(abs(log2)>=0.5,TRUE,FALSE),]
CNV_all_cb[,variant:=ifelse(significant==TRUE & log2>0,"amplification",ifelse(significant==TRUE & log2<0,"deletion","nc")),]
CNV_all_red=CNV_all_cb[significant==TRUE,list(total_size=sum(seg_width),mean_log2=mean(log2)),by=c("chromosome","chromArm","chromArm_length","variant","sample")]
CNV_all_red[,variant_ratio:=total_size/chromArm_length,]
CNV_all_red[variant_ratio>0.6]


#plot
ggplot(CNV_all_cb[sample=="N1020_10"],aes(y=log2,yend=log2))+geom_segment(aes(x = start, xend = end,col=variant))+facet_grid(~chromosome,scales="free_x")+xlab("")+ylab("")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+scale_color_manual(values=c("amplification"="red","deletion"="green","nc"="black"))

