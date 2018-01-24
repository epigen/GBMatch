library(project.init)
project.init2("GBMatch")
library(gsubfn)
library(gtools)
library(ggrepel)
library(pROC)



CNA_dir="CNAprofiles_single_100kb"
wd=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/03-CopywriteR/results",CNA_dir,"validation")
dir.create(wd)
setwd(wd)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

##get cytoband
load(file="../summary/cytoband_annot.RData")
cytoband_red=unique(cytoband_annot$dt[,list(chromStart=min(chromStart),chromEnd=max(chromEnd)),c("chrom","chromArm","chromArm_length"),])
setnames(cytoband_red,"chrom","chromosome")
cytoband_red_gr=with(cytoband_red,GRanges(seqnames = Rle(chromosome), IRanges(start=chromStart, end=chromEnd),strand=Rle("*")))


#get RRBS CNV segments and annotate with cytoband
load(file="../summary/segments_ctr_annot.RData")
ov_RRBS_seg=as.data.frame(findOverlaps(segments_ctr_annot$gr,cytoband_red_gr,type="any",select="all"))
RRBS_segments_all_cb=unique(cbind(segments_ctr_annot$dt[ov_RRBS_seg$queryHits],cytoband_red[ov_RRBS_seg$subjectHits,c("chromArm","chromStart","chromEnd","chromArm_length"),with=FALSE]))
RRBS_segments_cb=RRBS_segments_all_cb[category=="GBmatch_valF"]
setnames(RRBS_segments_cb,c("chrom","seg.mean"),c("chromosome","log2"))
RRBS_segments_cb[,frag_len_RRBS:=loc.end-loc.start,by=c("sample","chromosome","chromArm","loc.start","loc.end")]
RRBS_segments_cb[,frag_len_corr_RRBS:=min(loc.end,chromEnd)-max(loc.start,chromStart),by=c("sample","chromosome","chromArm","loc.start","loc.end")]
RRBS_segments_cb[,variant:=ifelse(log2<0,"deletion",ifelse(log2>0,"amplification","nc")),]
RRBS_segments_cb_red=RRBS_segments_cb[,list(frag_len_corr_RRBS=sum(frag_len_corr_RRBS),log2_mean_RRBS=weighted.mean(log2,frag_len_corr_RRBS)),by=c("sample","chromosome","chromArm","chromArm_length","variant")]



##get WGS CNV calls
WGS_CNV_dir=file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/cnvkit_vs_pga")
WGS_CNV_files=list.files(WGS_CNV_dir,pattern="call.cns",full.names=TRUE)

WGS_segments_all=data.table()
for (path in WGS_CNV_files){
  sample=gsub(".*/|_S[0-9]+.*","",path)
  CNV=fread(path,drop="gene")
  CNV[,sample:=sample]
  WGS_segments_all=rbindlist(list(WGS_segments_all,CNV))
}


#annotate WGS CNV with cytoband (chromosome arms)
WGS_segments_all_gr=with(WGS_segments_all,GRanges(seqnames = Rle(chromosome), IRanges(start=start, end=end),strand=Rle("*"),sample=sample))
ov_CNV_seg=as.data.frame(findOverlaps(WGS_segments_all_gr,cytoband_red_gr,type="any",select="all"))
WGS_segments_cb=unique(cbind(WGS_segments_all[ov_CNV_seg$queryHits],cytoband_red[ov_CNV_seg$subjectHits,c("chromArm","chromStart","chromEnd","chromArm_length"),with=FALSE]))

WGS_segments_cb[,frag_len_WGS:=end-start,by=c("sample","chromosome","chromArm","start","end")]
WGS_segments_cb[,frag_len_corr_WGS:=min(end,chromEnd)-max(start,chromStart),by=c("sample","chromosome","chromArm","start","end")]
WGS_segments_cb[,variant:=ifelse(log2<0,"deletion",ifelse(log2>0,"amplification","nc")),]
WGS_segments_cb_red=WGS_segments_cb[,list(frag_len_corr_WGS=sum(frag_len_corr_WGS),log2_mean_WGS=weighted.mean(log2,frag_len_corr_WGS)),by=c("sample","chromosome","chromArm","chromArm_length","variant")]
WGS_segments_cb_red[,sample:=paste0(sample,"_fv")]


#combine RRBS and WGS results on the basis of chromosome arms and variant
segments_cb_combined=merge(WGS_segments_cb_red,RRBS_segments_cb_red[sample%in%WGS_segments_cb_red$sample],by=c("sample","chromosome","chromArm","chromArm_length","variant"),all=TRUE)
segments_cb_combined[is.na(segments_cb_combined)]=0
segments_cb_combined[,c("var_fraction_WGS","var_fraction_RRBS"):=list(frag_len_corr_WGS/chromArm_length,frag_len_corr_RRBS/chromArm_length),]
segments_cb_combined[,var_fraction_max:=pmax(var_fraction_WGS,var_fraction_RRBS),]
segments_cb_combined[,chrom_correlation:=round(cor(log2_mean_WGS,log2_mean_RRBS),3),by="chromosome"]
segments_cb_combined[,chromosome_annot:=factor(paste0("chr ",chromosome),levels=paste0("chr ",c(as.character(1:22),"X","Y"))),]

pdf("chromArm_foldChange_cor.pdf",height=8,width=15)
ggplot(segments_cb_combined,aes(x=log2_mean_WGS,y=log2_mean_RRBS,col=chromArm,fill=chromArm,size=var_fraction_max))+geom_point(shape=21,alpha=0.4)+facet_wrap(~chromosome_annot+paste0("r=",chrom_correlation),scale="free",nrow=4)+scale_size_continuous(range=c(1,3))
dev.off()


##ROC analysis

# for selected variants and significance thresholds
segments_cb_combined[,c("sig_0.1_WGS","sig_0.3_WGS","sig_0.5_WGS","sig_0.1_RRBS","sig_0.3_RRBS","sig_0.5_RRBS"):=list(abs(log2_mean_WGS)>=0.1,abs(log2_mean_WGS)>=0.3,abs(log2_mean_WGS)>=0.5,abs(log2_mean_RRBS)>=0.1,abs(log2_mean_RRBS)>=0.3,abs(log2_mean_RRBS)>=0.5),]

sub=segments_cb_combined[chromosome=="19"&variant=="deletion"]

plot(roc(sub$sig_0.1_WGS,sub$log2_mean_RRBS))
plot(roc(sub$sig_0.3_WGS,sub$log2_mean_RRBS))
plot(roc(sub$sig_0.5_WGS,sub$log2_mean_RRBS))

#more automated way

create_roc_family=function(data,chrom,var){
  sub=data[chromosome==chrom&variant==var]
  range=with(sub,c(min(abs(log2_mean_WGS[log2_mean_WGS!=0])),max(abs(log2_mean_WGS))))
  sig_thres=seq(from=range[1],to=range[2],length.out=10)
  
  exp=as.data.table(expand.grid(sample=unique(sub$sample),chromArm=unique(sub$chromArm)))
  sub_exp=merge(sub,exp,by=c("sample","chromArm"),all=TRUE)
  sub_exp[is.na(log2_mean_WGS),log2_mean_WGS:=0]
  sub_exp[is.na(log2_mean_RRBS),log2_mean_RRBS:=0]
  
  all_ROCs=data.table()
  for (thres in sig_thres){
    print(thres)
    ROC=roc(as.numeric(abs(sub_exp$log2_mean_WGS)>=thres),sub_exp$log2_mean_RRBS)
    ROC_dt=data.table(Sensitivity=ROC$sensitivities,`1-Specificity`=1-ROC$specificities,AUC=ROC$auc,N_true=length(ROC$cases),N_false=length(ROC$controls),sig_thres=thres,variant=paste0("chrom ",chrom," ",var))
    ROC_dt=ROC_dt[order(c(Sensitivity))]
    all_ROCs=rbindlist(list(all_ROCs,ROC_dt))
  }
return(all_ROCs)

}

ROC_chr10_del=create_roc_family(data=segments_cb_combined,chrom="10",var="deletion")
ROC_chr7_ampl=create_roc_family(data=segments_cb_combined,chrom="7",var="amplification")
ROC_chr19_del=create_roc_family(data=segments_cb_combined,chrom="19",var="deletion")
ROC_chr1_del=create_roc_family(data=segments_cb_combined,chrom="1",var="deletion")


pdf("chromArm_foldChange_ROC.pdf",height=6,width=11)
ggplot(ROC_chr10_del,aes(x=`1-Specificity`,y=Sensitivity,col=sig_thres))+geom_line()+facet_wrap(~paste0("log2FC=",round(sig_thres,3))+paste0("AUC=",round(AUC,3),"\nture=",N_true," false=",N_false),nrow=2)+ggtitle("Chromosome 10 deletion")

ggplot(ROC_chr7_ampl,aes(x=`1-Specificity`,y=Sensitivity,col=sig_thres))+geom_line()+facet_wrap(~paste0("log2FC=",round(sig_thres,3))+paste0("AUC=",round(AUC,3),"\nture=",N_true," false=",N_false),nrow=2)+ggtitle("Chromosome 7 amplification")

ggplot(ROC_chr19_del,aes(x=`1-Specificity`,y=Sensitivity,col=sig_thres))+geom_line()+facet_wrap(~paste0("log2FC=",round(sig_thres,3))+paste0("AUC=",round(AUC,3),"\nture=",N_true," false=",N_false),nrow=2)+ggtitle("Chromosome 19 deletion")

ggplot(ROC_chr1_del,aes(x=`1-Specificity`,y=Sensitivity,col=sig_thres))+geom_line()+facet_wrap(~paste0("log2FC=",round(sig_thres,3))+paste0("AUC=",round(AUC,3),"\nture=",N_true," false=",N_false),nrow=2)+ggtitle("Chromosome 1 deletion")
dev.off()


