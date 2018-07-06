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
centromere=cytoband_red[chromArm=="p",c("chromosome","chromEnd"),with=FALSE]
setnames(centromere,"chromEnd","centromerePos")


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
WGS_CNV_dir=file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/cnvkit_vs_pga_hg38")
WGS_CNV_files=list.files(WGS_CNV_dir,pattern="call.cns",full.names=TRUE)

WGS_segments_all=data.table()
for (path in WGS_CNV_files){
  sample=gsub(".*/|_S[0-9]+.*","",path)
  CNV=fread(path,drop="gene")
  CNV[,sample:=sample]
  WGS_segments_all=rbindlist(list(WGS_segments_all,CNV))
}
WGS_segments_all[,chromosome:=gsub("chr","",chromosome),]

#annotate WGS CNV with cytoband (chromosome arms)
WGS_segments_all_gr=with(WGS_segments_all,GRanges(seqnames = Rle(chromosome), IRanges(start=start, end=end),strand=Rle("*"),sample=sample))
ov_CNV_seg=as.data.frame(findOverlaps(WGS_segments_all_gr,cytoband_red_gr,type="any",select="all"))
WGS_segments_cb=unique(cbind(WGS_segments_all[ov_CNV_seg$queryHits],cytoband_red[ov_CNV_seg$subjectHits,c("chromArm","chromStart","chromEnd","chromArm_length"),with=FALSE]))

WGS_segments_cb[,frag_len_WGS:=end-start,by=c("sample","chromosome","chromArm","start","end")]
WGS_segments_cb[,frag_len_corr_WGS:=min(end,chromEnd)-max(start,chromStart),by=c("sample","chromosome","chromArm","start","end")]
WGS_segments_cb[,variant:=ifelse(log2<0,"deletion",ifelse(log2>0,"amplification","nc")),]
WGS_segments_cb_red=WGS_segments_cb[,list(frag_len_corr_WGS=sum(frag_len_corr_WGS),log2_mean_WGS=weighted.mean(log2,frag_len_corr_WGS)),by=c("sample","chromosome","chromArm","chromArm_length","variant")]
WGS_segments_cb_red[,sample:=paste0(sample,"_fv")]


# combine RRBS and WGS results as segments to plot into one chromosome plot
WGS_segments_cb_uni=WGS_segments_cb[,c("sample","chromosome","chromArm","chromStart","chromEnd", "chromArm_length","start","end","log2","frag_len_WGS", "frag_len_corr_WGS"),with=FALSE]
setnames(WGS_segments_cb_uni,c("start","end","frag_len_WGS", "frag_len_corr_WGS","log2"),c("seg_start","seg_end","seg_len", "seg_len_corr","log2FC"))
WGS_segments_cb_uni[,sample:=paste0(sample,"_fv"),]
WGS_segments_cb_uni[,assay:="WGS",]


RRBS_segments_cb_uni=RRBS_segments_cb[,c("sample","chromosome","chromArm","chromStart","chromEnd", "chromArm_length","loc.start","loc.end","log2","frag_len_RRBS", "frag_len_corr_RRBS"),with=FALSE]
setnames(RRBS_segments_cb_uni,c("loc.start","loc.end","frag_len_RRBS", "frag_len_corr_RRBS","log2"),c("seg_start","seg_end","seg_len", "seg_len_corr","log2FC"))
RRBS_segments_cb_uni[,assay:="RRBS",]

segments_uni=rbindlist(list(WGS_segments_cb_uni,RRBS_segments_cb_uni),use.names=TRUE)
segments_uni=segments_uni[sample!="N1390_13_fv",,]
segments_uni=merge(segments_uni,centromere,by="chromosome")
segments_uni[,meanLog2FC:=weighted.mean(log2FC,seg_len),by=c("chromosome","sample","assay")]

segments_uni=merge(segments_uni,setnames(annotation[,c("patID","N_number_seq","surgery.x"),with=FALSE],"N_number_seq","sample"),by="sample")
segments_uni[is.na(surgery.x),surgery.x:=0,]
segments_uni[,sample_ID:=paste0(patID,"_",surgery.x),]

####Figure S4a
pdf("chrom_foldChange.pdf",height=6,width=6)
for (chr in c(as.character(1:22),"X","Y")){
  sub=segments_uni[chromosome==chr]
  sub[,sample_ID:=factor(sample_ID,levels=unique(sample_ID[assay=="WGS"][order(meanLog2FC[assay=="WGS"])])),]
  h_col="red";l_col="blue";m_col="white"

  range=max(sub$log2FC)-min(sub$log2FC)
  mini=min(sub$log2FC)
  lowValue=-1.0
  midValue=0.0
  higValue=1.0
  low=(lowValue-mini)/range
  mid=(midValue-mini)/range
  hig=(higValue-mini)/range
  print(c(low, mid, hig))
  if(low<=0 & hig>=1) {
    colors=c("green","gray96","gold")
    values=c(0,mid,1)
    print("01")
  } else if(hig>=1) {
    colors=c("blue","green","gray96","gold")
    values=c(0,low,mid,1)
    print("?1")
  } else if(low<=0) {
    colors=c("green","gray96","gold","red")
    values=c(0,mid,hig,1)
    print("0?")
  } else {
    colors=c("blue","green","gray96","gold","red")
    values=c(0,low,mid,hig,1)
    print("??")
  } 
  
#  pl1=ggplot(sub)+geom_segment(aes(y=sample,x=seg_start/1000000,yend=sample,xend=seg_end/1000000,col=log2FC,group=assay),size=3)+geom_point(aes(y=0,x=unique(centromerePos/1000000)),shape=4,size=3)+facet_wrap(~assay)+scale_color_gradient2(high=h_col,mid=m_col,low=l_col)+ylab("")+xlab("Chromosome position [Mb]")+ggtitle(paste0("Chromosome ", chr))
  
  pl2=ggplot(sub)+geom_segment(aes(y=sample_ID,x=seg_start/1000000,yend=sample_ID,xend=seg_end/1000000,col=log2FC,group=assay),size=3)+geom_vline(aes(xintercept=unique(centromerePos/1000000)),size=2)+facet_wrap(~assay)+scale_color_gradientn(colours=colors,values=values)+ylab("")+xlab("Chromosome position [Mb]")+ggtitle(paste0("Chromosome ", chr))


#  print(pl1)
  print(pl2)
}
dev.off()


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

# for selected variants and significance thresholds (focus on sig_0.5 because that is also what CNVkit selects)

#expand sothat all variants in all chromosome arms are represented in all samples with log2FC=0
exp=as.data.table(with(segments_cb_combined,expand.grid(sample=unique(sample),chromosome=unique(chromosome),chromArm=unique(chromArm),variant=unique(variant))))
segments_cb_combined_exp=merge(segments_cb_combined,exp,by=c("sample","chromArm","chromosome","variant"),all=TRUE)
segments_cb_combined_exp[is.na(log2_mean_WGS),log2_mean_WGS:=0]
segments_cb_combined_exp[is.na(log2_mean_RRBS),log2_mean_RRBS:=0]

segments_cb_combined_exp[,c("sig_0.1_WGS","sig_0.3_WGS","sig_0.5_WGS","sig_0.1_RRBS","sig_0.3_RRBS","sig_0.5_RRBS"):=list(abs(log2_mean_WGS)>=0.1,abs(log2_mean_WGS)>=0.3,abs(log2_mean_WGS)>=0.5,abs(log2_mean_RRBS)>=0.1,abs(log2_mean_RRBS)>=0.3,abs(log2_mean_RRBS)>=0.5),]

segments_cb_combined_exp[,c("sig_0.5_true","sig_0.5_false"):=list(sum(sig_0.5_WGS==TRUE),sum(sig_0.5_WGS==FALSE)),by=c("chromosome","variant")]

get_ROC=function(target,values){
  ROC=roc(target,values)
  dt=data.table(TPR=ROC$sensitivities,FPR=1-ROC$specificities,AUC=ROC$auc,run=0,rand=FALSE)
  dt=dt[order(FPR)]
  dt=dt[order(TPR)]
  
  set.seed(1234)
  rand_labs=lapply(seq_len(10), function(x) sample(target))
  i=1
  for (rand_lab in rand_labs){
    print(i)
    ROC=roc(rand_lab,values)
    rand_dt=data.table(TPR=ROC$sensitivities,FPR=1-ROC$specificities,AUC=ROC$auc,run=i,rand=TRUE)
    rand_dt=rand_dt[order(FPR)]
    rand_dt=rand_dt[order(TPR)]
    i=i+1
    dt=rbindlist(list(dt,rand_dt))
  }
  dt_int=dt[,approx(FPR,TPR,xout=seq(0,1,0.01),yleft=0,yright=1,method="constant",ties=max),by=c("AUC","run","rand")]
  setnames(dt_int,c("x","y"),c("FPR","TPR"))
  dt_int[FPR==0,TPR:=0,]
  dt_int[FPR==1,TPR:=1,]
  
  return(list(roc_dt=dt,roc_dt_int=dt_int))
}

ROC_sig_0.5=segments_cb_combined_exp[sig_0.5_true>0&sig_0.5_false>0,get_ROC(sig_0.5_WGS,log2_mean_RRBS)$roc_dt,by=c("chromosome","variant","sig_0.5_true","sig_0.5_false")]
ROC_sig_0.5_int=segments_cb_combined_exp[sig_0.5_true>0&sig_0.5_false>0,get_ROC(sig_0.5_WGS,log2_mean_RRBS)$roc_dt_int,by=c("chromosome","variant","sig_0.5_true","sig_0.5_false")]

ROC_sig_0.5[,facet_label:=paste0("Chr",chromosome," ",variant,"\nAUC=",unique(round(AUC[rand==FALSE],2)),"/",round(mean(AUC[rand==TRUE]),2),"\nture=",sig_0.5_true," false=",sig_0.5_false),by=c("chromosome","variant")]
ROC_sig_0.5_int[,facet_label:=paste0("Chr",chromosome," ",variant,"\nAUC=",unique(round(AUC[rand==FALSE],2)),"/",round(mean(AUC[rand==TRUE]),2),"\nture=",sig_0.5_true," false=",sig_0.5_false),by=c("chromosome","variant")]

####Figure S4b
pdf("chromArm_foldChange_ROC_sig_0.5_rand.pdf",height=10,width=27)
ggplot(ROC_sig_0.5,aes(x=FPR,y=TPR,col=rand,group=run,alpha=rand))+geom_line()+facet_wrap(~facet_label,nrow=4,scale="free")+scale_color_manual(values=c("TRUE"="grey","FALSE"="blue"))+scale_alpha_manual(values=c("TRUE"=0.5,"FALSE"=1))

ggplot(ROC_sig_0.5_int,aes(x=FPR,y=TPR,group=rand))+stat_summary(geom="ribbon", fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), fill="lightgrey",alpha=0.4)+stat_summary(geom="line",aes(col=rand), fun.y=mean)+facet_wrap(~facet_label,nrow=4,scale="free")+scale_color_manual(values=c("TRUE"="grey","FALSE"="blue"))
dev.off()

#plot individually
sub=segments_cb_combined[chromosome=="19"&variant=="deletion"]
plot(roc(sub$sig_0.1_WGS,sub$log2_mean_RRBS))
plot(roc(sub$sig_0.3_WGS,sub$log2_mean_RRBS))
plot(roc(sub$sig_0.5_WGS,sub$log2_mean_RRBS))

#more automated way (select range of significance levels based on actual lo2FC values)

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
    if (length(unique(abs(sub_exp$log2_mean_WGS)>=thres))<2){print("Skipping");next}
    
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


