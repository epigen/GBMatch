#NOTE: This script performs the MIRA analysis. (Inference of transcription factor/chromatin 
#protein activity through local DNA methylation depletion)
library(project.init)
project.init2("GBMatch")
source(file.path(getOption("PROJECT.DIR"),"src/99-MIRA.R"))
library(LOLA)

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/09-dipScore")
dir.create(out_dir)
setwd(out_dir)

eload(loadAnnotation())
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
transc_st=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.1-GBM_classifier/class_probs_annot_27_noNeural_predRRBS.tsv"))
transc_st=transc_st[!is.na(sub_group)]
transc_st[,sub_group_prob:=get(sub_group),by=sample]
setnames(transc_st, "sample","N_number_seq")
annotation=merge(annotation,transc_st[,-c("patID","IDH","surgery.x","Follow-up_years","WHO2016_classification","category"),],by="N_number_seq")

setkey(annotation, N_number_seq)
setkey(SV$psa,sample_name)
SV$psa[,cell_type:=NULL,]
SV$psa[annotation,cell_type:=ifelse(category=="control","white matter","glioblastoma"),allow=TRUE]
SV$psa[annotation,category:=category,allow=TRUE]

#get DNA methylation data
loadCaches("rrbsCg")

regionDB = loadRegionDB(file.path(Sys.getenv("RESOURCES"),"regions/LOLACore/hg38/"),collections=c("encode_tfbs","codex"))
regionDB_seg = loadRegionDB(file.path(Sys.getenv("RESOURCES"),"regions/customRegionDB/hg38/"),collections=list("roadmap_segmentation"))

cellType_anno=fread(file.path(getOption("PROJECT.DIR"),"metadata/LOLA_annot/CellTypes.tsv"))
regionDB$regionAnno[,GRL_id:=1:nrow(regionDB$regionAnno),]
regionDB$regionAnno=merge(regionDB$regionAnno,cellType_anno,by="cellType",all.x=TRUE)

cellType_anno_seg=fread(file.path(Sys.getenv("RESOURCES"),"regions/customRegionDB/hg38/roadmap_segmentation/index.txt"),select=c("filename","EID", "seg_code","seg_explanation","Epigenome name (from EDACC Release 9 directory)","ANATOMY"))
setnames(cellType_anno_seg,"Epigenome name (from EDACC Release 9 directory)","cellType_corr")
regionDB_seg$regionAnno[,GRL_id:=1:nrow(regionDB_seg$regionAnno),]
regionDB_seg$regionAnno=merge(regionDB_seg$regionAnno,cellType_anno_seg,by="filename",all.x=TRUE)

#select region sets to focus on (also used for follow-up analysis --> data export/last section)
select=c("CTCF__Gliobla_None_62750","CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465","CHD1_(A301-218A)__H1-hESC_None_7245","NRSF__H1-hESC_None_13281","NRSF__U87_None_11741","NANOG__Embryonic Stem Cell_NA_18532","SOX2__Embryonic Stem Cell_NA_5204","Pol2__Gliobla_None_17405","Pol2__H1-hESC_None_20311","EZH2_(39875)__NH-A_None_6488","EZH2__Embryonic Stem Cell_NA_3204","POU5F1__Embryonic Stem Cell_NA_10374","POU5F1_(SC-9081)__H1-hESC_None_3996","RBBP5_(A300-109A)__H1-hESC_None_16143","POLR2A__Embryonic Stem Cell_NA_14498","KDM4A__Embryonic Stem Cell_NA_20621")
select_seg=c("Active_TSS__H1_Cell_Line__E003_17209","Active_TSS__NH-A_Astrocytes__E125_20684","Bivalent_enhancers__H1_Cell_Line__E003_15759","Bivalent_enhancers__NH-A_Astrocytes__E125_5604","Bivalent_Poised_TSS__H1_Cell_Line__E003_9125","Bivalent_Poised_TSS__NH-A_Astrocytes__E125_3212","Enhancers__H1_Cell_Line__E003_98851","Enhancers__NH-A_Astrocytes__E125_100086","Heterochromatin__H1_Cell_Line__E003_35494","Heterochromatin__NH-A_Astrocytes__E125_19974","Quiescent__H1_Cell_Line__E003_82316","Quiescent__NH-A_Astrocytes__E125_63717","Repressed_polycomb__H1_Cell_Line__E003_17100","Repressed_polycomb__NH-A_Astrocytes__E125_19527")

total_width=5000
#analysis for TFs
#select single factors/cell types
factor="Astrocyte|ESC"
mylist=lapply(c(regionDB$regionAnno[size>1000][grep(factor,antibody)]$GRL_id,regionDB$regionAnno[size>1000][grep(factor,cellType_corr)]$GRL_id),function(x){name=paste0(regionDB$regionAnno[GRL_id==x]$antibody,"__",regionDB$regionAnno[GRL_id==x]$cellType,"_",regionDB$regionAnno[GRL_id==x]$treatment,"_",regionDB$regionAnno[GRL_id==x]$size);dt=as.data.table(as.data.frame(regionDB$regionGRL[[x]]));dt[,start:=ceiling(start-(total_width-width)/2),];dt[,end:=ceiling(end+(total_width-width)/2),];dt[,width:=end-start,];return(list("name"=name,"dt"=dt))})

rangesList=lapply(mylist,function(x){x$dt})
names(rangesList)=as.character(unlist(lapply(mylist,function(x){x$name})))
attach(rangesList)

binCount=21
binResults=binProcess(rangesList, list("rrbsCg"), annotation, binCount)

#analysis for segmentation
#select single factors/cell types
factor_seg="Astrocyte|H1_Cell_Line|H9_Cell_Line"
mylist_seg=lapply(c(regionDB_seg$regionAnno[size>1000][grep(factor_seg,antibody)]$GRL_id,regionDB_seg$regionAnno[size>1000][grep(factor_seg,cellType_corr)]$GRL_id),function(x){name=paste0(regionDB_seg$regionAnno[GRL_id==x]$seg_explanation,"__",regionDB_seg$regionAnno[GRL_id==x]$cellType_corr,"__",regionDB_seg$regionAnno[GRL_id==x]$EID,"_",regionDB_seg$regionAnno[GRL_id==x]$size);dt=as.data.table(as.data.frame(regionDB_seg$regionGRL[[x]]));dt[,start:=ceiling(start-(total_width-width)/2),];dt[,end:=ceiling(end+(total_width-width)/2),];dt[,width:=end-start,];return(list("name"=name,"dt"=dt))})

rangesList_seg=lapply(mylist_seg,function(x){x$dt})
names(rangesList_seg)=gsub("/| ","_",as.character(unlist(lapply(mylist_seg,function(x){x$name}))))
attach(rangesList_seg)

binCount=21
binResults_seg=binProcess(rangesList_seg, list("rrbsCg"), annotation, binCount)


#plot profiles by subgroup overview (all data)
pdf(paste0("aggregatedMeth_subgroup_",factor,".pdf"),height=6,width=8)
for(assay in names(binResults)){
  pl=ggplot(binResults[[assay]]$binnedBSDT[(readCount>1000&(auc>0.8&sub_group_prob>0.8)|category=="control")&category%in%c("control","GBMatch","GBmatch_val","GBmatch_valF")],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+facet_wrap(~category,scale="free")+ggtitle(assay)+scale_x_discrete(labels=labelBinCenter(binCount))+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))
  print(pl)
}
dev.off()

pdf(paste0("aggregatedMeth_subgroup_",factor_seg,"_seg.pdf"),height=6,width=8)
for(assay in names(binResults_seg)){
  pl=ggplot(binResults_seg[[assay]]$binnedBSDT[(readCount>1000&(auc>0.8&sub_group_prob>0.8)|category=="control")&category%in%c("control","GBMatch","GBmatch_val","GBmatch_valF")],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+facet_wrap(~category,scale="free")+ggtitle(assay)+scale_x_discrete(labels=labelBinCenter(binCount))+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))
  print(pl)
}
dev.off()


#plot profiles by subgroup only GBMatch IDH=="wt"
####Figure 2i
pdf(paste0("aggregatedMeth_subgroup_GBMatch_",factor,".pdf"),height=3.5,width=4)
for(assay in names(binResults)){
  spl=unlist(strsplit(gsub("_\\([0-9]*\\)","",assay),"_"))
  ab=spl[1]
  cellType=ifelse(grepl("embryonic|ESC",spl[3],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[3],ignore.case=TRUE),"Astrocytes",spl[3]))
  annot=paste0(ab," (",cellType,")")
  data=binResults[[assay]]$binnedBSDT[readCount>1000&auc>0.8&sub_group_prob>0.8&category=="GBMatch"&IDH=="wt",list(regionGroupID,methyl,sub_group,id)]
  data[,label:=annot,]
  data[,title:=assay,]
  pl=ggplot(data,aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
  if(assay %in% select){
    write.table(data,paste0("aggregatedMeth_subgroup_GBMatch_",assay,".csv"),sep=";",row.names=FALSE,quote=FALSE)
  }
  
}
dev.off()

#plot profiles by subgroup only GBMatch IDH=="wt" for segmentation
####Figure S7e
pdf(paste0("aggregatedMeth_subgroup_GBMatch_",factor_seg,"_seg.pdf"),height=3.5,width=4)
for(assay in names(binResults_seg)){
  spl=unlist(strsplit(assay,"__"))
  ab=spl[1]
  cellType=ifelse(grepl("embryonic|ESC",spl[2],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[2],ignore.case=TRUE),"Astrocytes",spl[2]))
  annot=paste0(ab,"\n(",cellType,")")
  pl=ggplot(binResults_seg[[assay]]$binnedBSDT[readCount>1000&auc>0.8&sub_group_prob>0.8&category=="GBMatch"&IDH=="wt"],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
}
dev.off()

#plot profiles by surgery only GBMatch IDH=="wt"
pdf(paste0("aggregatedMeth_surgery_GBMatch_",factor,".pdf"),height=3.5,width=4)
for(assay in names(binResults)){
  spl=unlist(strsplit(gsub("_\\([0-9]*\\)","",assay),"_"))
  ab=spl[1]
  cellType=ifelse(grepl("embryonic|ESC",spl[3],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[3],ignore.case=TRUE),"Astrocytes",spl[3]))
  annot=paste0(ab," (",cellType,")")
  pl=ggplot(binResults[[assay]]$binnedBSDT[readCount>1000&surgery.x%in%c(1,2)&category=="GBMatch"&IDH=="wt"],aes(x=factor(regionGroupID),y=methyl,col=as.factor(surgery.x)))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=as.factor(surgery.x)), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
}
dev.off()

#plot profiles by surgery only GBMatch IDH=="wt" for segmentation
pdf(paste0("aggregatedMeth_surgery_GBMatch_",factor_seg,"_seg.pdf"),height=3.5,width=4)
for(assay in names(binResults_seg)){
  spl=unlist(strsplit(assay,"__"))
  ab=spl[1]
  cellType=ifelse(grepl("H1_|H9_",spl[2],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[2],ignore.case=TRUE),"Astrocytes",spl[3]))
  annot=paste0(ab,"\n(",cellType,")")
  pl=ggplot(binResults_seg[[assay]]$binnedBSDT[readCount>1000&surgery.x%in%c(1,2)&category=="GBMatch"&IDH=="wt"],aes(x=factor(regionGroupID),y=methyl,col=as.factor(surgery.x)))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=as.factor(surgery.x)), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
}
dev.off()


#plot profiles by subgroup only GBMatch_val IDH=="wt"
pdf(paste0("aggregatedMeth_subgroup_GBMatch_val_",factor,".pdf"),height=3.5,width=4)
for(assay in names(binResults)){
  spl=unlist(strsplit(gsub("_\\([0-9]*\\)","",assay),"_"))
  ab=spl[1]
  cellType=ifelse(grepl("embryonic|ESC",spl[3],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[3],ignore.case=TRUE),"Astrocytes",spl[3]))
  annot=paste0(ab," (",cellType,")")
  pl=ggplot(binResults[[assay]]$binnedBSDT[readCount>1000&auc>0.8&sub_group_prob>0.8&category=="GBmatch_val"&IDH=="wt"],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
}
dev.off()

#plot profiles by subgroup only GBMatch_val IDH=="wt" segmentation
pdf(paste0("aggregatedMeth_subgroup_GBMatch_val_",factor_seg,"_seg.pdf"),height=3.5,width=4)
for(assay in names(binResults_seg)){
  spl=unlist(strsplit(assay,"__"))
  ab=spl[1]
  cellType=ifelse(grepl("embryonic|ESC",spl[3],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[3],ignore.case=TRUE),"Astrocytes",spl[3]))
  annot=paste0(ab,"\n(",cellType,")")
  pl=ggplot(binResults_seg[[assay]]$binnedBSDT[readCount>1000&auc>0.8&sub_group_prob>0.8&category=="GBmatch_val"&IDH=="wt"],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
}
dev.off()


#collect dip scores for TFs
combined_dipScores=binResults[[1]]$dipScores[,c("id","surgery","sub_group","sub_group_prob","category","IDH"),with=FALSE]

for(assay in names(binResults)){
  print(assay)
  dipScores=binResults[[assay]]$dipScores[,c("id","dipScore"),with=FALSE]
  setnames(dipScores,"dipScore",assay)
  if (length(combined_dipScores)==0){combined_dipScores=dipScores}else{
    combined_dipScores=merge(combined_dipScores,dipScores,by="id")
  }
}

sort(names(combined_dipScores))
write.table(combined_dipScores[,c("id",select),with=FALSE],paste0("dipScores_",factor,"_sel.tsv"),sep="\t",row.names=FALSE,quote=FALSE)


#collect dip scores for segmentation
combined_dipScores_seg=binResults_seg[[1]]$dipScores[,c("id","surgery","sub_group","sub_group_prob","category","IDH"),with=FALSE]

for(assay in names(binResults_seg)){
  print(assay)
  dipScores_seg=binResults_seg[[assay]]$dipScores[,c("id","dipScore"),with=FALSE]
  setnames(dipScores_seg,"dipScore",assay)
  if (length(combined_dipScores_seg)==0){combined_dipScores_seg=dipScores_seg}else{
    combined_dipScores_seg=merge(combined_dipScores_seg,dipScores_seg,by="id")
  }
}

sort(names(combined_dipScores_seg))
write.table(combined_dipScores_seg[,c("id",select_seg),with=FALSE],paste0("dipScores_",factor_seg,"_sel_seg.tsv"),sep="\t",row.names=FALSE,quote=FALSE)




