#NOTE: This script checks for associations between various annotation covariates.
#Mostly categorical vs continuous using the wilcoxon rank sum test.
library(project.init)
project.init2("GBMatch")
library(gtools)


out = '12.1-Annot_Associations_follow-up/'
dir.create(dirout(out))
setwd(dirout(out))

##### PREPARE ANNOTATION DATA

annotation = fread(dirout("01.1-combined_annotation/", "annotation_combined_final.tsv"))
load(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/column_annotation_combined.RData"))

annotation[,long_survivor:=ifelse(`Follow-up_years`>=3,TRUE,FALSE),]
annotation[,StuppComplete:=ifelse(StuppComplete=="yes",TRUE,FALSE),]
annotation[,VPM:=all_count/bg_calls*1000000,]
#check progression_type overlap --> some patients have more than one type --> resolve
table(annotation[,psum(`classic T1`,`cT1 flare up`,`T2 diffus`,`T2 circumscribed`),])
#merge progression types into 1 column (T2 diffuse dominates cT1 flare up)
annotation[,progression_types:=ifelse(`T2 diffus`==1,"T2 diffus",ifelse(`classic T1`==1,"classic T1",ifelse(`cT1 flare up`==1,"cT1 flare up",ifelse(`T2 circumscribed`==1,"T2 circumscribed",NA)))),]

annotation[,ShiftPhenotype:=TumorPhenotype,]
annotation[Shape_shift=="stable",ShiftPhenotype:="stable",]

##### Functions
test_difference=function(values,groups,ntests){
  if(length(na.omit(unique(groups)))<2){
    combinations=as.data.table(permutations(n=2,r=2,v=c("d1","d2"),repeats.allowed=FALSE))
  }else{
    combinations=as.data.table(permutations(n=length(na.omit(unique(groups))),r=2,v=as.character(unique(groups)),repeats.allowed=FALSE))
  }
  combinations[,id:=paste0(sort(c(substr(V1,1,3),substr(V2,1,3))),collapse=","),1:nrow(combinations)]
  combinations=combinations[!duplicated(id)]
  
  if (length(na.omit(values))>5){
    combinations[,p.value:=ifelse(sum(!is.na(values[groups==V2]))>0&sum(!is.na(values[groups==V1]))>0,try(signif(as.numeric(wilcox.test(x=values[groups==V1],y=values[groups==V2])$p.value),3),silent=TRUE),1),by=1:nrow(combinations)]}else{
      combinations[,p.value:=1]
    }
  combinations[,p.value.round:=ifelse(p.value<0.001,"<0.001",paste0("=",as.character(signif(p.value,2)))),]
  combinations[,annot:=paste0("p-val(",id,")=",p.value),]
  combinations[,N_V1:=length(na.omit(values[groups==V1])),by=1:nrow(combinations)]
  combinations[,N_V2:=length(na.omit(values[groups==V2])),by=1:nrow(combinations)]
  combinations[,meanV1:=signif(mean(values[groups==V1],na.rm=TRUE),3),by=1:nrow(combinations)]
  combinations[,meanV2:=signif(mean(values[groups==V2],na.rm=TRUE),3),by=1:nrow(combinations)]
  combinations[,sdV1:=signif(sd(values[groups==V1],na.rm=TRUE),3),by=1:nrow(combinations)]
  combinations[,sdV2:=signif(sd(values[groups==V2],na.rm=TRUE),3),by=1:nrow(combinations)]
  return(combinations)
}

plot_boxplots=function(sub,test,test_category,transc_stats,class,by){
  
  message(test)
  ties=table(sub[,get(test),])[table(sub[,get(test),])>1]
  if(length(ties>1)){print(ties)}
  
  sig=sub[,test_difference(get(test),get(class),ntests=length(test_list)),by=by]
  sig[,test:=test,]
  sig[,test_category:=test_category,]
  transc_stats=rbindlist(list(transc_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.2*(ypos-min(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE))
  sig_annot=sig[,list(annot=paste0(annot,collapse="\n")),by=by]
  sig_annot[,ypos:=ypos]  
  sig_annot[,xpos:=0.5]
  
  pl=ggplot()+geom_point(data=sub,aes(x=substr(get(class),1,3),y=get(test),group=get(class)),shape=21,color="grey",alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_annot,lineheight=0.9,aes(x=xpos,y=ypos,label=annot),hjust=0,vjust=1,size=2.5)+geom_boxplot(data=sub,aes(x=substr(get(class),1,3),y=get(test),group=get(class)),outlier.size=NA,fill="transparent")+ylab(test)+xlab("")+facet_wrap(as.formula(paste0("~",paste0(by,collapse="+"))))
  print(pl)
  source_data=sub
  setnames(source_data,c(class,test),c("varX","varY"))
  source_data[,varY_name:=test,]
  source_data[,varX_name:=class,]
  return(list(transc_stats[N_V1>0&N_V2>0],source_data))
}

summarize_stats=function(stats,test_category,height=6,width=6,by){
  transc_stats[,test:=factor(test,levels=unique(test[order(get(by))])),]
  transc_stats[,p.value.top:=ifelse(p.value<=0.05,ifelse(p.value<=0.01,ifelse(p.value<=0.001,"<=0.001","<=0.01"),"<=0.05"),">0.05"),]
  
  pdf(paste0(test_category,"_statsSummary.pdf"),height=height,width=width)
  print(ggplot(transc_stats)+geom_tile(aes(x=id,y=test,fill=p.value.top),col="grey")+geom_point(aes(x=0,y=test,col=test_category),size=4,shape=15)+scale_fill_manual(values=c("<=0.001"="#f03b20","<=0.01"="#feb24c","<=0.05"="#ffeda0",">0.05"="white"))+ylab("")+xlab("")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(as.formula(paste0("~",paste0(by,collapse="+")))))
  dev.off()
}


##############dip score###########

test_list=c("CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465" ,"CHD1_(A301-218A)__H1-hESC_None_7245","NANOG__Embryonic Stem Cell_NA_18532"  ,"SOX2__Embryonic Stem Cell_NA_5204", "Pol2__H1-hESC_None_20311", "EZH2_(39875)__NH-A_None_6488" ,"EZH2__Embryonic Stem Cell_NA_3204","POU5F1__Embryonic Stem Cell_NA_10374","RBBP5_(A300-109A)__H1-hESC_None_16143","POLR2A__Embryonic Stem Cell_NA_14498","KDM4A__Embryonic Stem Cell_NA_20621","Active_TSS__H1_Cell_Line__E003_17209","Active_TSS__NH-A_Astrocytes__E125_20684","Bivalent_enhancers__H1_Cell_Line__E003_15759","Bivalent_enhancers__NH-A_Astrocytes__E125_5604","Bivalent_Poised_TSS__H1_Cell_Line__E003_9125","Bivalent_Poised_TSS__NH-A_Astrocytes__E125_3212","Enhancers__H1_Cell_Line__E003_98851","Enhancers__NH-A_Astrocytes__E125_100086","Heterochromatin__H1_Cell_Line__E003_35494","Heterochromatin__NH-A_Astrocytes__E125_19974","Quiescent__H1_Cell_Line__E003_82316","Quiescent__NH-A_Astrocytes__E125_63717","Repressed_polycomb__H1_Cell_Line__E003_17100","Repressed_polycomb__NH-A_Astrocytes__E125_19527")

class_prob=0.8  
classes=list(c("sub_group",3,4),c("surgery",3,4),c("progression_types",4,4),c("ShiftPhenotype",3,4),c("category",3,3.5))

for(class in classes){
  message(class)
  transc_stats=data.table()
  source_data=data.table()
  pdf(paste0(class[1],"_dip-score.pdf"),height=as.numeric(class[2]),width=as.numeric(class[3]))
  for (test in test_list){
    if (class[1]=="sub_group"){
      ####Figure 2i; S7a; S7f
      sub=unique(annotation[sub_group_prob>=class_prob&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="category"
    }
    if (class[1]=="surgery"){
      sub=unique(annotation[surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="progression_types"){
      sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="ShiftPhenotype"){
      sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="category"){
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery"),with=FALSE])
      by="surgery"
    }
    
    res=plot_boxplots(sub,test,"dip-score",transc_stats,class=class[1],by=by)
    transc_stats=res[[1]]
    source_data=rbindlist(list(source_data,res[[2]]))
  }
  dev.off()
  transc_stats[,p.value.BH:=signif(p.adjust(p.value,"BH"),3),by=c(by,"id")]
  transc_stats[,p.value.BF:=signif(p.adjust(p.value,"bonferroni"),3),by=c(by,"id")]
  write.table(transc_stats,paste0(class[1],"_dip-score_stats.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(source_data,paste0(class[1],"_dip-score_stats_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  summarize_stats(transc_stats,paste0(class[1],"_dip-score"),by=by)
}

#correlation between class probability and MIRA score
####Figure S7b
auc_thres=0.8
classes=c("Classical","Mesenchymal","Proneural")
pdf("transcSG-dip_cor.pdf",height=6,width=7)
for (test in test_list){
  sub=unique(annotation[auc>=auc_thres&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_seq",test,classes,"surgery","category"),with=FALSE])
  sub_long=melt(sub,id.vars=c("N_number_seq",test,"category"),measure.vars=c("Classical","Mesenchymal","Proneural"))
  my_cors=sub_long[,list(xpos=min(value,na.rm=TRUE)+0.1,ypos=max(get(test),na.rm=TRUE)+0.1,cor=cor(get(test),value,use="complete.obs"),cor.test=cor.test(get(test),value,use="complete.obs")$p.value,N=nrow(na.omit(cbind(get(test),value)))),by=c("variable","category")]
  
  pl=ggplot(sub_long, aes(x=value,y=get(test)))+geom_point(shape=21,color="grey",alpha=0.7,size=2)+geom_text(data=my_cors,hjust=0,vjust=1,aes(y=ypos,x=xpos,label=paste0("r= ",round(cor,3),"\np=",signif(cor.test,3),"\nN= ",N)))+geom_smooth(method="lm",aes(fill=variable,col=variable))+facet_wrap(category~variable,scale="free")+ggtitle(test)+xlim(c(0,1))+xlab("Class probability")+ylab("MIRA score")+theme(aspect.ratio=1)
  print(pl)
}
dev.off()


############## histo_immuno #################################
test_list=c("CD163","CD3","CD34","CD45ro","CD68","CD8","cell","EZH2","FOXP3","HLA-DR","MIB")

classes=list(c("sub_group",3,3.5),c("surgery",3,3.5),c("progression_types",3,3),c("ShiftPhenotype",3,3),c("long_survivor",3,3.5),c("category",3,3.5))
for(class in classes){
  message(class)
  transc_stats=data.table()
  source_data=data.table()
  pdf(paste0(class[1],"_histo_immune.pdf"),height=as.numeric(class[2]),width=as.numeric(class[3]))
  for (test in test_list){
    if (class[1]=="sub_group"){
      ####Figure 3a; 4a; S7d; S8a; S8c
      sub=unique(annotation[sub_group_prob>=class_prob&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="category"
    }
    if (class[1]=="surgery"){
      ####Figure 3d; S7d; S8c; S11a
      sub=unique(annotation[surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="progression_types"){
      ####Figure 3f
      sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="ShiftPhenotype"){
      ####Figure 4h
      sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="long_survivor"){
      sub=unique(annotation[!is.na(long_survivor)&surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="category"){
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery"),with=FALSE])
      by="surgery"
    }
    res=plot_boxplots(sub,test,"histo_immuno",transc_stats,class=class[1],by=by)
    transc_stats=res[[1]]
    source_data=rbindlist(list(source_data,res[[2]]))
  }
  dev.off()
  transc_stats[,p.value.BH:=signif(p.adjust(p.value,"BH"),3),by=c(by,"id")]
  transc_stats[,p.value.BF:=signif(p.adjust(p.value,"bonferroni"),3),by=c(by,"id")]
  write.table(transc_stats,paste0(class[1],"_histo_immuno_stats.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(source_data,paste0(class[1],"_histo_immuno_stats_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  summarize_stats(transc_stats,paste0(class[1],"_histo_immuno"),by=by)
}

#correlateion between EZH2+ cells and EZH2 dip score
sub=annotation[category%in%c("GBMatch")&IDH=="wt",c("N_number_st","EZH2","EZH2_(39875)__NH-A_None_6488","sub_group","Mesenchymal","sub_group_prob","surgery"),with=FALSE]
cors=unique(sub[,cor.test(EZH2,`EZH2_(39875)__NH-A_None_6488`,use="complete.obs"),by="sub_group"][,c("sub_group","p.value","estimate"),with=FALSE])
Ns=unique(sub[,list(N=nrow(na.omit(cbind(EZH2,`EZH2_(39875)__NH-A_None_6488`)))),by="sub_group"])
cors=merge(cors,Ns,by="sub_group")

####Figure S7c
pdf("EZH2_validation.pdf",height=2.5,width=4.5)
ggplot(sub,aes(x=EZH2,y=`EZH2_(39875)__NH-A_None_6488`,fill=sub_group,col=sub_group))+geom_point(shape=21,alpha=0.7)+geom_smooth(method="lm")+geom_text(data=cors,x=0.4,y=c(0.7,0.63,0.56),aes(label=paste0("r=",round(estimate,2)," p=",round(p.value,2)," N=",N)))+ylab("MIRA score EZH2 (Astrocytes)")+xlab("Fraction EZH2+ cells")
dev.off()

##################histo_segmentation#######################################
#comment: ties warning and equel values for different samples is ok.
test_list=grep("Mean |Mode |File|Image|Block| STD |mm|background|Area|Necrosis|Whole|Scars",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)
classes=list(c("sub_group",3,3.5),c("surgery",3,3.5),c("progression_types",3,3),c("ShiftPhenotype",3,3),c("long_survivor",3,3.5),c("category",3,3.5))
for(class in classes){
  message(class)
  transc_stats=data.table()
  source_data=data.table()
  pdf(paste0(class[1],"_histo_segmentation.pdf"),height=as.numeric(class[2]),width=as.numeric(class[3]))
  for (test in test_list){
    if (class[1]=="sub_group"){
      ####Figure S8d
      sub=unique(annotation[sub_group_prob>=class_prob&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="category"
    }
    if (class[1]=="surgery"){
      ####Figure S8e; S11f
      sub=unique(annotation[surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="progression_types"){
      sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="ShiftPhenotype"){
      ####Figure 4f,h
      sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="long_survivor"){
      sub=unique(annotation[!is.na(long_survivor)&surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="category"){
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery"),with=FALSE])
      by="surgery"
    }
    res=plot_boxplots(sub,test,"histo_segmentation",transc_stats,class=class[1],by=by)
    transc_stats=res[[1]]
    source_data=rbindlist(list(source_data,res[[2]]))
  }
  dev.off()
  transc_stats[,p.value.BH:=signif(p.adjust(p.value,"BH"),3),by=c(by,"id")]
  transc_stats[,p.value.BF:=signif(p.adjust(p.value,"bonferroni"),3),by=c(by,"id")]
  write.table(transc_stats,paste0(class[1],"_histo_segmentation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(source_data,paste0(class[1],"_histo_segmentation_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  
  summarize_stats(transc_stats,paste0(class[1],"_histo_segmentation"),by=by)
}

################imaging_segmentation###################################
#comment: ties warning and equel values for different samples is ok.
test_list=grep("ID|surgery|VASARI",column_annotation_combined_clean$imaging_segmentation,invert=TRUE,value=TRUE)
classes=list(c("sub_group",3,3.5),c("surgery",3,3.5),c("progression_types",3,3),c("ShiftPhenotype",3,3),c("long_survivor",3,3.5),c("category",3,3.5))
for(class in classes){
  message(class)
  transc_stats=data.table()
  source_data=data.table()
  pdf(paste0(class[1],"_imaging_segmentation.pdf"),height=as.numeric(class[2]),width=as.numeric(class[3]))
  for (test in test_list){
    if (class[1]=="sub_group"){
      ####Figure S9b,c
      sub=unique(annotation[sub_group_prob>=class_prob&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="category"
    }
    if (class[1]=="surgery"){
      ####Figure S9b,d
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="progression_types"){
      sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="ShiftPhenotype"){
      sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="long_survivor"){
      sub=unique(annotation[!is.na(long_survivor)&surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="category"){
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c("N_number_st",test,class[1],"surgery"),with=FALSE])
      by="surgery"
    }
    res=plot_boxplots(sub,test,"imaging_segmentation",transc_stats,class=class[1],by=by)
    transc_stats=res[[1]]
    source_data=rbindlist(list(source_data,res[[2]]))
  }
  dev.off()
  transc_stats[,p.value.BH:=p.adjust(p.value,"BH"),by=c(by,"id")]
  transc_stats[,p.value.BF:=p.adjust(p.value,"bonferroni"),by=c(by,"id")]
  write.table(transc_stats,paste0(class[1],"_imaging_segmentation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(transc_stats,paste0(class[1],"_imaging_segmentation_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  summarize_stats(transc_stats,paste0(class[1],"_imaging_segmentation"),by=by)
}


###########molecular data##################################
test_list=c(column_annotation_combined_clean$meth_heterogeneity,"VPM")

classes=list(c("sub_group",3,3.5),c("surgery",3,3.5),c("progression_types",3,3),c("ShiftPhenotype",3,3),c("long_survivor",3,3.5),c("StuppComplete",3,4.5),c("category",3,3.5))
for(class in classes){
  message(class)
  transc_stats=data.table()
  source_data=data.table()
  pdf(paste0(class[1],"_meth_het.pdf"),height=as.numeric(class[2]),width=as.numeric(class[3]))
  for (test in test_list){
    if (class[1]=="sub_group"){
      sub=unique(annotation[sub_group_prob>=class_prob&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="category"
    }
    if (class[1]=="surgery"){
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }
    if (class[1]=="progression_types"){
      sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="ShiftPhenotype"){
      sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by="surgery"
    }
    if (class[1]=="long_survivor"){
      sub=unique(annotation[!is.na(long_survivor)&surgery%in%c(1,2,3)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"category"),with=FALSE])
      by="category"  
    }   
    if (class[1]=="StuppComplete"){
      sub=unique(annotation[!is.na(StuppComplete)&surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"surgery","category"),with=FALSE])
      by=c("surgery","category")
    }
    if (class[1]=="category"){
      sub=unique(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&enrichmentCycles>12&enrichmentCycles<16,c("N_number_st",test,class[1],"surgery"),with=FALSE])
      by="surgery"
    }
    res=plot_boxplots(sub,test,"meth_het",transc_stats,class=class[1],by=by)
    transc_stats=res[[1]]
    source_data=rbindlist(list(source_data,res[[2]]))
  }
  dev.off()
  transc_stats[,p.value.BH:=p.adjust(p.value,"BH"),by=c(by,"id")]
  transc_stats[,p.value.BF:=p.adjust(p.value,"bonferroni"),by=c(by,"id")]
  write.table(transc_stats,paste0(class[1],"_meth_het.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(transc_stats,paste0(class[1],"_meth_het_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  summarize_stats(transc_stats,paste0(class[1],"_meth_het"),by=by)
}

###########Meth heterogeneity###################################
annotation[,immune_cells:=sum(c(CD163,CD3,CD34,CD45ro,CD68,CD8,CD80,FOXP3),na.rm=TRUE),by=1:nrow(annotation)]
sub=annotation[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery%in%c(1,2)&enrichmentCycles>12&enrichmentCycles<16]

factors=c("Area total [mm²]","Area tumor [mm²]","Total (mm3)","Enhancing (mm3)","Relative share tumor","immune_cells")

####Figure S12b
pdf(paste0("heterogeneity_size.pdf"),height=2.5,width=5.5)
for (factor in factors){
  cor=sub[,list(p.value=cor.test(get(factor),mean_entropy,use="pairwise.complete.obs")$p.value,cor=cor.test(get(factor),mean_entropy,use="pairwise.complete.obs")$estimate,N=nrow(na.omit(cbind(get(factor),mean_entropy)))),by=c("surgery","category")]
  x=sub[,max(get(factor),na.rm=TRUE)/3,]
  ggpl=ggplot(sub,aes(x=get(factor),y=mean_entropy))+geom_point(shape=21,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cor,aes(x=x,y=45,label=paste0("r=",signif(cor,2)," p=",signif(p.value,2)," N=",N)),size=2.5)+facet_wrap(~surgery+category,scale="free")+xlab(factor)
  print(ggpl)
  
  cor=sub[,list(p.value=cor.test(get(factor),mean_pdr,use="pairwise.complete.obs")$p.value,cor=cor.test(get(factor),mean_pdr,use="pairwise.complete.obs")$estimate,N=nrow(na.omit(cbind(get(factor),mean_entropy)))),by=c("surgery","category")]
  x=sub[,max(get(factor),na.rm=TRUE)/3,]
  ggpl=ggplot(sub,aes(x=get(factor),y=mean_pdr*100))+geom_point(shape=21,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cor,aes(x=x,y=30,label=paste0("r=",signif(cor,2)," p=",signif(p.value,2)," N=",N)),size=2.5)+facet_wrap(~surgery+category,scale="free")+xlab(factor)
  print(ggpl)
}
dev.off()

