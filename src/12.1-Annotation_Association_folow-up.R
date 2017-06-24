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
#check progression_type overlap --> some patients have more than one type --> resolve
table(annotation[,psum(`classic T1`,`cT1 flare up`,`T2 diffus`,`T2 circumscribed`),])
#merge progression types into 1 column (T2 diffuse dominates cT1 flare up)
annotation[,progression_types:=ifelse(`T2 diffus`==1,"T2 diffus",ifelse(`classic T1`==1,"classic T1",ifelse(`cT1 flare up`==1,"cT1 flare up",ifelse(`T2 circumscribed`==1,"T2 circumscribed",NA)))),]

annotation[,ShiftPhenotype:=TumorPhenotype,]
annotation[Shape_shift=="stable",ShiftPhenotype:="stable",]

##### Functions

test_difference=function(values,groups,ntests){
  combinations=as.data.table(permutations(n=length(na.omit(unique(groups))),r=2,v=unique(groups),repeats.allowed=FALSE))
  combinations[,id:=paste0(sort(c(substr(V1,1,3),substr(V2,1,3))),collapse="_"),1:nrow(combinations)]
  combinations=combinations[!duplicated(id)]
  
  combinations[,p.value:=as.numeric(p.adjust(wilcox.test(x=values[groups==V1],y=values[groups==V2])$p.value),n=ntests,method="BH"),by=1:nrow(combinations)]
  combinations[,p.value.round:=ifelse(p.value<0.001,"<0.001",as.character(round(p.value,3))),]
  combinations[,annot:=paste0(id,": ",p.value.round),]
  combinations[,meanV1:=mean(values[groups==V1],na.rm=TRUE),by=1:nrow(combinations)]
  combinations[,meanV2:=mean(values[groups==V2],na.rm=TRUE),by=1:nrow(combinations)]
  combinations[,sdV1:=sd(values[groups==V1],na.rm=TRUE),by=1:nrow(combinations)]
  combinations[,sdV2:=sd(values[groups==V2],na.rm=TRUE),by=1:nrow(combinations)]
  return(combinations)
}


##############transc subgroups###########
transc_stats=data.table()

#transc sugroups vs dip scores
test_list=c("CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465" ,"CHD1_(A301-218A)__H1-hESC_None_7245","NANOG__Embryonic Stem Cell_NA_18532"  ,"SOX2__Embryonic Stem Cell_NA_5204", "Pol2__H1-hESC_None_20311", "EZH2_(39875)__NH-A_None_6488" ,"EZH2__Embryonic Stem Cell_NA_3204","POU5F1__Embryonic Stem Cell_NA_10374","RBBP5_(A300-109A)__H1-hESC_None_16143","POLR2A__Embryonic Stem Cell_NA_14498","KDM4A__Embryonic Stem Cell_NA_20621")

class_prob=0.8
classes="sub_group"

pdf("transcSG-dip.pdf",height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[sub_group_prob>class_prob&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntests=length(test_list)))
  sig[,test:=test,]
  sig[,category:="dip-score",]
  transc_stats=rbindlist(list(transc_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)+0.2
  pl=ggplot(sub, aes(x=substr(get(classes),1,3),y=get(test),group=get(classes)))+geom_point(shape=21,color="grey",alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


auc_thres=0.8
classes=c("Classical","Mesenchymal","Proneural")
pdf("transcSG-dip_cor.pdf",height=3,width=7)
for (test in test_list){
  sub=unique(annotation[auc>auc_thres&category=="GBMatch"&IDH=="wt",c("N_number_seq",test,classes,"surgery"),with=FALSE])
  sub_long=melt(sub,id.vars=c("N_number_seq",test),measure.vars=c("Classical","Mesenchymal","Proneural"))
  my_cors=sub_long[,list(xpos=min(value,na.rm=TRUE)+0.1,ypos=max(get(test),na.rm=TRUE)+0.1,cor=cor(get(test),value,use="complete.obs"),cor.test=cor.test(get(test),value,use="complete.obs")$p.value),by=variable]
  
  pl=ggplot(sub_long, aes(x=value,y=get(test)))+geom_point(shape=21,color="grey",alpha=0.7,size=2)+geom_text(data=my_cors,hjust=0,vjust=1,aes(y=ypos,x=xpos,label=paste0("r= ",round(cor,3),"\np= 10^",round(log10(cor.test),0))))+geom_smooth(method="lm",aes(fill=variable,col=variable))+facet_wrap(~variable,scale="free")+ggtitle(test)+xlim(c(0,1))+xlab("Class probability")+ylab("MIRA score")+theme(aspect.ratio=1)
  print(pl)
}
dev.off()



#subgroup vs histo_immuno
test_list=column_annotation_combined_clean$histo_immuno
classes="sub_group"

pdf(paste0("transcSG-histo_immune.pdf"),height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[sub_group_prob>class_prob&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="histo_immuno",]
  transc_stats=rbindlist(list(transc_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=substr(get(classes),1,3),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#subgroup vs histo_segmentation
test_list=grep("Mean |Mode |File|Image|Block| STD |mm",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)
classes="sub_group"

pdf(paste0("transcSG-histo_segmentation.pdf"),height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[sub_group_prob>class_prob&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="histo_segmentation",]
  transc_stats=rbindlist(list(transc_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=substr(get(classes),1,3),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#subgroup vs imaging_segmentation
test_list=grep("ID|surgery",column_annotation_combined_clean$imaging_segmentation,invert=TRUE,value=TRUE)
classes="sub_group"

pdf(paste0("transcSG-imaging_segmentation.pdf"),height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[sub_group_prob>class_prob&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="imaging_segmentation",]
  transc_stats=rbindlist(list(transc_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=substr(get(classes),1,3),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()

write.table(transc_stats,"transc_stats.tsv",row.names=FALSE,quote=FALSE,sep="\t")

#summarize transc_stats
transc_stats[,test:=factor(test,levels=unique(test[order(category)])),]
transc_stats[,p.value.top:=ifelse(p.value<=0.05,ifelse(p.value<=0.01,ifelse(p.value<=0.001,"<=0.001","<=0.01"),"<=0.05"),">0.05"),]

pdf(paste0("transc_stats.pdf"),height=9,width=6)
ggplot(transc_stats)+geom_tile(aes(x=id,y=test,fill=p.value.top),col="grey")+geom_point(aes(x=0,y=test,col=category),size=4,shape=15)+scale_fill_manual(values=c("<=0.001"="#f03b20","<=0.01"="#feb24c","<=0.05"="#ffeda0",">0.05"="white"))+ylab("")+xlab("")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



###########surgery###################################
surgery_stats=data.table()

#surgery vs dip scores
test_list=c("CTCF__Gliobla_None_62750","CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465" ,"CHD1_(A301-218A)__H1-hESC_None_7245"  ,"NRSF__H1-hESC_None_13281" ,"NRSF__U87_None_11741" ,"NANOG__Embryonic Stem Cell_NA_18532"  ,"SOX2__Embryonic Stem Cell_NA_5204","Pol2__Gliobla_None_17405", "Pol2__H1-hESC_None_20311", "EZH2_(39875)__NH-A_None_6488" ,"EZH2__Embryonic Stem Cell_NA_3204")

classes="surgery"

pdf("surgery-dip.pdf",height=3,width=3)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2,3)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntests=length(test_list)))
  sig[,test:=test,]
  sig[,category:="dip-score",]
  surgery_stats=rbindlist(list(surgery_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)+0.2
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(shape=21,color="grey",alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()



#surgery vs histo_immuno
test_list=column_annotation_combined_clean$histo_immuno
classes="surgery"

pdf(paste0("surgery-histo_immune.pdf"),height=3,width=2)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2,3)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="histo_immuno",]
  surgery_stats=rbindlist(list(surgery_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#surgery vs histo_segmentation
test_list=grep("Mean |Mode |File|Image|Block| STD |mm",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)
classes="surgery"

pdf(paste0("surgery-histo_segmentation.pdf"),height=3,width=2)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2,3)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="histo_segmentation",]
  surgery_stats=rbindlist(list(surgery_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#surgery vs imaging_segmentation
test_list=grep("ID|surgery",column_annotation_combined_clean$imaging_segmentation,invert=TRUE,value=TRUE)
classes="surgery"

pdf(paste0("surgery-imaging_segmentation.pdf"),height=3,width=2)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="imaging_segmentation",]
  surgery_stats=rbindlist(list(surgery_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=as.factor(get(classes)),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#surgery vs molecular data
test_list=column_annotation_combined_clean$meth_heterogeneity
classes="surgery"

pdf(paste0("surgery-meth_het.pdf"),height=3,width=2)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  sig[,test:=test,]
  sig[,category:="meth_het",]
  surgery_stats=rbindlist(list(surgery_stats,sig))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.1*ypos
  pl=ggplot(sub, aes(x=as.factor(get(classes)),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()

write.table(surgery_stats,"surgery_stats.tsv",row.names=FALSE,quote=FALSE,sep="\t")

#summarize surgery_stats
surgery_stats[,test:=factor(test,levels=unique(test[order(category)])),]
surgery_stats[,p.value.top:=ifelse(p.value<=0.05,ifelse(p.value<=0.01,ifelse(p.value<=0.001,"<=0.001","<=0.01"),"<=0.05"),">0.05"),]

pdf(paste0("surgery_stats.pdf"),height=8,width=6)
ggplot(surgery_stats)+geom_tile(aes(x=id,y=test,fill=p.value.top),col="grey")+geom_point(aes(x=0,y=test,col=category),size=4,shape=15)+scale_fill_manual(values=c("<=0.001"="#f03b20","<=0.01"="#feb24c","<=0.05"="#ffeda0",">0.05"="white"))+ylab("")+xlab("Surgery comparison")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

############long_survivor#################################

#long_survivor vs histo_immuno
test_list=column_annotation_combined_clean$histo_immuno
classes="long_survivor"

pdf(paste0("longSurv-histo_immuno.pdf"),height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2,3)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.5*ypos
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(color="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#long_survivor vs histo_segmentation
test_list=grep("Mean |Mode |File|Image|Block| STD ",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)
classes="long_survivor"

pdf(paste0("longSurv-histo_segmentation.pdf"),height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.5*ypos
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#long_survivor vs imaging_segmentation
test_list=grep("ID|surgery",column_annotation_combined_clean$imaging_segmentation,invert=TRUE,value=TRUE)
classes="long_survivor"

pdf(paste0("longSurv-imaging_segmentation.pdf"),height=3,width=2.5)
for (test in test_list){
  sub=unique(annotation[surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  sig=with(sub,test_difference(get(test),get(classes),ntest=length(test_list)))
  ypos=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)
  ypos=ypos+0.5*ypos
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+annotate(geom="text",x=0.5,y=ypos,hjust=0,vjust=1,label=paste0(sig$annot,collapse="\n"))+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()

###########ShiftPhenotype###################################
ShiftPhenotype_stats=data.table()

#ShiftPhenotype vs dip scores
test_list=c("CTCF__Gliobla_None_62750","CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465" ,"CHD1_(A301-218A)__H1-hESC_None_7245"  ,"NRSF__H1-hESC_None_13281" ,"NRSF__U87_None_11741" ,"NANOG__Embryonic Stem Cell_NA_18532"  ,"SOX2__Embryonic Stem Cell_NA_5204","Pol2__Gliobla_None_17405", "Pol2__H1-hESC_None_20311", "EZH2_(39875)__NH-A_None_6488" ,"EZH2__Embryonic Stem Cell_NA_3204")

classes="ShiftPhenotype"

pdf("ShiftPhenotype-dip.pdf",height=3,width=4)
for (test in test_list){
  sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="dip",]
  ShiftPhenotype_stats=rbindlist(list(ShiftPhenotype_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1,]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","ShiftPhenotype")]
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(shape=21,color="grey",alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")
  print(pl)
}
dev.off()

#ShiftPhenotype vs histo_immuno
test_list=column_annotation_combined_clean$histo_immuno
classes="ShiftPhenotype"

pdf(paste0("ShiftPhenotype-histo_immune.pdf"),height=3,width=4)
for (test in test_list){
  sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="histo_immune",]
  ShiftPhenotype_stats=rbindlist(list(ShiftPhenotype_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1,]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","ShiftPhenotype")]
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#ShiftPhenotype vs histo_segmentation
test_list=grep("Mean |Mode |File|Image|Block| STD |mm",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)
classes="ShiftPhenotype"

pdf(paste0("ShiftPhenotype-histo_segmentation.pdf"),height=3,width=4)
for (test in test_list){
  sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="histo_segmentation",]
  ShiftPhenotype_stats=rbindlist(list(ShiftPhenotype_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1,]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","ShiftPhenotype")]
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery)+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#ShiftPhenotype vs imaging_segmentation
test_list=grep("ID|surgery",column_annotation_combined_clean$imaging_segmentation,invert=TRUE,value=TRUE)
classes="ShiftPhenotype"

pdf(paste0("ShiftPhenotype-imaging_segmentation.pdf"),height=3,width=4)
for (test in test_list){
  sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="imaging_segmentation",]
  ShiftPhenotype_stats=rbindlist(list(ShiftPhenotype_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1,]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","ShiftPhenotype")]
  pl=ggplot(sub, aes(x=as.factor(get(classes)),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")
  print(pl)
}
dev.off()


#ShiftPhenotype vs molecular data
test_list=column_annotation_combined_clean$meth_heterogeneity
classes="ShiftPhenotype"

pdf(paste0("ShiftPhenotype-meth_het.pdf"),height=3,width=4)
for (test in test_list){
  sub=unique(annotation[ShiftPhenotype%in%c("classic to sarcoma","stable")&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="meth_het",]
  ShiftPhenotype_stats=rbindlist(list(ShiftPhenotype_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1,]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","ShiftPhenotype")]

  pl=ggplot(sub, aes(x=as.factor(get(classes)),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+facet_wrap(~surgery)+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")
  print(pl)
}
dev.off()

write.table(ShiftPhenotype_stats,"ShiftPhenotype_stats.tsv",row.names=FALSE,quote=FALSE,sep="\t")

#summarize surgery_stats
ShiftPhenotype_stats[,test:=factor(test,levels=unique(test[order(category)])),]
ShiftPhenotype_stats[,p.value.top:=ifelse(p.value<=0.05,ifelse(p.value<=0.01,ifelse(p.value<=0.001,"<=0.001","<=0.01"),"<=0.05"),">0.05"),]

pdf(paste0("ShiftPhenotype_stats.pdf"),height=9,width=5.5)
ggplot(ShiftPhenotype_stats)+geom_tile(aes(x=id,y=test,fill=p.value.top),col="grey")+geom_point(aes(x=0,y=test,col=category),size=4,shape=15)+scale_fill_manual(values=c("<=0.001"="#f03b20","<=0.01"="#feb24c","<=0.05"="#ffeda0",">0.05"="white"))+ylab("")+xlab("Surgery comparison")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

###########progression_types###################################
progression_types_stats=data.table()

#progression_types vs dip scores
test_list=c("CTCF__Gliobla_None_62750","CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465" ,"CHD1_(A301-218A)__H1-hESC_None_7245"  ,"NRSF__H1-hESC_None_13281" ,"NRSF__U87_None_11741" ,"NANOG__Embryonic Stem Cell_NA_18532"  ,"SOX2__Embryonic Stem Cell_NA_5204","Pol2__Gliobla_None_17405", "Pol2__H1-hESC_None_20311", "EZH2_(39875)__NH-A_None_6488" ,"EZH2__Embryonic Stem Cell_NA_3204")

classes="progression_types"

pdf("progression_types-dip.pdf",height=4,width=4)
for (test in test_list){
  sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="dip",]
  progression_types_stats=rbindlist(list(progression_types_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1[1],]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","progression_types")]
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(shape=21,color="grey",alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  print(pl)
}
dev.off()

#progression_types vs histo_immuno
test_list=column_annotation_combined_clean$histo_immuno
classes="progression_types"

pdf(paste0("progression_types-histo_immune.pdf"),height=3,width=4)
for (test in test_list){
  sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="histo_immune",]
  progression_types_stats=rbindlist(list(progression_types_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1[1],]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","progression_types")]
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  print(pl)
}
dev.off()


#progression_types vs histo_segmentation
test_list=grep("Mean |Mode |File|Image|Block| STD |mm",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)
classes="progression_types"

pdf(paste0("progression_types-histo_segmentation.pdf"),height=4,width=4)
for (test in test_list){
  sub=sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="histo_segmentation",]
  progression_types_stats=rbindlist(list(progression_types_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1[1],]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","progression_types")]
  pl=ggplot(sub, aes(x=get(classes),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  print(pl)
}
dev.off()


#progression_types vs imaging_segmentation
test_list=grep("ID|surgery",column_annotation_combined_clean$imaging_segmentation,invert=TRUE,value=TRUE)
classes="progression_types"

pdf(paste0("progression_types-imaging_segmentation.pdf"),height=4,width=4)
for (test in test_list){
  sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="imaging_segmentation",]
  progression_types_stats=rbindlist(list(progression_types_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1[1],]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","progression_types")]
  pl=ggplot(sub, aes(x=as.factor(get(classes)),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+geom_boxplot(outlier.size=NA,fill="transparent")+facet_wrap(~surgery,scale="free")+ylab(test)+xlab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  print(pl)
}
dev.off()


#progression_types vs molecular data
test_list=column_annotation_combined_clean$meth_heterogeneity
classes="progression_types"

pdf(paste0("progression_types-meth_het.pdf"),height=4,width=4)
for (test in test_list){
  sub=unique(annotation[!is.na(progression_types)&surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt",c("N_number_st",test,classes,"surgery"),with=FALSE])
  sig=sub[,test_difference(get(test),get(classes),ntest=length(test_list)),by="surgery"]
  sig[,test:=test,]
  sig[,category:="meth_het",]
  progression_types_stats=rbindlist(list(progression_types_stats,sig))
  sig[,ypos:=max(as.numeric(unlist(sub[,test,with=FALSE])),na.rm=TRUE)*1.1,by="surgery"]
  sig[,eval(classes):=V1[1],]
  sig_red=sig[,list(annot=paste0(annot,collapse="\n")),by=c("surgery","ypos","progression_types")]
  
  pl=ggplot(sub, aes(x=as.factor(get(classes)),y=get(test),group=get(classes)))+geom_point(col="grey",shape=21,alpha=1,size=3,position=position_jitter(width=0.2,height=0))+geom_text(data=sig_red,hjust=0,vjust=1,aes(y=ypos,label=annot))+facet_wrap(~surgery,scale="free")+geom_boxplot(outlier.size=NA,fill="transparent")+ylab(test)+xlab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  print(pl)
}
dev.off()

write.table(progression_types_stats,"progression_types_stats.tsv",row.names=FALSE,quote=FALSE,sep="\t")

#summarize surgery_stats
progression_types_stats[,test:=factor(test,levels=unique(test[order(category)])),]
progression_types_stats[,p.value.top:=ifelse(p.value<=0.05,ifelse(p.value<=0.01,ifelse(p.value<=0.001,"<=0.001","<=0.01"),"<=0.05"),">0.05"),]

pdf(paste0("progression_types_stats.pdf"),height=9,width=5.5)
ggplot(progression_types_stats)+geom_tile(aes(x=id,y=test,fill=p.value.top),col="grey")+geom_point(aes(x=0,y=test,col=category),size=4,shape=15)+scale_fill_manual(values=c("<=0.001"="#f03b20","<=0.01"="#feb24c","<=0.05"="#ffeda0",">0.05"="white"))+ylab("")+xlab("Surgery comparison")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


###########Meth heterogeneity###################################
annotation[,immune_cells:=sum(c(CD163,CD3,CD34,CD45ro,CD68,CD8,CD80,FOXP3),na.rm=TRUE),by=1:nrow(annotation)]
sub=annotation[category=="GBMatch"&IDH=="wt"&surgery%in%c(1,2)&enrichmentCycles>12&enrichmentCycles<16]

factors=c("Area total [mm²]","Area tumor [mm²]","Total (mm3)","Enhancing (mm3)","Relative share tumor","immune_cells")

pdf(paste0("heterogeneity_size.pdf"),height=2.5,width=4.5)
for (factor in factors){
cor=sub[,cor(get(factor),mean_entropy,use="pairwise.complete.obs"),by=surgery]
x=sub[,max(get(factor),na.rm=TRUE)/2,]
ggpl=ggplot(sub,aes(x=get(factor),y=mean_entropy))+geom_point(shape=21,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cor,aes(x=x,y=45,label=paste0("r=",round(V1,3))))+facet_wrap(~surgery,scale="free")+xlab(factor)
print(ggpl)

cor=sub[,cor(get(factor),mean_pdr,use="pairwise.complete.obs"),by=surgery]
x=sub[,max(get(factor),na.rm=TRUE)/2,]
ggpl=ggplot(sub,aes(x=get(factor),y=mean_pdr*100))+geom_point(shape=21,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cor,aes(x=x,y=30,label=paste0("r=",round(V1,3))))+facet_wrap(~surgery,scale="free")+xlab(factor)
print(ggpl)
}
dev.off()

