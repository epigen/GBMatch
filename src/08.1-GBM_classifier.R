library(project.init)
project.init2("GBMatch")
library(simpleCache)
library(riverplot)
library(pheatmap)

source(file.path(getOption("PROJECT.DIR"),"src/99-liblinearFunctions.R"))

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.1-GBM_classifier")
dir.create(out_dir)
setwd(out_dir)

annot=fread(file.path(extData,"extData/TCGA/annotation_Brennanetal.txt"))

annot[,beta_file:=head(list.files(path=file.path(extData,"extData/TCGA/"),pattern=paste0(".*",`Case ID`,".*hg38.txt"),recursive = TRUE,full.names=TRUE),n=1),by=1:nrow(annot)]

#load annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))



########Model training and testing on microarray data########################################

annot_sub=annot[grepl("beta_27",beta_file)&`IDH1 status`=="WT"&`G-CIMP methylation`=="non-G-CIMP"&`Expression Subclass`!="Neural"]
sub="27_noNeural"

collect_beta=function(beta_annot){
for(i in (1:nrow(beta_annot))){
  if (i==1){
    beta_mat=fread(beta_annot[i]$beta_file,select = c(1:2))
    setnames(beta_mat,names(beta_mat),c("id",beta_annot[i]$`Case ID`))
  }else{
    temp=fread(beta_annot[i]$beta_file,select = c(1:2))
    setnames(temp,names(temp),c("id",beta_annot[i]$`Case ID`))
    beta_mat=merge(beta_mat,temp,by="id")
  }
}
return(beta_mat)
}

beta_mat=simpleCache(paste0("beta_",sub),collect_beta(beta_annot=annot_sub),cacheSubDir="GBM_classifier",recreate=FALSE)

beta_mat_t=t(beta_mat[,-c(1),with=FALSE])
colnames(beta_mat_t)=beta_mat$id
na_cols=apply(beta_mat_t,2,function(x){any(is.na(x))})
beta_mat_t=beta_mat_t[,!na_cols]

#get array annotation
array_loc=fread(annot_sub$beta_file[1],select = c(1,3,4,5))
setnames(array_loc,names(array_loc),c("id","chr","start","end"))

table(row.names(beta_mat_t)==annot_sub$`Case ID`)

data=beta_mat_t
labels=annot_sub$`Expression Subclass`
samples=annot_sub$`Case ID`


pl=check_prediction(data=data,labels=labels,samples=samples,type=0,cost=1)

pdf(file.path(paste0("roc_auc_",sub,".pdf")),height=3,width=4)
pl$plot
dev.off()


################Classification of GBMatch###############################################################
#GBMatch
loadCaches("rrbsCg")
setnames(rrbsCg,"id","sample")
array_loc=fread(annot_sub$beta_file[1],select = c(1,3,4,5))
setnames(array_loc,names(array_loc),c("id","chr","start","end"))
rrbs_annot=merge(rrbsCg,array_loc,by=c("chr","start"))

#test if classification still works only on the RRBS covered cpgs
data=beta_mat_t[,colnames(beta_mat_t)%in%unique(rrbs_annot$id)]
labels=annot_sub$`Expression Subclass`
samples=annot_sub$`Case ID`

pl=check_prediction(data=data,labels=labels,samples=samples,type=0,cost=1)

pdf(file.path(paste0("roc_auc_",sub,"inRRBS.pdf")),height=3,width=4)
pl$plot
dev.off()

#now create model for each RRBS sample individually, test it and predict
predict_RRBS=function(cpgID,RRBS_meth,type=NULL,cost=NULL,train_mat,scaleAndCenter=FALSE){
  data=train_mat[,colnames(train_mat)%in%cpgID]
  labels=annot_sub$`Expression Subclass`
  samples=annot_sub$`Case ID`
  
  #build and check model
  res=check_prediction(data=data,labels=labels,samples=samples,type=type,scaleAndCenter=scaleAndCenter,cost=cost)
  #order cpgs to match the order of the model
  RRBS_meth=RRBS_meth[cpgID%in%colnames(data)]
  cpgID=cpgID[cpgID%in%colnames(data)]
  RRBS_data=matrix(RRBS_meth[order(match(cpgID,grep("Bias",colnames(res$model$W),invert=TRUE,value = TRUE)))],nrow = 1)
  RRBS_data_names=cpgID[order(match(cpgID,grep("Bias",colnames(res$model$W),invert=TRUE,value = TRUE)))]
  colnames(RRBS_data)=RRBS_data_names
  
  pr=FALSE
  dv=FALSE
  if(res$model$Type==0 || res$model$Type==7) pr=TRUE
  if(res$model$Type%in%c(0:7)) dv=TRUE
  
  if (scaleAndCenter==TRUE){
    p=predict(res$model,scale(RRBS_data,res$attr_center,res$attr_scale),proba=pr,decisionValues=dv)
    p_val=predict(res$model,scale(data,res$attr_center,res$attr_scale),proba=pr,decisionValues=dv)
  }else{
    p=predict(res$model,RRBS_data,proba=pr,decisionValues=dv)
    p_val=predict(res$model,data,proba=pr,decisionValues=dv)
  }
  
  #see how model performs on same data that were used for training  
  val=table(p_val$predictions==labels)
  val_rat=max(na.rm=TRUE,c(val["TRUE"],0))/(max(na.rm=TRUE,c(val["TRUE"],0))+max(na.rm=TRUE,c(val["FALSE"],0)))
  
  return(list(val_rat=val_rat,prediction=p,check=res,auc=res$auc,auc_rand=res$auc_rand))
}

#warning probably due to column check being of length 4 and column prediction being of length 3
#beta without scaling
simpleCache(paste0("res_",sub,"_predRRBS_notScaled"),"rrbs_annot[,predict_RRBS(cpgID=id,RRBS_meth=methyl,type = 0,cost=1,train_mat=beta_mat_t),by=c('sampleName','sample')]", cacheSubDir="GBM_classifier")


#use not scaled beta values
test_res_prob=res_27_noNeural_predRRBS_notScaled

#calculate CpGs used for each prediction
CpG_tab=test_res_prob[,ncol(check[[3]]$W),by="sampleName"]
write.table(as.matrix(summary(CpG_tab$V1)),file=paste0("used_features_",sub,"_predRRBS.tsv"),sep="\t",quote=FALSE)

#plot ROC curves
print_res=function(plots,title){
  for (plot in plots){
    print(plot+ggtitle(title)) 
  }
}

pdf(file.path(paste0("roc_auc_",sub,"_predRRBS.pdf")),height=4,width=4.5)
test_res_prob[,print_res(check[1][[1]],paste0(sample[1]," ",prediction[[1]][1]," prob:",round(prediction[[2]][which(colnames(prediction[[2]])==prediction[[1]][1])],3)," auc:",round(auc,2))),by="sampleName"]
dev.off()

#make sure to analyse only samples that are in the annotation
test_res_prob=test_res_prob[sample%in%annotation$N_number_seq]


class_probs=test_res_prob[,list(sub_group=prediction[[1]][1],prob=prediction[[2]][which(colnames(prediction[[2]])==prediction[[1]][1])],auc=auc[1],auc_rand=auc_rand[1],val_rat=val_rat[1]),c("sampleName","sample")]
write.table(class_probs,paste0("class_probs",sub,"_predRRBS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)


class_probs_all=test_res_prob[,list(sub_group=prediction[[1]][1],Classical=prediction[[2]][1],Proneural=prediction[[2]][2],Mesenchymal=prediction[[2]][3],auc=auc[1],auc_rand=auc_rand[1]),c("sampleName","sample")]

class_probs_all_annot=class_probs_all[annotation[,c("N_number_seq","patID","category","cohort","IDH","surgery.x","WHO2016_classification","Follow-up_years"),with=FALSE], on=c(sample="N_number_seq")]

class_probs_all_annot[is.na(sub_group)]
class_probs_all_annot=class_probs_all_annot[!is.na(sub_group)]

class_probs_all_annot[,sub_group:=as.character(sub_group),]
class_probs_all_annot[,sub_group_prob:=get(sub_group),by=1:nrow(class_probs_all_annot)]
class_probs_all_annot[,switching:=ifelse(length(unique(sub_group[category=="GBMatch"&surgery.x%in%c(1,2)]))==1,FALSE,TRUE),by="patID"]

write.table(class_probs_all_annot,paste0("class_probs_annot_",sub,"_predRRBS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

pdf(file.path(paste0("class_probs_",sub,"_predRRBS.pdf")),height=3,width=5)
ggplot(class_probs_all_annot[category%in%c("GBMatch","GBmatch_val")],aes(y=sub_group_prob,x=substr(sub_group,1,3),group=sub_group))+geom_point(aes(fill=auc,col=auc),position=position_jitter(width=0.2,height=0),alpha=0.5,pch=21)+ylim(c(0,1))+geom_boxplot(fill="transparent")+scale_fill_continuous(low="blue",high="red")+scale_color_continuous(low="blue",high="red")+facet_wrap(~cohort)

ggplot(class_probs_all_annot[category%in%c("GBMatch","GBmatch_val")],aes(y=auc,x=substr(sub_group,1,3),group=sub_group))+geom_point(aes(fill=sub_group,col=sub_group),position=position_jitter(width=0.2,height=0),alpha=0.2,pch=21)+ylim(c(0,1))+geom_boxplot(outlier.shape=NA,fill="transparent")+geom_hline(yintercept=0.8,lty=20,col="grey")+xlab("")+ylab("ROC AUC")+facet_wrap(~cohort)
dev.off()


class_probs_all_annot_long=melt(class_probs_all_annot,id.vars=c("sampleName","sample","patID","category","sub_group","sub_group_prob","switching","auc","auc_rand","surgery.x","WHO2016_classification","Follow-up_years","IDH"),value.name="prob")


class_probs_all_annot_long[,variable:=factor(variable,levels=c("Classical","Mesenchymal","Proneural")),]
pdf(file.path(paste0("multisel_",sub,"_predRRBS.pdf")),height=3,width=8.5)
ggplot(class_probs_all_annot_long[category=="multiselector"&!is.na(prob)],aes(x=paste0(surgery.x,"_",gsub(".*[_BCA675]","",sample)),y=prob,fill=variable,alpha=auc))+geom_bar(stat="identity")+facet_grid(~.~patID,space="free",scale="free")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+xlab("")+ylab("Class probability")+scale_alpha_continuous(range=c(0.35,1))
dev.off()

#class_probs_all_annot_long[,sub_group_prob:=prob[sub_group==variable],by="sample"]
class_probs_all_annot_long[,sample:=factor(sample,levels=unique(sample[order(sub_group_prob,decreasing=TRUE)])),]
#class_probs_all_annot_long[,switching:=ifelse(length(unique(sub_group[category=="GBMatch"&surgery.x%in%c(1,2)]))==1,FALSE,TRUE),by="patID"]

pdf(file.path(paste0("class_probs_stack",sub,"_predRRBS_prim.pdf")),height=3,width=6)
ggplot(class_probs_all_annot_long[!is.na(prob)&category=="GBMatch"&IDH=="wt",],aes(x=sample,y=prob,fill=variable))+geom_bar(stat="identity",position="stack")+facet_grid(.~sub_group,scales="free",space="free")+ ylab("Class probability")+theme(axis.text.x = element_blank(),axis.ticks = element_blank())
dev.off()

pdf(file.path(paste0("class_probs_stack",sub,"_predRRBS_bySurg_prim.pdf")),height=3,width=7)
ggplot(class_probs_all_annot_long[!is.na(prob)&category=="GBMatch"&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=sample,y=prob,fill=variable))+geom_bar(stat="identity",position="stack")+facet_grid(~sub_group+surgery.x,scales="free",space="free")+ ylab("Class probability")+theme(axis.text.x = element_blank(),axis.ticks = element_blank())
dev.off()

pdf(file.path(paste0("class_probs_stack",sub,"_predRRBS_val.pdf")),height=3,width=6)
ggplot(class_probs_all_annot_long[!is.na(prob)&category=="GBmatch_val"&IDH=="wt",],aes(x=sample,y=prob,fill=variable))+geom_bar(stat="identity",position="stack")+facet_grid(.~sub_group,scales="free",space="free")+ ylab("Class probability")+theme(axis.text.x = element_blank(),axis.ticks = element_blank())
dev.off()

#compare subtype composition between primary and validation cohort
pdf(file.path(paste0("sub_group_composition",sub,"_predRRBS.pdf")),height=3,width=4.5)
ggplot(class_probs_all_annot[!is.na(sub_group)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=paste0(category,"\n","Surgery:",surgery.x),fill=sub_group))+geom_bar()+xlab("")+ggtitle("All")
ggplot(class_probs_all_annot[!is.na(sub_group)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=paste0(category,"\n","Surgery:",surgery.x),fill=sub_group))+geom_bar(position="fill")+xlab("")+ggtitle("All")

ggplot(class_probs_all_annot[!is.na(sub_group)&auc>0.8&sub_group_prob>0.65&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=paste0(category,"\n","Surgery:",surgery.x),fill=sub_group))+geom_bar()+xlab("")+ggtitle("High fidelity assignments")
ggplot(class_probs_all_annot[!is.na(sub_group)&auc>0.8&sub_group_prob>0.65&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=paste0(category,"\n","Surgery:",surgery.x),fill=sub_group))+geom_bar(position="fill")+xlab("")+ggtitle("High fidelity assignments")
dev.off()



#for IDH mut
pdf(file.path(paste0("class_probs_stack_IDHmut",sub,"_predRRBS.pdf")),height=3,width=7.5)
ggplot(class_probs_all_annot_long[!is.na(prob)&IDH=="mut"&category%in%c("GBMatch","GBmatch_add","GBmatch_val")],aes(x=paste0(patID,"_",surgery.x),y=prob,fill=variable,alpha=auc))+geom_bar(stat="identity",position="stack")+facet_grid(~.~patID,scales="free",space="free")+ ylab("Class probability")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+scale_alpha_continuous(range=c(0.6,1))
dev.off()

class_probs_all_annot_stat=class_probs_all_annot[surgery.x%in%c(1,2)&IDH=="mut"&category%in%c("GBMatch","GBmatch_add","GBmatch_val"),.N,by=c("sub_group","surgery.x")]
pdf(file.path(paste0("class_probs_pie_IDHmut",sub,"_predRRBS.pdf")),height=3,width=3.5)
ggplot(class_probs_all_annot_stat[surgery.x==1],aes(x="", y=N, fill=sub_group))+geom_bar(stat="identity",width=1)+ ylab("Class probability")+ coord_polar("y", start=0)+ggtitle("Primary tumor")+scale_fill_manual(values=c("Classical"="#F8766D","Mesenchymal"="#00BA38","Proneural"="#619CFF"))
ggplot(class_probs_all_annot_stat[surgery.x==2],aes(x="", y=N, fill=sub_group))+geom_bar(stat="identity",width=1)+ ylab("Class probability")+ coord_polar("y", start=0)+ggtitle("Relapsed tumor")+scale_fill_manual(values=c("Classical"="#F8766D","Mesenchymal"="#00BA38","Proneural"="#619CFF"))
dev.off()

pdf(file.path(paste0("class_probs_switch",sub,"_predRRBS.pdf")),height=3,width=6)
ggplot(class_probs_all_annot[!is.na(sub_group)&auc>0.8&category=="GBMatch"&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=surgery.x,y=sub_group_prob,group=patID))+geom_line()+geom_point(aes(fill=sub_group),pch=21)+facet_wrap(~switching)
dev.off()

pdf(file.path(paste0("sub_group_follow-up",sub,"_predRRBS.pdf")),height=3,width=9)
ggplot(class_probs_all_annot[!is.na(sub_group)&auc>0.8&sub_group_prob>0.65&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2),],aes(x=sub_group,y=`Follow-up_years`))+geom_boxplot(outlier.shape=NA)+geom_point(aes(fill=sub_group),alpha=0.6,position=position_jitter(width=0.25),pch=21)+facet_wrap(surgery.x~category)
dev.off()


#transit analysis (river plot)
class_probs_all_annot_long_transit=unique(class_probs_all_annot_long[surgery.x%in%c(1,2)&category=="GBMatch"&IDH=="wt",list(auc.1=auc[surgery.x==1],auc.2=auc[surgery.x==2],prob.1=sub_group_prob[surgery.x==1],prob.2=sub_group_prob[surgery.x==2],transit=paste0(sub_group[surgery.x==1],"_",sub_group[surgery.x==2])),by=patID])


auc_trs=0.8
transits=class_probs_all_annot_long_transit[auc.1>auc_trs&auc.2>auc_trs,.N,by=transit]

transits[,N1:=gsub(".*_","",transit),]
transits[,N2:=gsub("_.*","",transit),]
transits=transits[!N1==""&!N2==""]

edges=list()
for (type in unique(transits$N1)){
  for (trans in unique(transits[N1==type]$N2)){
    edges[[paste0(type,".1")]][paste0(trans,".2")]=transits[N1==type&N2==trans]$N
    
  }
}

r <- makeRiver(nodes=sort(c(names(edges), names(edges[[1]]))),edges=edges,node_xpos=c(1,2,1,2,1,2),node_ypos=c(1,1,2,2,3,3),node_styles= list( Mesenchymal.2= list( col= "#00BA38" ),Mesenchymal.1= list( col= "#00BA38" ),Classical.2= list( col= "#F8766D" ),Classical.1= list( col= "#F8766D" ),Proneural.2= list( col= "#619CFF" ),Proneural.1= list( col= "#619CFF"  )))

annot=merge(as.data.table(r$edges),setnames(as.data.table(r$nodes),"ID","N1"),by="N1",all.x=TRUE)
annot=merge(annot,setnames(as.data.table(r$nodes),"ID","N2"),by="N2",all.x=TRUE)
annot[,total.1:=sum(Value),by=N1]
annot[,ratio:=Value/total.1,]

annot_1=annot[,sum(Value),by=c("N1","x.x","y.x")]
annot_2=annot[,sum(Value),by=c("N2","x.y","y.y")]
annot_1[,x:=1.05,]
annot_2[,x:=1.95,]
annot_all=rbindlist(list(annot_1,annot_2))

pdf(file.path(paste0("river_auc_",auc_trs,"_predRRBS.pdf")),height=6,width=6)
plot(r)
text(annot_all$x,annot_all$y.x,labels=annot_all$V1)
dev.off()

