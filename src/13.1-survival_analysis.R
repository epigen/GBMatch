#NOTE: This script performs survival analysis using Kaplan-Meier plots and the log-rank test.
library(project.init)
project.init2("GBMatch")
library(survival)
library(broom)
library("survminer")
require(bit64)

##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/13.1-survival/")
dir.create(out_dir)
setwd(out_dir)

#functions to binarize continuous covariates
binarize=function(vector,quantile){
  ret=ifelse(vector<quantile(vector,quantile,na.rm=TRUE),"low","high")
  return(as.character(ret))
}

binarize2=function(vector,low,high){
  ret=ifelse(vector<=quantile(vector,low,na.rm=TRUE),"low",ifelse(vector>=quantile(vector,high,na.rm=TRUE),"high",NA))
  return(as.character(ret))
}


##function for survival plotting
plot_surv=function(data,data_name,type="surv"){
  pl=list()
  pl2=list()
  ret=data.table()
  #function for p-value calculation
  get_pval=function (fit) 
  {
    if (length(levels(summary(fit)$strata)) == 0) 
      return(NULL)
    sdiff <- survival::survdiff(eval(fit$call$formula), data = eval(fit$call$data),rho=0) #rho=0 means log-rank test
    pvalue <- stats::pchisq(sdiff$chisq, length(sdiff$n) - 1,lower.tail = FALSE)
    return(pvalue)
  }
  
  categories=grep("patID|event|follow_up|category|annoclass",names(data),invert=TRUE,value=TRUE)
  
  source_data=data.table()
  stats=data.table()
  for (category in categories){
    print(category)
    sub=na.omit(data[,grepl(paste0("patID|event|follow_up|^",category,"$"),names(data)),with=FALSE])
    if (length(unique(unlist(sub[,category,with=FALSE])))<2){
      print("Too fiew categories!")
      next
    }
    survfit_1=survfit(Surv(time = follow_up, event = event) ~ get(category), data=sub)
    pvalue=signif(get_pval(survfit_1),2)
    print(pvalue)
    
    #prepare source_data
    sub[,category_name:=category,]
    sub[,analysis:=data_name,]
    sub[,analysis_type:=type,]
    setnames(sub,category,"category")
    source_data=rbindlist(list(source_data,sub))
    
    #calculate sample number
    sample_N=sub[,.N,by="category"]
    sample_N[,annot:=paste0("N(",category,")=",N),]
    
    #collect stats
    forStats=sub[,list(N=.N,events=sum(event),mean_follow_up=signif(mean(follow_up),3),sd_follow_up=signif(sd(follow_up),3)),by=category]
    stats=rbindlist(list(stats,data.table(V1=forStats$category[1],V2=forStats$category[2],V3=forStats$category[3],N_V1=forStats$N[1],N_V2=forStats$N[2],N_V3=forStats$N[3],events_V1=forStats$events[1],events_V2=forStats$events[2],events_V3=forStats$events[3],mean_follow_up_V1=forStats$mean_follow_up[1],mean_follow_up_V2=forStats$mean_follow_up[2],mean_follow_up_V3=forStats$mean_follow_up[3],sd_follow_up_V1=forStats$sd_follow_up[1],sd_follow_up_V2=forStats$sd_follow_up[2],sd_follow_up_V3=forStats$sd_follow_up[3],p.value=pvalue,category_name=category,analysis=data_name,analysis_type=type  )))
    
    if (type=="surv"){ylab="Overall survival probability"}else if (type=="rel"){ylab="Progression-free survival probability"}
    pl[category]=ggsurvplot(fit=survfit_1,conf.int=TRUE,main=paste0(category,"\np-value=",pvalue,"\n",paste0(sample_N$annot,collapse=" ")),palette="Set2",font.main=c(10, "plain", "black"),ylab = ylab,xlab="Months")
    ret=rbindlist(list(ret,data.table(category=category,p.value=pvalue)))
    
  }
  write.table(source_data,paste0(data_name,"_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  write.table(stats,paste0(data_name,"_stats.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  return(list(pl,ret))
}

plot_dots=function(data){
  data_long=melt(data,id.vars=c("event", "follow_up","patID","category"))
  pl=list()
  for(sub in unique(data_long$variable)){
    if(sum(!is.na(data_long[variable==sub]$value))==0){next}
    dotplot=ggplot(data_long[variable==sub&!is.na(value)],aes(x=value,y=follow_up))+geom_point(shape=21,size=2.5,aes(fill=value),alpha=0.6,position=position_jitter(width=0.2))+geom_boxplot(fill="transparent",outlier.size=NA)+ylab("Months")+scale_fill_brewer(palette="Set2")+xlab(sub)
    pl[[sub]]=dotplot
  }
  return(pl)
}

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined_final.tsv"))
load(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/column_annotation_combined.RData"))

annot_ASC=annotation[IDH=="wt",list(cycles.1=enrichmentCycles[surgery==1][1],cycles.2=enrichmentCycles[surgery==2][1],Age=min(Age),Sex=unique(Sex),annoclass="cycles"),by=patID]

annot_surv=annotation[surgery%in%c(1)&IDH=="wt"&category%in%c("GBMatch","GBmatch_val"),list(event=unique(VitalStatus),follow_up=unique(`Follow-up_years`)*12,category=unique(category),annoclass="survival"),by=patID]
annot_surv[,event:=ifelse(event=="dead",1,0)]

annot_relapse=unique(annotation[surgery%in%c(1,2)&IDH=="wt"&category%in%c("GBMatch","GBmatch_val"),list(event=1,follow_up=(unique(timeToFirstProg)/365)*12,category=category,annoclass="relapse"),by=patID])

annot_relapse_est=unique(annotation[surgery%in%c(1,2)&IDH=="wt"&category%in%c("GBMatch","GBmatch_val"),list(event=ifelse(is.na(timeToFirstProg),unique(VitalStatus),"dead"),follow_up=ifelse(is.na(timeToFirstProg),unique(`Follow-up_years`)*12-2,(unique(timeToFirstProg)/365)*12),category=category,annoclass="relapse"),by=patID])
annot_relapse_est[,event:=ifelse(event=="dead",1,0)]

combi=rbindlist(list(annot_surv,annot_relapse),use.names=TRUE)
combi[,long_surv:=ifelse(follow_up[annoclass=="survival"]>=36,TRUE,FALSE),by="patID"]

#check if survival time is longer that time to relapse. Should be TRUE or NA
surv_rel_check=combi[,round(follow_up[annoclass=="survival"],2)>=round(follow_up[annoclass=="relapse"],2),by="patID"]
surv_rel_check[V1==FALSE]

surv_stats=combi[,round(median(follow_up,na.rm=TRUE),2),by=c("category","annoclass")]

#survival distribution
pdf("follow_up_overview.pdf",height=4,width=6)
ggplot(combi,aes(x=annoclass,y=follow_up))+geom_point(size=2.5,shape=21,alpha=0.6,aes(fill=long_surv),position=position_jitter(width=0.3))+geom_boxplot(outlier.size=NA,fill="transparent")+annotate(geom="text",x=1.5,y=114,label=paste0("Median survival: ",surv_stats[annoclass=="survival"]$V1,"\n","Median time to relapse: ",surv_stats[annoclass=="relapse"]$V1))+scale_fill_manual(values=c("FALSE"="grey","TRUE"="red"))+scale_y_continuous(breaks=seq(from=0,to=125,by=6))+xlab("")+ylab("months")+facet_wrap(~category)
dev.off()


################################
##prom diff-meth stratification
################################
prom_diff=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/11.2-diffMeth_single/sel_recurring_trend_pat.tsv"))
prom_diff_enrichr_term=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/11.2-diffMeth_single/promoter_diff_meth_enrichr_term.tsv"))

prom_diff_combined=merge(prom_diff[Ngenes>0],prom_diff_enrichr_term,by="patID",all=TRUE)

prom_diff_bin=prom_diff_combined[,list(patID=patID,dist_bin=binarize2(diff_trend_dist_norm,0.2,0.8),trend_thres=ifelse(diff_trend_dist_norm<=0.7,"trend",ifelse(diff_trend_dist_norm>=1.15,"anti-trend",NA)),trend_thres_3=ifelse(diff_trend_dist_norm<=0.7,"trend",ifelse(diff_trend_dist_norm>=1.15,"anti-trend","moderate")),trend_split=ifelse(diff_trend_dist_norm<1,"trend",ifelse(diff_trend_dist_norm>=1,"anti-trend",NA)),mean_diffmeth.Development=binarize2(mean_diffmeth.Development,0.3,0.7),`mean_diffmeth.Apoptosis`=binarize2(`mean_diffmeth.Apoptosis`,0.3,0.7),`mean_diffmeth.Wnt signalling`=binarize2(`mean_diffmeth.Wnt signalling`,0.3,0.7),`mean_diffmeth.Immune response`=binarize2(`mean_diffmeth.Immune response`,0.3,0.7)),]

prom_diff_bin[,`mean_diffmeth.Wnt signalling`:=factor(`mean_diffmeth.Wnt signalling`,levels=c("low","high")),]

prom_diff_surv=merge(prom_diff_bin,annot_surv,by="patID")
prom_diff_relapse=merge(prom_diff_bin,annot_relapse,by="patID")

####Figure 6d,f
pdf("prom_diff_annotation_surv.pdf",height=4,width=3)
print(plot_surv(prom_diff_surv,"prom_diff_surv")[[1]])
dev.off()

pdf("prom_diff_annotation_surv_dp.pdf",height=4,width=3)
print(plot_dots(prom_diff_surv))
dev.off()

####Figure 6d,f
pdf("prom_diff_annotation_relapse.pdf",height=4,width=3)
print(plot_surv(prom_diff_relapse,"prom_diff_relapse",type="rel")[[1]])
dev.off()

pdf("prom_diff_annotation_relapse_dp.pdf",height=4,width=3)
print(plot_dots(prom_diff_relapse))
dev.off()



################################
##methclone stratification
################################
methclone=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/06-methclone/summary_minReads40/EPM_epy_cutoff-80/EPM_1vs2.tsv"))
methclone=methclone[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"]
setnames(methclone,"patient","patID")

#for 1vs2
methclone[cycles.1>=16|cycles.1<=12,entropy.1:=NA,]
methclone[cycles.2>=16|cycles.2<=12,entropy.2:=NA,]
methclone[cycles.2>=16|cycles.2<=12|cycles.1>=16|cycles.1==12,mean_dentropy:=NA,]
methclone[cycles.2>=16|cycles.2<=12|cycles.1>=16|cycles.1==12,EPM:=NA,]

entropy=methclone[,list(patID=patID,entropy_1=binarize2(mean_entropy1,0.2,0.8),entropy_2=binarize2(mean_entropy2,0.2,0.8),d_entropy=binarize2(mean_dentropy,0.5,0.5),EPM= binarize2(EPM,0.5,0.5)),by="category"]

entropy_surv=merge(entropy,annot_surv,by=c("patID","category"))
entropy_relapse=merge(entropy,annot_relapse,by=c("patID","category"))
entropy_relapse_est=merge(entropy,annot_relapse_est,by=c("patID","category"))


categories=c("GBMatch")

for (sel_category in categories){
  ####Figure S12e
  pdf(paste0("methclone_annotation_surv_",sel_category,".pdf"),height=3.5,width=3.0)
  print(plot_surv(entropy_surv[category==sel_category],paste0("methclone_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("methclone_annotation_surv_dp_",sel_category,".pdf"),height=2,width=2)
  print(plot_dots(entropy_surv[category==sel_category]))
  dev.off()
  
  ####Figure S12e
  pdf(paste0("methclone_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3.0)
  print(plot_surv(entropy_relapse[category==sel_category],paste0("methclone_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("methclone_annotation_relapse_dp_",sel_category,".pdf"),height=2,width=2)
  print(plot_dots(entropy_relapse[category==sel_category]))
  dev.off()
  
  pdf(paste0("methclone_annotation_relapse_est_",sel_category,".pdf"),height=3.5,width=3.0)
  print(plot_surv(entropy_relapse_est[category==sel_category],paste0("methclone_annotation_relapse_est_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("methclone_annotation_relapse_est_dp_",sel_category,".pdf"),height=2,width=2)
  print(plot_dots(entropy_relapse_est[category==sel_category]))
  dev.off()
}


########################################
#combined heterogeineity stratification
########################################
simpleCache("combined_heterogeneity",assignToVariable="combined_heterogeneity")
setnames(combined_heterogeneity,"id","N_number_seq")

het_annot=merge(combined_heterogeneity,annotation[,c("category","N_number_seq","patID","surgery","IDH","enrichmentCycles"),with=FALSE],by="N_number_seq")
het_annot=het_annot[category%in%c("GBMatch","GBmatch_val")&surgery%in%c(1,2)&IDH=="wt"]

het_annot[enrichmentCycles>=16|enrichmentCycles<=12,mean_entropy:=NA,]
het_annot[enrichmentCycles>=16|enrichmentCycles<=12,mean_pdr:=NA,]
het_annot[enrichmentCycles>=16|enrichmentCycles<=12,enrichmentCycles:=NA,]

het_bin=het_annot[,list(patID=patID,cycles=binarize2(enrichmentCycles,0.2,0.8),mean_entropy=binarize2(mean_entropy,0.2,0.8),mean_pdr=binarize2(mean_pdr,0.2,0.8)),by=c("surgery","category")]

het_wide=reshape(het_bin,idvar=c("patID","category"),timevar="surgery",direction="wide")


het_surv=merge(het_wide,annot_surv,by=c("patID","category"))
het_relapse=merge(het_wide,annot_relapse,by=c("patID","category"))
het_relapse_est=merge(het_wide,annot_relapse_est,by=c("patID","category"))

categories=c("GBMatch","GBmatch_val")

for (sel_category in categories){
  ####Figure 5d; S12c
  pdf(paste0("heterogeneity_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(het_surv[category==sel_category],paste0("heterogeneity_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("heterogeneity_annotation_surv_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(het_surv[category==sel_category]))
  dev.off()
  
  ####Figure 5d
  pdf(paste0("heterogeneity_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(het_relapse[category==sel_category],paste0("heterogeneity_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("heterogeneity_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(het_relapse[category==sel_category]))
  dev.off()
  
  ####Figure S12c
  pdf(paste0("heterogeneity_annotation_relapse_est_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(het_relapse_est[category==sel_category],paste0("heterogeneity_annotation_relapse_est_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("heterogeneity_annotation_relapse_est_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(het_relapse_est[category==sel_category]))
  dev.off() 
}


###########################
##bissnp stratification
###########################
bissnp=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/04-bissnp/bissnp_var_pat.tsv"))
bissnp_norm=bissnp[surgery%in%c(1,2),list(patID=patID,category=category,surgery=surgery,all_count=all_count/bg_calls*1000000,H_count=H_count/bg_calls*1000000,M_count=M_count/bg_calls*1000000),]

bissnp_bin=bissnp_norm[patID%in%annotation[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"]$patID,list(patID=patID,surgery=surgery,all_count=binarize2(all_count,0.2,0.8),H_count=binarize2(H_count,0.2,0.8),M_count=binarize2(M_count,0.2,0.8)),by=c("surgery","category")]

bissnp_wide=reshape(bissnp_bin[surgery%in%c(1,2)],timevar="surgery",idvar=c("patID","category"),direction="wide")
bissnp_wide[,switch:=ifelse(H_count.1==H_count.2|(is.na(H_count.1)&is.na(H_count.2)),FALSE,TRUE),]
bissnp_wide[,toHigh:=ifelse((H_count.1!="high"|is.na(H_count.1))&H_count.2=="high",TRUE,FALSE),]
bissnp_wide[,toLow:=ifelse((H_count.1!="low"|is.na(H_count.1))&H_count.2=="low",TRUE,FALSE),]

bissnp_surv=merge(bissnp_wide,annot_surv,by=c("patID","category"))
bissnp_relapse=merge(bissnp_wide,annot_relapse,by=c("patID","category"))

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  pdf(paste0("bissnp_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(bissnp_surv[category==sel_category],paste0("bissnp_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("bissnp_annotation_surv_dp_",sel_category,".pdf"),height=2,width=2)
  print(plot_dots(bissnp_surv[category==sel_category]))
  dev.off()
  
  
  pdf(paste0("bissnp_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(bissnp_relapse[category==sel_category],paste0("bissnp_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("bissnp_annotation_relapse_dp_",sel_category,".pdf"),height=2,width=2)
  print(plot_dots(bissnp_relapse[category==sel_category]))
  dev.off()
}

###########################
##transcriptional subtype stratification
###########################
subtypes=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.1-GBM_classifier/class_probs_annot_27_noNeural_predRRBS.tsv"))
subtypes=subtypes[!is.na(sub_group)]
subtypes[,sub_type_prob:=get(sub_group),by=1:nrow(subtypes)]
subtypes[,isMesenchymal:=ifelse(sub_group=="Mesenchymal",TRUE,FALSE),]
subtypes[,isClassical:=ifelse(sub_group=="Classical",TRUE,FALSE),]
subtypes[,isProneural:=ifelse(sub_group=="Proneural",TRUE,FALSE),]

subtypes_wide=reshape(subtypes[order(auc,decreasing=TRUE)][category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2)&auc>0.8,c("surgery.x","patID","category","sub_group","isMesenchymal","isClassical","isProneural"),with=FALSE],timevar="surgery.x",idvar=c("patID","category"),direction="wide")
subtypes_wide[,switch:=ifelse(sub_group.1==sub_group.2,FALSE,TRUE),]
subtypes_wide[,toMesenchymal:=ifelse(sub_group.1!="Mesenchymal"&sub_group.2=="Mesenchymal",TRUE,FALSE),]
subtypes_wide[,toProneural:=ifelse(sub_group.1!="Proneural"&sub_group.2=="Proneural",TRUE,FALSE),]
subtypes_wide[,toClassical:=ifelse(sub_group.1!="Classical"&sub_group.2=="Classical",TRUE,FALSE),]

subtypes_surv=merge(subtypes_wide,annot_surv,by=c("patID","category"))
subtypes_relapse=merge(subtypes_wide,annot_relapse,by=c("patID","category"))

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  
  ####Figure 2e; S6h
  pdf(paste0("subtype_annotation_surv_",sel_category,".pdf"),height=4,width=4)
  print(plot_surv(subtypes_surv[category==sel_category],paste0("subtype_annotation_surv_",sel_category),type="surv")[[1]])
  dev.off()
  
  pdf(paste0("subtype_annotation_surv_dp_",sel_category,".pdf"),height=4,width=4)
  print(plot_dots(subtypes_surv[category==sel_category]))
  dev.off()
  
  ####Figure 2e; S6h
  pdf(paste0("subtype_annotation_relapse_",sel_category,".pdf"),height=4,width=4)
  print(plot_surv(subtypes_relapse[category==sel_category],paste0("subtype_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("subtype_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=4)
  print(plot_dots(subtypes_relapse[category==sel_category]))
  dev.off()
}

##############################################
#stratification from copywriteR 
##############################################

copywriter=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/03-CopywriteR/results/CNAprofiles_single_100kb/summary/CNA_Chromosome.tsv"))
setnames(copywriter,"sample","N_number_seq")

copywriter=merge(copywriter,unique(annotation[,c("N_number_seq","patID"),with=FALSE]),by="N_number_seq")
copywriter=copywriter[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery%in%c(1,2)]

#now convert info about lenth of deletions/ amplifications into mere presence (1) absence (0) information
temp_patID=copywriter$patID
temp_surgery=copywriter$surgery
temp_category=copywriter$category

copywriter[,N_number_seq:=NULL,]
copywriter[,date:=NULL,]
copywriter[,surgery:=NULL,]
copywriter[,category:=NULL,]
copywriter[,IDH:=NULL,]

copywriter[copywriter!=0]=1

copywriter[,patID:=temp_patID,]
copywriter[,surgery:=temp_surgery,]
copywriter[,category:=temp_category,]

copywriter[,all_del:=as.integer(rowSums(copywriter[,grep("deletion",names(copywriter)),with=FALSE])),]
copywriter[,all_ampl:=as.integer(rowSums(copywriter[,grep("amplification",names(copywriter)),with=FALSE])),]
copywriter[,all_CNV:=as.integer(rowSums(copywriter[,grep("deletion|amplification",names(copywriter)),with=FALSE])),]

copywriter[,all_del_bin:=binarize(all_del,0.5),by=c("category","surgery")]
copywriter[,all_ampl_bin:=binarize(all_ampl,0.5),by=c("category","surgery")]
copywriter[,all_CNV_bin:=binarize(all_CNV,0.5),by=c("category","surgery")]

copywriter_wide=reshape(copywriter,timevar="surgery",idvar=c("patID","category"),direction="wide")

copywriter_date_surv=merge(copywriter_wide,annot_surv,by=c("patID","category"))
copywriter_date_relapse=merge(copywriter_wide,annot_relapse,by=c("patID","category"))

#chr10q deletion overview
chr10q_del_surv=copywriter_date_surv[,list(N=.N,mean_follow_up=mean(follow_up,na.rm=TRUE),sd_follow_up=sd(follow_up,na.rm=TRUE)),by=c("10_q_deletion.1","10_q_deletion.2","category")]
chr10q_del_rel=copywriter_date_relapse[,list(N=.N,mean_follow_up=mean(follow_up,na.rm=TRUE),sd_follow_up=sd(follow_up,na.rm=TRUE)),by=c("10_q_deletion.1","10_q_deletion.2","category")]
write.table(chr10q_del_surv,"chr10q_del_surv.tsv",quote=FALSE,sep="\t",row.names=FALSE)
write.table(chr10q_del_rel,"chr10q_del_rel.tsv",quote=FALSE,sep="\t",row.names=FALSE)
mat=matrix(c(35, 34, 16, 26),nrow = 2, dimnames = list(first = c("nodel", "del"),second = c("nodel", "del")))
fisher.test(mat,alternative="two.sided")

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  ####Figure S3b
  pdf(paste0("copywriter_date_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(copywriter_date_surv[category==sel_category],paste0("copywriter_date_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("copywriter_date_annotation_surv_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(copywriter_date_surv[category==sel_category]))
  dev.off()
  
  ####Figure S3b
  pdf(paste0("copywriter_date_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(copywriter_date_relapse[category==sel_category],paste0("copywriter_date_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("copywriter_date_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(copywriter_date_relapse[category==sel_category]))
  dev.off()
}

####################################
##stratification from annotation
####################################
annotation[,ShiftPhenotype:=TumorPhenotype,]
annotation[Shape_shift=="stable",ShiftPhenotype:="stable",]
annotation[ShiftPhenotype=="other",ShiftPhenotype:=NA,]
# warning probably due to missing surgery 2 in the validation cohort
annot_pat=annotation[,list(Bevacizumab=any(Bevacizumab!="no"),position_1=unique(position[which(surgery==1)]),border_1=unique(border[which(surgery==1)]),border_2=unique(border[which(surgery==2)]),contrast_1=unique(contrast[which(surgery==1)]),contrast_2=unique(contrast[which(surgery==2)]),age=min(age),sex=unique(sex),StuppComplete=unique(StuppComplete),center=unique(Center),vasc_fibrosis=unique(`Vascular fibrosis`[which(surgery==2)]),pseudopal_necrosis=unique(`Pseudopalisading necrosis`[which(surgery==2)]),rad_necrosis=unique(`Radiation necrosis`[which(surgery==2)]),ShiftPhenotype=ShiftPhenotype[surgery==1],Shape_shift=Shape_shift[surgery==1],Shape_shift_type=Shape_shift_type[surgery==1],TumorPhenotype=TumorPhenotype[surgery==1]),by=patID]

annot_pat[,NoStupp:=ifelse(StuppComplete=="not applicable/other treatment",1,0),]
annot_pat[,StuppComplete:=ifelse(StuppComplete=="yes",1,0),]
annot_pat[,age:=ifelse(age<50,"young","old"),]
annot_pat[,contrast_1:=gsub("mixed solid.*","solid",contrast_1),]
annot_pat[,contrast_1:=gsub("mixed necrotic.*","necrotic",contrast_1),]
annot_pat[,contrast_2:=gsub("mixed solid.*","solid",contrast_2),]
annot_pat[,contrast_2:=gsub("mixed necrotic.*","necrotic",contrast_2),]
annot_pat[,vasc_fibrosis_red:=ifelse(vasc_fibrosis=="Abundant","Present",vasc_fibrosis),]

annot_pat_surv=merge(annot_pat,annot_surv,by="patID")
annot_pat_relapse=merge(annot_pat,annot_relapse,by="patID")

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  ####Figure 4j
  pdf(paste0("general_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_surv[category==sel_category],paste0("general_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("general_annotation_surv_dp_",sel_category,".pdf"),height=4,width=4)
  print(plot_dots(annot_pat_surv[category==sel_category]))
  dev.off()
  
  ####Figure 4j
  pdf(paste0("general_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_relapse[category==sel_category],paste0("general_annotation_relapse_",sel_category,".pdf"),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("general_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=4)
  print(plot_dots(annot_pat_relapse[category==sel_category]))
  dev.off()
}


##############################################
##stratification from histo segmentation
##############################################
segmentation_cols=grep("Mean |Mode |File|Image|Block| STD |mm|background|Area|Necrosis|Whole|Scars",column_annotation_combined_clean$histo_segmentation,invert=TRUE,value=TRUE)

sub=melt(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c(segmentation_cols,"patID","category","surgery"),with=FALSE],id.vars=c("patID","surgery","category"),)
sub[,value.bin:=binarize(value,0.5),by=c("surgery","variable","category")]
sub[,variable:=gsub("|-|/|\\+|\\[|\\|\\]|#|\xb2","",variable),]
sub=unique(sub)

annot_pat=reshape(sub,idvar=c("patID","surgery","category"),timevar="variable",direction="wide",drop="value")
annot_pat=reshape(annot_pat,idvar=c("patID","category"),timevar="surgery",direction="wide")
setnames(annot_pat,names(annot_pat),gsub("value.bin.","",names(annot_pat)))

annot_pat_surv=merge(annot_pat,annot_surv,by=c("patID","category"))
annot_pat_relapse=merge(annot_pat,annot_relapse,by=c("patID","category"))

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  
  pdf(paste0("segmentation_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_surv[category==sel_category],paste0("segmentation_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("segmentation_annotation_surv_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_surv[category==sel_category]))
  dev.off()
  
  pdf(paste0("segmentation_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_relapse[category==sel_category],paste0("segmentation_annotation_relapse_",sel_category,".pdf"),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("segmentation_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_relapse[category==sel_category]))
  dev.off()
}


#######################################
##MGMT priomoter meth stratification###
#######################################

mgmt_cols=column_annotation_combined_clean$mgmt_status
mgmt_cols=mgmt_cols[!mgmt_cols%in%c("mgmt_methyl","meth_max","meth_min" , "mgmt_readCount", "mgmt_CpGcount")]

sub=annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&!is.na(mgmt_conf4),c(mgmt_cols,"patID","surgery"),with=FALSE]
sub=unique(sub)

annot_pat=reshape(sub,idvar="patID",timevar="surgery",direction="wide")

annot_pat_surv=merge(annot_pat,annot_surv,by="patID")
annot_pat_relapse=merge(annot_pat,annot_relapse,by="patID")
annot_pat_relapse_est=merge(annot_pat,annot_relapse_est,by="patID")

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  ####Figure S2d
  pdf(paste0("mgmt_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_surv[category==sel_category],paste0("mgmt_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("mgmt_annotation_surv_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_surv[category==sel_category]))
  dev.off()
  
  ####Figure S2d
  pdf(paste0("mgmt_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_relapse[category==sel_category],paste0("mgmt_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("mgmt_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_relapse[category==sel_category]))
  dev.off()
  
  ####Figure S2d
  pdf(paste0("mgmt_annotation_relapse_est_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_relapse_est[category==sel_category],paste0("mgmt_annotation_relapse_est_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("mgmt_annotation_relapse_est_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_relapse_est[category==sel_category]))
  dev.off()
}

#progression and extension cohort together

pdf(paste0("mgmt_annotation_surv_",paste0(categories,collapse="_"),".pdf"),height=3.5,width=3)
print(plot_surv(annot_pat_surv[category%in%categories],paste0("mgmt_annotation_surv_",paste0(categories,collapse="_")))[[1]])
dev.off()

pdf(paste0("mgmt_annotation_surv_dp_",paste0(categories,collapse="_"),".pdf"),height=4,width=3)
print(plot_dots(annot_pat_surv[category%in%categories]))
dev.off()


pdf(paste0("mgmt_annotation_relapse_est_",paste0(categories,collapse="_"),".pdf"),height=3.5,width=3)
print(plot_surv(annot_pat_relapse_est[category%in%categories],paste0("mgmt_annotation_relapse_est_",paste0(categories,collapse="_")),type="rel")[[1]])
dev.off()

pdf(paste0("mgmt_annotation_relapse_est_dp_",paste0(categories,collapse="_"),".pdf"),height=4,width=3)
print(plot_dots(annot_pat_relapse_est[category%in%categories]))
dev.off()

##############################################
##stratification from annotation immune cells
##############################################
histo_cols=c("CD163","CD3","CD34","CD45ro","CD68","CD8","CD80","FOXP3","MIB","EZH2","Galectin","HLA-DR","PD1","Tim3","cell")

sub=melt(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c(histo_cols,"patID","surgery","category"),with=FALSE],id.vars=c("patID","surgery","category"),)
sub[,value.bin:=binarize(value,0.5),by=c("surgery","variable","category")]

annot_pat=reshape(sub,idvar=c("patID","surgery","category"),timevar="variable",direction="wide",drop="value")
annot_pat=reshape(annot_pat,idvar=c("patID","category"),timevar="surgery",direction="wide")
setnames(annot_pat,names(annot_pat),gsub("value.bin.","",names(annot_pat)))

annot_pat_surv=merge(annot_pat,annot_surv,by=c("patID","category"))
annot_pat_relapse=merge(annot_pat,annot_relapse,by=c("patID","category"))

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  ####Figure 3c; 4b; S8f; S11b
  pdf(paste0("immune_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_surv[category==sel_category],paste0("immune_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("immune_annotation_surv_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_surv[category==sel_category]))
  dev.off()
  
  ####Figure 3c; 4b; S8f; S11b
  pdf(paste0("immune_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_relapse[category==sel_category],paste0("immune_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("immune_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=3)
  print(plot_dots(annot_pat_relapse[category==sel_category]))
  dev.off()
}

##############################################
##stratification from annotation immaging
##############################################
annotation[,progression_types:=ifelse(`T2 diffus`==1,"T2 diffus",ifelse(`classic T1`==1,"classic T1",ifelse(`cT1 flare up`==1,"cT1 flare up",ifelse(`T2 circumscribed`==1,"T2 circumscribed",NA)))),]

img_cols=c("location","position","border","contrast","Numberoflesions","Siteofsurgery","Side","Extentofresection","classic T1","cT1 flare up","T2 diffus","T2 circumscribed","primary nonresponder","progression_types")
img_cols_cont=c("Necrotic  (mm3)","Edema  (mm3)","Cystic (mm3)","Enhancing (mm3)","Total (mm3)","Proportion Necrotic (%)","Proportion Edema (%)","Proportion Cystic (%)","Proportion Enhancing (%)")

sub=melt(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c(img_cols,"patID","surgery","category"),with=FALSE],id.vars=c("patID","surgery","category"),)
setnames(sub,"value","value.bin")

sub_cont=melt(annotation[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt",c(img_cols_cont,"patID","surgery","category"),with=FALSE],id.vars=c("patID","surgery","category"))
sub_cont[,value.bin:=binarize(value,0.5),by=c("surgery","variable","category")]

sub=rbindlist(list(sub,sub_cont),use.names=TRUE,fill=TRUE)
sub[,variable:=gsub("|-|/|\\+|\\[|\\|\\]|#|\\(|\\)|%","",variable),]
sub=unique(sub)

annot_pat=reshape(sub,idvar=c("patID","surgery","category"),timevar="variable",direction="wide",drop="value")
annot_pat=reshape(annot_pat,idvar=c("patID","category"),timevar="surgery",direction="wide")
setnames(annot_pat,names(annot_pat),gsub("value.bin.","",names(annot_pat)))

annot_pat_surv=merge(annot_pat,annot_surv,by=c("patID","category"))
annot_pat_relapse=merge(annot_pat,annot_relapse,by=c("patID","category"))

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){
  ####Figure S9e,g
  pdf(paste0("imaging_annotation_surv_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_surv[category==sel_category],paste0("imaging_annotation_surv_",sel_category))[[1]])
  dev.off()
  
  pdf(paste0("imaging_annotation_surv_dp_",sel_category,".pdf"),height=4,width=4)
  print(plot_dots(annot_pat_surv[category==sel_category]))
  dev.off()
  
  ####Figure S9e,g
  pdf(paste0("imaging_annotation_relapse_",sel_category,".pdf"),height=3.5,width=3)
  print(plot_surv(annot_pat_relapse[category==sel_category],paste0("imaging_annotation_relapse_",sel_category),type="rel")[[1]])
  dev.off()
  
  pdf(paste0("imaging_annotation_relapse_dp_",sel_category,".pdf"),height=4,width=4)
  print(plot_dots(annot_pat_relapse[category==sel_category]))
  dev.off()
}