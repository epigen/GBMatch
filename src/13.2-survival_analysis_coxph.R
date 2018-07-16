#NOTE: This script performs survival regression analysis using Cox Proportional-Hazards Models.
library(project.init)
project.init2("GBMatch")
library(data.table)
library(ggplot2)
library(survival)
theme_set(theme_bw())

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/13.2-survival_coxph/")
dir.create(out_dir)
setwd(out_dir)

##functons
binarize=function(vector,low,high){
  ret=ifelse(vector<=quantile(vector,low,na.rm=TRUE),"low",ifelse(vector>=quantile(vector,high,na.rm=TRUE),"high",NA))
  return(as.character(ret))
}

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined_final.tsv"))
#column name lists
load(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/column_annotation_combined.RData"))

annotation_sub=annotation[IDH=="wt"&category%in%c("GBMatch","GBmatch_val")&surgery%in%c(1,2)]
annotation_sub[,Age_1:=min(Age),by="patID"]
annotation_sub[,mean_entropy_conf:=ifelse(enrichmentCycles>12&enrichmentCycles<16,mean_entropy,NA),]
annotation_sub[,mean_pdr_conf:=ifelse(enrichmentCycles>12&enrichmentCycles<16,mean_pdr,NA),]
annotation_sub[,DPM_conf:=ifelse(enrichmentCycles>12&enrichmentCycles<16,DPM,NA),]
annotation_sub[,paste0(grep("EPM_",names(annotation_sub),value=TRUE),"_conf"):=mget(grep("EPM_",names(annotation_sub),value=TRUE)),]
annotation_sub[!(enrichmentCycles>12&enrichmentCycles<16),grep("EPM_.*_conf",names(annotation_sub),value=TRUE):=NA,]


#strange thing: reshape drops rows that contain a missing value in an idvar column --> don't use timeToFirstProg in idvar (although it should be there)
annotation_wide=reshape(annotation_sub,idvar=c("patID","category","Follow-up_years","VitalStatus","Sex","Age_1"),timevar="surgery",drop=grep("N_|N-",names(annotation_sub),value=TRUE),direction="wide")
setnames(annotation_wide,c("timeToFirstProg.1","Follow-up_years"),c("PFS","OS"))
annotation_wide[,PFS:=PFS/30.4,]
annotation_wide[,OS:=OS*12,]
annotation_wide[,VitalStatus:=ifelse(VitalStatus=="dead",1,0)]
annotation_wide[,event:=1]
annotation_wide[,Age_bin:=ifelse(Age_1>60,"old","young"),]
annotation_wide[,StuppComplete:=ifelse(StuppComplete.1=="yes","yes","no"),]

#univariate models
covariates=c("Age_1","Sex","mgmt_methyl.1","mgmt_methyl.2","10_q_deletion.1","10_q_deletion.2","7_p_amplification.1","7_p_amplification.2","StuppComplete","CD163.1","CD163.2","CD68.1","CD68.2","MIB.1","MIB.2","mean_entropy_conf.1","mean_entropy_conf.2","mean_pdr_conf.1","mean_pdr_conf.2","DPM_conf.1",paste0(grep("EPM_.*_conf",names(annotation_sub),value=TRUE),".1"),"EZH2_(39875)__NH-A_None_6488.1","EZH2_(39875)__NH-A_None_6488.2","NANOG__Embryonic Stem Cell_NA_18532.1","NANOG__Embryonic Stem Cell_NA_18532.2","mean_diffmeth.Wnt signalling.1","diff_trend_dist_norm.1","Necrotic  (mm3).1","Necrotic  (mm3).2","Enhancing (mm3).1","Enhancing (mm3).2")

sel_categories=c("GBmatch_val","GBMatch") 
modes=c("OS","PFS")
for (sel_category in sel_categories){
  for (mode in modes){  
    
    uni_models=data.table()
    for (covariate in covariates){
      print(covariate)
      data=annotation_wide[category%in%sel_category]
      if(length(na.omit(data[,get(covariate),]))==0){message("No data. Skipping");next}
      if(mode=="OS"){
        x=summary(coxph(Surv(time = get(mode), event = VitalStatus) ~ get(covariate), data=data))
      }
      if(mode=="PFS"){
        x=summary(coxph(Surv(time = get(mode), event = event) ~ get(covariate), data=data))
      }
      
      nevent=x$nevent
      pvalue=signif(x$coef[5],digits=2)
      pvalue_wald=signif(x$wald["pvalue"], digits=2)
      wald.test=signif(x$wald["test"], digits=2)
      beta=signif(x$coef[1], digits=2)
      HR =signif(x$coef[2], digits=2)
      HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
      HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
      HR = paste0(HR, " (", 
                  HR.confint.lower, "-", HR.confint.upper, ")")
      res=data.table(covariate=covariate,events=nevent,beta=beta, `HR (95% CI for HR)`=HR, wald.test=wald.test, pvalue_wald=pvalue_wald,pvalue=pvalue)
      uni_models=rbindlist(list(uni_models,res)) 
    }
    uni_models[,code:=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*",""))),]
    
    write.table(uni_models,paste0(paste0(sel_category,collapse="_"),"_",mode,"_uni_models.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    
    #multivariate model with significant covariates
    multi_formula=as.formula(paste0("Surv(time = ",mode,", event = ",ifelse(mode=="PFS","event","VitalStatus"),") ~ ",paste0("`",uni_models[pvalue<0.05]$covariate,"`",collapse="+")))
    multi_cox=summary(coxph(multi_formula , data=data))
    
    capture.output(multi_cox,file=paste0(paste0(sel_category,collapse="_"),"_",mode,"_multi_model_signif.txt"))
    
    #check out multivariate model for DNA meth heterogeneity measures
    cv_methhet_list=list(cv_1=c("mean_entropy_conf.1","mean_pdr_conf.1","EZH2_(39875)__NH-A_None_6488.1","NANOG__Embryonic Stem Cell_NA_18532.1"),cv_2=c("mean_entropy_conf.2","mean_pdr_conf.2","EZH2_(39875)__NH-A_None_6488.2","NANOG__Embryonic Stem Cell_NA_18532.2"))
    
    for (cv_methhet in names(cv_methhet_list)){
      if(sel_category=="GBmatch_val"&cv_methhet=="cv_2"){next}
      
      multi_formula_methhet=as.formula(paste0("Surv(time = ",mode,", event = ",ifelse(mode=="PFS","event","VitalStatus"),") ~ ",paste0("`",cv_methhet_list[[cv_methhet]],"`",collapse="+")))
      multi_cox_methhet=summary(coxph(multi_formula_methhet , data=data))
      
      capture.output(multi_cox_methhet,file=paste0(paste0(sel_category,collapse="_"),"_",mode,"_multi_model_methhet_",cv_methhet,".txt"))
    }
  }
}

