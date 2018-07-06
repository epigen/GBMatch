library(project.init)
project.init2("GBMatch")
library(simpleCache)
library("impute")
library(pheatmap)
library(sva)
library(LOLA)

source(file.path(getOption("PROJECT.DIR"),"src/99-liblinearFunctions.R"))

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.2-meth_pred")
dir.create(out_dir)
setwd(out_dir)

#get meth data from caches
cacheName="rrbsTiled5ksub"
loadCaches(cacheName)
meth_data_dt_all=get(cacheName)
rm(rrbsTiled5ksub)

#load annotation
annotation_all=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
annotation_all[,cohort_surv_median:=median(`Follow-up_years`,na.rm=TRUE),by=c("category","IDH","surgery.x")]
annotation_all[,survival_class:=ifelse(category=="GBMatch"&`Follow-up_years`>=cohort_surv_median,"long",ifelse(category=="GBmatch_val"&`Follow-up_years`<=cohort_surv_median,"short",NA)),]
annotation_all[,timeToFirstProg_est:=ifelse(is.na(timeToFirstProg),`Follow-up_years`*12-2,timeToFirstProg),]


#column name lists
load(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/column_annoation.RData"))
column_annotation=lapply(column_annotation,function(x){grep("N_number|N-number|ID|read_length",x,value=TRUE,invert=TRUE)})

meth_data_dt_all[,region:=paste0(chr,"_",start,"_",end),]

#select samples to use
#prepare all data together and split after (makes cross prediction possible)
meth_data_dt=meth_data_dt_all[id%in%annotation_all[category%in%c("GBMatch","GBmatch_val","GBmatch_add")]$N_number_seq]
set_scale=FALSE
prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9)
meth_sel="_bval_notscaled_0.9_cross"
#remove unneeded objects: needed when running in parallel (cores>1)
rm(meth_data_dt_all)
rm(meth_data_dt)
gc()

#analyze only selected features in primary cohort
####Figure S10b
anno_sub=prepped_data$annotation[category%in%c("GBMatch","GBmatch_add")]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_selected"),set_targets=c("CD163","CD3","CD68","CD45ro","CD8","EZH2","HLA-DR","CD34","cell","MIB","Mean AVG Eccentricity Tumor","Mean COV Eccentricity Tumor","Mean Nuclei # Tumor","Relative share necrosis","IDH","Sex","age"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10,nReps=100)
#comment: age is readily predictable here because this data set also includes the IDH mut/ oligodendroglioma samples, the patients of which are younger on average. So actually the classifier finds the IDH signature and not actually a signature for age.

#analyze only selected features in GBMatch IDH=wt
####Figure 3g; 4c,i; S8g; S10b,d

anno_sub=prepped_data$annotation[category%in%c("GBMatch")&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_wt_selected"),set_targets=c("CD163","CD3","CD68","CD45ro","CD8","EZH2","HLA-DR","CD34","cell","MIB","Mean AVG Eccentricity Tumor","Mean COV Eccentricity Tumor","Mean Nuclei # Tumor","Relative share necrosis","Sex","age"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10,nReps=100,cores=32)
#comment: set cores > 1 to seed up. But don't do this on login node!!


#analyze only selected features in validation cohort
anno_sub=prepped_data$annotation[category%in%c("GBmatch_val")&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("val_wt_selected"),set_targets=c("CD163","CD68","MIB","Mean AVG Eccentricity Tumor","Mean COV Eccentricity Tumor","Mean Nuclei # Tumor","Relative share necrosis","StuppComplete","Follow-up_years","timeToFirstProg","Sex","age"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10,nReps=100)


#only surgery 1 in primary
anno_sub=prepped_data$annotation[category%in%c("GBMatch")&surgery.x==1&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_s1_selected"),set_targets=c("StuppComplete","Shape_shift","TumorPhenotype","Follow-up_years","timeToFirstProg","timeToSecSurg","Sex","age"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10,nReps=100)

#only surgery 2 in primary
anno_sub=prepped_data$annotation[category%in%c("GBMatch")&surgery.x==2&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_s2_selected"),set_targets=c("StuppComplete","Shape_shift","TumorPhenotype","Follow-up_years","timeToFirstProg","timeToSecSurg","Sex","age"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10,nReps=100)

#selected in progression and validation cohort only surgery 1, wt
#Figure 6g
anno_sub=prepped_data$annotation[category%in%c("GBMatch","GBmatch_val")&surgery.x==1&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim&val_s1_selected"),set_targets=c("StuppComplete","survival_class","Follow-up_years","timeToFirstProg_est","Sex","age"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10,nReps=100,cores=32)
#note core >1 --> don't run on login node!! (srun -n1 --cpus-per-task=32 --mem-per-cpu=8000 --part=longq --time=4-00:00:00 R --no-save)


##full run on all data and all the different annotations and data in loo-cv --> takes ages (adapt e.g. to run only on primary or validation cohort)
#meth_sel="_bval_notscaled_0.9_cross"
#meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)


#############################folowup on genomic regions (features) selected in the final models################################
###functions#####
get_features=function(pl,prepped_data,factor,rank_cutoff=50){ 
  weights_df=as.data.frame(t(pl[[factor]]$model$W))
  weights_df$region=row.names(weights_df)
  weights_dt=as.data.table(weights_df)
  
  annotation_sel=prepped_data$annotation[!is.na(get(factor))]
  if (is.numeric(unlist(annotation_sel[,factor,with=FALSE]))){
    annotation_sel[,bin:=binarize(get(factor),0.2,0.8),by="category"]
  }else{
    annotation_sel[,bin:=get(factor),]
  }
  comp=sort(annotation_sel[,unique(bin)])
  #annotation_sel[,bin:=binarize(get(factor),0.5,0.5)]
  
  weights_dt[,rank:=rank(-abs(get(comp[1])),ties.method="random"),]
  
  #now continue the analysis
  #sel_weights=weights_dt[rank<=rank_cutoff]
  sel_weights=weights_dt[rank<=rank_cutoff&region%in%colnames(prepped_data$meth_data_imputed)]
  sel=prepped_data$meth_data_imputed[row.names(prepped_data$meth_data_imputed)%in%annotation_sel[!is.na(bin)]$N_number_seq,sel_weights$region]
  
  annot_row=data.frame(cont=annotation_sel[!is.na(bin),factor,with=FALSE],bin=annotation_sel[!is.na(bin)]$bin,category=annotation_sel[!is.na(bin)]$category,row.names=annotation_sel[!is.na(bin)]$N_number_seq)
  
  #order samples by annotation
  sel=sel[rownames(annot_row[order(annot_row[,3],annot_row[,1]),]),]
  sel=sel[,sel_weights[order(get(comp[1]))]$region]
  annot_col=data.frame(weight=weights_dt[region%in%colnames(sel),get(comp[1]),],row.names=weights_dt[region%in%colnames(sel)]$region)
  
  return(list(annot_col=annot_col,annot_row=annot_row,data=sel,comp=comp))
}

plot_heatmap=function(file_name=NULL,pl,prepped_data,factor,rank_cutoff,title=factor){
  features_res=get_features(pl=pl,prepped_data=prepped_data,factor=factor,rank_cutoff=rank_cutoff)
  bin=c()
  bin[features_res$comp]=c("#ff9289","#00d8e0")
  colors=list(bin=bin,weight=c("#af8dc3","#f7f7f7","#7fbf7b"),category=c("GBMatch"="grey","GBmatch_val"="black"))
  if(!is.null(file_name)){
    pdf(paste0(file_name,".pdf"),height=5,width=9)
  }
  pheatmap(t(features_res$data),clustering_distance_rows=dist(t(scale(features_res$data,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(features_res$data,center=TRUE,scale=TRUE)),cluster_rows=FALSE,cluster_cols=TRUE,show_rownames=FALSE,show_colnames=FALSE,annotation_col=features_res$annot_row,annotation_row=features_res$annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=paste0(title," ", rank_cutoff))
  if(!is.null(file_name)){
    dev.off()
  }
  if(!is.null(file_name)){
    if(all(features_res$annot_row[,factor]==features_res$annot_row[,"bin"])){features_res$annot_row[,factor]=NULL}
    png(paste0(file_name,".png"),height=250,width=270)
  pheatmap(t(features_res$data),clustering_distance_rows=dist(t(scale(features_res$data,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(features_res$data,center=TRUE,scale=TRUE)),annotation_names_row=FALSE,annotation_names_col=FALSE,treeheight_col=25,cluster_rows=FALSE,cluster_cols=TRUE,show_rownames=FALSE,show_colnames=FALSE,annotation_col=features_res$annot_row,annotation_row=features_res$annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=paste0(title," top ", rank_cutoff," features"),fontsize=11,legend=FALSE,annotation_legend=FALSE)
    dev.off()
  #write heatmap data
  write.table(as.data.table(t(features_res$data),keep.rownames="pos"),paste0(file_name,"_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
  }
}


#set up analysis
#load previously built (on original cohort) and crossvalidated classifiers
meth_sel="bval_notscaled_0.9_cross"
load("dat_prim_wt_selected_bval_notscaled_0.9_cross.RData")


####Matching methylation matrix on which cross validation was performed##############
#like above (only rerun when not available)
#meth_data_dt=meth_data_dt_all[id%in%annotation_all[category%in%c("GBMatch","GBmatch_val","GBmatch_add")]$N_number_seq]
#set_scale=FALSE
#prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9)

#create subsets of samples to analyze
prepped_data_prim=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBMatch"&IDH=="wt"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBMatch"&IDH=="wt"])
stopifnot(row.names(prepped_data_prim$meth_data_imputed)==prepped_data_prim$annotation$N_number_seq)

prepped_data_val=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBmatch_val"&IDH=="wt"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBmatch_val"&IDH=="wt"])
stopifnot(row.names(prepped_data_val$meth_data_imputed)==prepped_data_val$annotation$N_number_seq)

prepped_data_wt=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[IDH=="wt"]$N_number_seq,],annotation=prepped_data$annotation[IDH=="wt"])
stopifnot(row.names(prepped_data_wt$meth_data_imputed)==prepped_data_wt$annotation$N_number_seq)


#make results directory
dir.create("feature_analysis")

#run this to try different rank_cutoffs
factor="MIB"
factor="CD45ro"
factor="CD163"
factor="CD3"
factor="CD68"
factor="CD8"
factor="CD34"
factor="HLA-DR"
factor="EZH2"
factor="Sex"
factor="Relative share necrosis"
factor="Mean AVG Eccentricity Tumor" #196
factor="Mean COV Eccentricity Tumor" #160
#only for survival prediction in prim&val
factor="survival_class" #rank_cutoff=200

for (rank_cutoff in seq(150,200,by=2)){
plot_heatmap(pl=pl,prepped_data=prepped_data_prim,factor=factor,rank_cutoff=rank_cutoff)
}

#set factor and rank_cutoff for complrehensive analysis
#complete val
#factor_list=list(c("MIB",150),c("CD45ro",150),c("CD163",190),c("CD3",196),c("CD68",164),c("CD8",190),c("CD34",196),c("HLA-DR",188),c("Sex",196))
#without val plate 2 including IDH mut
#factor_list=list(c("MIB",190),c("CD45ro",184),c("CD163",182),c("CD3",184),c("CD68",174),c("CD8",174),c("CD34",150),c("HLA-DR",184),c("Sex",188))
#without val plate 2 without IDH mut
factor_list=list(c("MIB",194),c("CD45ro",196),c("CD163",160),c("CD3",158),c("CD68",178),c("CD8",178),c("CD34",164),c("HLA-DR",194),c("EZH2",190),c("Sex",186),c("Relative share necrosis",164),c("Mean AVG Eccentricity Tumor" ,196),c("Mean COV Eccentricity Tumor" ,160))


####Figure 4d,e; S8g; S10c,e,f,g; S11c-e

for (factor in factor_list){
  rank_cutoff=as.numeric(factor[2])
  factor=factor[1]
  message(factor)
  
  ##write data used for roc curve
  #select the level that was chosen for the roc plots in the paper
  roc_level=sort(names(pl[[factor]]$plot_int))
  roc_level=roc_level[roc_level%in%c("high","m")]
  roc_data=pl[[factor]]$plot_int[[roc_level]]$data
  write.table(roc_data,paste0("feature_analysis/roc_data_",factor,".csv"),sep=";",quote=FALSE,row.names=FALSE)
  
  val_labels=sum(prepped_data_val$annotation[,!is.na(get(factor)),])
  #######use predictor built on primary dataset to predict validation dataset  
  if (val_labels>=5){
    p_val=predict(pl[[factor]]$model,prepped_data_val$meth_data_imputed,proba=FALSE,decisionValues=TRUE)
    
    if (is.numeric(unlist(prepped_data_val$annotation[,factor,with=FALSE]))){
      true_label=binarize(prepped_data_val$annotation[,factor,with=FALSE],low=0.2,high=0.8)
    }else{
      true_label=unlist(prepped_data_val$annotation[,factor,with=FALSE])
    }
    
    prediction_val=data.table(sample=prepped_data_val$annotation$N_number_seq,predicted_label=p_val$predictions,true_label=true_label,prepped_data_val$annotation[,factor,with=FALSE],p_val$decisionValues)
    prediction_val_nona=prediction_val[!is.na(true_label)]
    comp=sort(prediction_val_nona[,unique(true_label)])
    
    sample_N=unique(prediction_val_nona[,list(sample,true_label),])[,.N,by=true_label]
    sample_N[,lab:=paste0("N(",true_label,")=",N),]
    
    roc_mat=multi_roc(class=comp[1],decisionValues=prediction_val_nona[,comp[1],with=FALSE],true_labels=prediction_val_nona$true_label)
    roc_mat[,c("mode","run"):=list("noRand",0),]
    
    for (i in 1:100){
      roc_mat_rand=multi_roc(class=comp[1],decisionValues=prediction_val_nona[,comp[1],with=FALSE],true_labels=sample(prediction_val_nona$true_label))
      roc_mat_rand[,c("mode","run"):=list("rand",i),]
      roc_mat=rbindlist(list(roc_mat,roc_mat_rand))
    }
    
    #interpolate and merge roc curves
    roc_mat_int=roc_mat[,approx(fpr,tpr,xout=seq(0,1,0.01),yleft=0,yright=1,method="constant",ties=max),by=c("auc","run","mode")]
    setnames(roc_mat_int,c("x","y"),c("fpr","tpr"))
    roc_mat_int[fpr==0,tpr:=0,]
    roc_mat_int[fpr==1,tpr:=1,]
    
    auc_annot=roc_mat[,list(mean_auc=mean(auc)),by=c("mode")]
    auc_annot[,x:=0.8,]
    auc_annot[,y:=ifelse(mode=="rand",0.03,0.1),]
    auc_annot[,label:=ifelse(mode=="rand",paste0("Control=",round(as.numeric(mean_auc),2)),paste0("AUC=",round(as.numeric(mean_auc),2))),]
    
    pdf(paste0("feature_analysis/val_cross_prediction_ROC_",factor,".pdf"),height=5,width=4.5)
    print(ggplot(roc_mat,aes(x=fpr,y=tpr,color=mode))+geom_line(aes(group=run,alpha=mode))+geom_text(data=auc_annot,alpha=1,aes(x=x,y=y,label=label))+scale_color_manual(values=c("rand"="darkgrey","noRand"="blue"))+scale_alpha_manual(values=c("rand"=0.4,"noRand"=1))+ylab("TPR")+xlab("FPR")+ggtitle(paste0(c(paste0(factor," ",paste0(comp,collapse=" vs. ")),sample_N$lab),collapse="\n"))+coord_fixed(ratio=1))
    print(ggplot(roc_mat_int,aes(x=fpr,y=tpr,group=mode))+stat_summary(geom="ribbon", fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), fill="lightgrey",alpha=0.4)+stat_summary(geom="line",aes(col=mode), fun.y=mean)+geom_text(data=auc_annot,alpha=1,aes(x=x,y=y,label=label))+scale_color_manual(values=c("rand"="darkgrey","noRand"="blue"))+ylab("TPR")+xlab("FPR")+ggtitle(paste0(c(paste0(factor," ",paste0(comp,collapse=" vs. ")),sample_N$lab),collapse="\n"))+coord_fixed(ratio=1)) 
    
    dev.off() 
  }else{print("Too few measurements in validation cohort. Skipping cross prediction.")}
  ######analyze features (heatmap etc)  
  #get features
  
  if (val_labels>=5){
    data_list=list(all=prepped_data_wt,GBMatch=prepped_data_prim,GBmatch_val=prepped_data_val)
  }else{
    data_list=list(GBMatch=prepped_data_prim)
    print("Too few measurements in validation cohort. Feature analysis only on primary cohort.")
  }
  
  
  for (i in 1:length(data_list)){
    sel_data=data_list[[i]]
    analysis=names(data_list[i])
    message(analysis)
    
    features_res=get_features(pl=pl,prepped_data=sel_data,factor=factor,rank_cutoff=rank_cutoff)
    
    sel=features_res$data
    annot_row=features_res$annot_row
    annot_col=features_res$annot_col
    
    annot_row_hm=annot_row
    if(all(annot_row_hm[,gsub("-| ",".",factor)]==annot_row_hm[,"bin"])){annot_row_hm[,gsub("-| ",".",factor)]=NULL}
    
    
    #heatmap analysis (clustering)
    bin=c()
    bin[features_res$comp]=c("#ff9289","#00d8e0")
    
    colors=list(bin=bin,weight=c("#af8dc3","#ffffff","#7fbf7b"),category=c("GBMatch"="grey","GBmatch_val"="black"))
    pdf(paste0("feature_analysis/heatmap_",factor,"_",meth_sel,"_",analysis,".pdf"),height=5,width=9)
    pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=FALSE,annotation_col=annot_row_hm,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=factor)
    
    pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=TRUE,annotation_col=annot_row_hm,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=factor)
    dev.off()
    
    png(paste0("feature_analysis/heatmap_",factor,"_",meth_sel,"_",analysis,".png"),height=250,width=270)
    pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,treeheight_col=25,clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=FALSE,annotation_col=annot_row_hm,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=paste0(factor," top ", rank_cutoff," features"),fontsize=11,legend=FALSE,annotation_legend=FALSE)
    dev.off()
    #write heatmap data
    write.table(as.data.table(t(sel),keep.rownames="pos"),paste0("feature_analysis/heatmap_",factor,"_",meth_sel,"_",analysis,"_data.csv"),sep=";",quote=FALSE,row.names=FALSE)
    
    #feature distribution
    sel_annot=merge(sel,annot_row,by=0)
    sel_long=as.data.table(melt(sel_annot,id.vars=c(paste0(gsub("-| ",".",factor)),"bin","Row.names","category"),value.name="methyl"))
    sel_long=merge(sel_long,data.table(weights=annot_col$weight,variable=row.names(annot_col)),by="variable")
    sel_long[,weight_group:=ifelse(weights>0,"weights +","weights -"),]
    
    all_annot=merge(sel_data$meth_data_imputed,annot_row,by=0,all.y=TRUE)
    all_long=as.data.table(melt(all_annot,id.vars=c(gsub("-| ",".",factor),"bin","Row.names","category"),value.name="methyl"))
    all_long[,weight_group:="core regions"]
    
    sel_all_long=rbindlist(list(sel_long,all_long),use.names=TRUE,fill=TRUE)
    
    spl=unlist(strsplit(sel_long$variable,"_"))
    sel_long[,chr:=spl[seq(1,length(spl),3)],]
    sel_long[,start:=spl[seq(2,length(spl),3)],]
    sel_long[,end:=spl[seq(3,length(spl),3)],]
    
    write.table(sel_long,paste0("feature_analysis/features_weight_meth_",factor,"_",meth_sel,"_",analysis,".tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    
    
    write.table(unique(sel_long[weight_group=="weights +",c("chr","start","end"),with=FALSE]),paste0("feature_analysis/weight_pos_",factor,"_",meth_sel,"_",analysis,".bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(unique(sel_long[weight_group=="weights -",c("chr","start","end"),with=FALSE]),paste0("feature_analysis/weight_neg_",factor,"_",meth_sel,"_",analysis,".bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    
    pdf(paste0("feature_analysis/density_",factor,"_",meth_sel,"_",analysis,".pdf"),height=2.5,width=7.5)
    print(ggplot(sel_all_long,aes(x=methyl,col=bin))+ggtitle(factor)+geom_density()+facet_wrap(~weight_group,scales="free_y"))
    dev.off()
    
    pdf(paste0("feature_analysis/qq_",factor,"_",meth_sel,"_",analysis,".pdf"),height=3,width=7.5)
    print(ggplot(sel_all_long) +stat_qq(aes(sample = methyl,col=bin),geom="line")+ggtitle(factor)+xlim(c(-4,4))+ coord_fixed(ratio=8)+facet_wrap(~weight_group)) 
    dev.off()
  }
}


#########for survival prediction
#load previously built (on combined cohort) and crossvalidated classifiers (for survival presiction)
#Figure 6g
load("dat_prim&val_s1_selected_bval_notscaled_0.9_cross.RData")

prepped_data_prim_val=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category%in%c("GBMatch","GBmatch_val")&surgery.x==1&IDH=="wt"]$N_number_seq,],annotation=prepped_data$annotation[category%in%c("GBMatch","GBmatch_val")&surgery.x==1&IDH=="wt"])
stopifnot(row.names(prepped_data_prim_val$meth_data_imputed)==prepped_data_prim_val$annotation$N_number_seq)
#plot and save only heatmap for selected factor and rank_cutoff
factor="survival_class" #rank_cutoff=200
plot_heatmap(file_name=paste0("feature_analysis/only_heatmap_",factor),pl=pl,prepped_data=prepped_data_prim_val,factor="survival_class",rank_cutoff=200,title="Survival class")

roc_data=pl[[factor]]$plot_int[["long"]]$data
write.table(roc_data,paste0("feature_analysis/roc_data_",factor,".csv"),sep=";",quote=FALSE,row.names=FALSE)


##########feature follow-up with LOLA
#load LOLA DBs
regionDB = loadRegionDB(file.path(Sys.getenv("RESOURCES"),"regions/LOLACore/hg38/"),collections=c("encode_tfbs","codex"))
regionDB_seg = loadRegionDB(file.path(Sys.getenv("RESOURCES"),"regions/customRegionDB/hg38/"),collections=list("roadmap_segmentation"))
cellType_anno=fread(file.path(getOption("PROJECT.DIR"),"metadata/LOLA_annot/CellTypes.tsv"))
cellType_anno_seg=fread(file.path(Sys.getenv("RESOURCES"),"regions/customRegionDB/hg38/roadmap_segmentation/index.txt"),select=c("filename","EID", "seg_code","seg_explanation","Epigenome name (from EDACC Release 9 directory)","ANATOMY"))
setnames(cellType_anno_seg,"Epigenome name (from EDACC Release 9 directory)","cellType_corr")

#get data
load("dat_prim&val_s1_selected_bval_notscaled_0.9_cross.RData")
sel_analysis="survival_class"
sel_level="long"

#convert features to granges
weight_mat=as.data.table(t(pl[[sel_analysis]]$model$W),keep.rownames=TRUE)
weight_mat=weight_mat[rn!="Bias"]
spl=unlist(strsplit(weight_mat$rn,"_"))
weight_mat[,chr:=spl[seq(1,length(spl),3)],]
weight_mat[,start:=as.numeric(spl[seq(2,length(spl),3)]),]
weight_mat[,end:=as.numeric(spl[seq(3,length(spl),3)]),]
weight_mat=weight_mat[order(get(sel_level))]
weight_mat[,weight_group:=ifelse(get(sel_level)>0,"pos","neg"),]

univ=with(weight_mat,GRanges(seqnames = Rle(chr), IRanges(start=start, end=end),strand=Rle("*"),weight=get(sel_level)))
uset=GRangesList()
uset[[paste0(sel_analysis,"_",sel_level,"_neg")]]=univ[1:1000]
uset[[paste0(sel_analysis,"_",sel_level,"_pos")]]=univ[(length(univ)-999):length(univ)]


locResults=runLOLA(userSets=uset,userUniverse=univ,regionDB=regionDB_seg)

locResults[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
locResults[,mlog10p.adjust:=-log10(p.adjust),]
locResults=merge(locResults,cellType_anno_seg,by="filename")

locResults_sel=locResults[p.adjust<0.01&ANATOMY%in%c("ESC","BRAIN")]
locResults_sel[,ANATOMY_cor:=ifelse(cellType_corr=="NH-A_Astrocytes","ASTRO",ANATOMY),]

####Figure 6h
pdf(paste0("feature_analysis/LOLA_",sel_analysis,".pdf"),height=3.5,width=6)
ggplot(locResults_sel,aes(y=-log10(p.adjust),x=seg_explanation,size=logOddsRatio,fill=ANATOMY_cor))+geom_point(alpha=0.6,shape=21,position=position_jitter(width=0.05))+geom_hline(yintercept=-log10(0.01),linetype="dashed",col="grey")+facet_wrap(~userSet)+coord_flip()+xlab("")+ theme(legend.position="bottom",legend.box = "vertical")
dev.off()

write.table(locResults_sel[,list(p.adjust,seg_explanation,logOddsRatio,ANATOMY_cor,userSet),],"feature_analysis/Source Data Figure6h.csv",sep=";",quote=FALSE,row.names=FALSE)




