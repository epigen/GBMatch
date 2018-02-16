library(project.init)
project.init2("GBMatch")
library(simpleCache)
library("impute")
library(pheatmap)
library(sva)

source(file.path(getOption("PROJECT.DIR"),"src/99-liblinearFunctions.R"))

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.2-meth_pred")
dir.create(out_dir)
setwd(out_dir)

#get meth data from caches
cacheName="rrbsTiled5ksub"
loadCaches(cacheName)
meth_data_dt_all=get(cacheName)

#load annotation
annotation_all=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
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


#analyze only selected features in primary cohort
anno_sub=prepped_data$annotation[category%in%c("GBMatch","GBmatch_add")]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_selected"),set_targets=c("CD163","CD3","CD68","CD45ro","CD8","EZH2","HLA-DR","CD34","cell","MIB","Mean AVG Eccentricity Tumor","Mean COV Eccentricity Tumor","Mean Nuclei # Tumor","Relative share necrosis","IDH","Sex"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10)


#analyze only selected features in validation cohort
anno_sub=prepped_data$annotation[category%in%c("GBmatch_val")]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("val_selected"),set_targets=c("CD163","CD3","CD68","CD45ro","CD8","EZH2","HLA-DR","CD34","cell","MIB","Mean AVG Eccentricity Tumor","Mean COV Eccentricity Tumor","Mean Nuclei # Tumor","Relative share necrosis","IDH","Sex"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10)


#only surgery 1 in primary
anno_sub=prepped_data$annotation[category%in%c("GBMatch")&surgery.x==1&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_s1_selected"),set_targets=c("StuppComplete","Shape_shift","TumorPhenotype","Follow-up_years","timeToFirstProg","timeToSecSurg","Sex"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10)

#only surgery 2 in primary
anno_sub=prepped_data$annotation[category%in%c("GBMatch")&surgery.x==2&IDH=="wt"]
meth_sub=prepped_data$meth_data_imputed[anno_sub$N_number_seq,]
stopifnot(row.names(meth_sub)==anno_sub$N_number_seq)

cost=heuristicC(meth_sub)
meth_pred_analysis(meth_data_imputed=meth_sub,annotation=anno_sub,column_annotation=column_annotation,to_analyse=c("prim_s2_selected"),set_targets=c("StuppComplete","Shape_shift","TumorPhenotype","Follow-up_years","timeToFirstProg","timeToSecSurg","Sex"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost,cross=10)


##full run on all data and all the different annotations and data in loo-cv --> takes ages (adapt e.g. to run only on primary or validation cohort)
meth_sel="_bval_notscaled_0.9_cross"
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)




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


#set up analysis
#load previously built (on original cohort) and crossvalidated classifiers
meth_sel="bval_notscaled_0.9_cross"
load("dat_prim_selected_bval_notscaled_0.9_cross.RData")


####Matching methylation matrix on which cross validation was performed##############
#like above
meth_data_dt=meth_data_dt_all[id%in%annotation_all[category%in%c("GBMatch","GBmatch_val","GBmatch_add")]$N_number_seq]
set_scale=FALSE
prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9)

#create subsets of samples to analyze
prepped_data_prim=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBMatch"&IDH=="wt"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBMatch"&IDH=="wt"])
stopifnot(row.names(prepped_data_prim$meth_data_imputed)==prepped_data_prim$annotation$N_number_seq)

prepped_data_val=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBmatch_val"&IDH=="wt"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBmatch_val"&IDH=="wt"])
stopifnot(row.names(prepped_data_val$meth_data_imputed)==prepped_data_val$annotation$N_number_seq)


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
factor="Sex"

for (rank_cutoff in seq(150,200,by=2)){
  
  features_res=get_features(pl=pl,prepped_data=prepped_data_prim,factor=factor,rank_cutoff=rank_cutoff)
  bin=c()
  bin[features_res$comp]=c("#ff9289","#00d8e0")
  colors=list(bin=bin,weight=c("#af8dc3","#f7f7f7","#7fbf7b"))
  ph=pheatmap(t(features_res$data),clustering_distance_rows=dist(t(scale(features_res$data,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(features_res$data,center=TRUE,scale=TRUE)),cluster_rows=FALSE,cluster_cols=TRUE,annotation_col=features_res$annot_row,annotation_row=features_res$annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=features_res$colors,main=paste0(factor," ", rank_cutoff))

  print(ph)
}


#set factor and rank_cutoff

factor_list=list(c("MIB",150),c("CD45ro",150),c("CD163",190),c("CD3",196),c("CD68",164),c("CD8",190),c("CD34",196),c("HLA-DR",188),c("Sex",196))


for (factor in factor_list){
  rank_cutoff=as.numeric(factor[2])
  factor=factor[1]
  message(factor)
  
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
  
  roc_mat=multi_roc(class=comp[1],decisionValues=prediction_val_nona[,comp[1],with=FALSE],true_labels=prediction_val_nona$true_label)
  roc_mat[,c("mode","run"):=list("noRand",0),]
  
  for (i in 1:10){
  roc_mat_rand=multi_roc(class=comp[1],decisionValues=prediction_val_nona[,comp[1],with=FALSE],true_labels=sample(prediction_val_nona$true_label))
  roc_mat_rand[,c("mode","run"):=list("rand",i),]
  roc_mat=rbindlist(list(roc_mat,roc_mat_rand))
  }

  auc_annot=roc_mat[,list(mean_auc=mean(auc)),by=c("mode")]
  auc_annot[,x:=0.8,]
  auc_annot[,y:=ifelse(mode=="rand",0.03,0.1),]
  auc_annot[,label:=ifelse(mode=="rand",paste0("Control=",round(as.numeric(mean_auc),2)),paste0("AUC=",round(as.numeric(mean_auc),2))),]
  
pdf(paste0("feature_analysis/val_cross_prediction_ROC_",factor,".pdf"),height=3,width=4)
  print(ggplot(roc_mat,aes(x=fpr,y=tpr,color=mode))+geom_line(aes(group=run,alpha=mode))+geom_text(data=auc_annot,alpha=1,aes(x=x,y=y,label=label))+scale_color_manual(values=c("rand"="darkgrey","noRand"="blue"))+scale_alpha_manual(values=c("rand"=0.4,"noRand"=1))+ggtitle(paste0(factor," ",paste0(comp,collapse=" vs. "))))
 dev.off() 
}else{print("Too few measurements in validation cohort. Skipping cross prediction.")}
######analyze features (heatmap etc)  
#get features

if (val_labels>=5){
data_list=list(all=prepped_data,GBMatch=prepped_data_prim,GBmatch_val=prepped_data_val)
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


#heatmap analysis (clustering)
bin=c()
bin[features_res$comp]=c("#ff9289","#00d8e0")
colors=list(bin=bin,weight=c("#af8dc3","#ffffff","#7fbf7b"))
pdf(paste0("feature_analysis/heatmap_",factor,"_",meth_sel,"_",analysis,".pdf"),height=5,width=9)
pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=FALSE,annotation_col=annot_row,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=factor)

pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=TRUE,annotation_col=annot_row,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=factor)
dev.off()

png(paste0("feature_analysis/heatmap_",factor,"_",meth_sel,"_",analysis,".png"),height=250,width=270)
pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,treeheight_col=25,clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=FALSE,annotation_col=annot_row,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=paste0(factor," top ", rank_cutoff," features"),fontsize=11,legend=FALSE,annotation_legend=FALSE)
dev.off()


#feature distribution
sel_annot=merge(sel,annot_row,by=0)
sel_long=as.data.table(melt(sel_annot,id.vars=c(paste0(gsub("-",".",factor)),"bin","Row.names","category"),value.name="methyl"))
sel_long=merge(sel_long,data.table(weights=annot_col$weight,variable=row.names(annot_col)),by="variable")
sel_long[,weight_group:=ifelse(weights>0,"weights +","weights -"),]

all_annot=merge(sel_data$meth_data_imputed,annot_row,by=0,all.y=TRUE)
all_long=as.data.table(melt(all_annot,id.vars=c(gsub("-",".",factor),"bin","Row.names","category"),value.name="methyl"))
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
}}