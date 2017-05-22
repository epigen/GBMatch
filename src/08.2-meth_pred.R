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
#beta value
meth_data_dt=meth_data_dt_all[id%in%annotation_all[category=="GBMatch"&IDH=="wt"]$N_number_seq]
#only for clinical (IDH and Sex discrimination)
#meth_data_dt=meth_data_dt_all[id%in%annotation_all[category%in%c("GBMatch","GBmatch_add")]$N_number_seq]

set_scale=FALSE
prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9)
cost=heuristicC(prepped_data$meth_data_imputed)


##############test stuff########################
meth_sel="_bval_notscaled_0.999"
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=1)

meth_sel="_bval_notscaled_0.999_varCost"
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=c(0.01,1,100))


################use this############################
meth_sel=paste0("_bval_notscaled_0.9_",round(cost,3))
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)

meth_sel=paste0("_bval_notscaled_IDH_sel_0.9_",round(cost,3))
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("clinical_annotation"),set_targets=c("IDH","Sex"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)

meth_sel=paste0("_bval_notscaled_sel_0.9_",round(cost,3))
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_segmentation"),set_targets=c("Median AVG Eccentricity Tumor","Median COV Eccentricity Tumor"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)



#m value
meth_data_dt=meth_data_dt_all[id%in%annotation_all[category=="GBMatch"&IDH=="wt"]$N_number_seq]
meth_data_dt[,methyl:=log2((methyl+0.001)/(1-methyl+0.001)),]
set_scale=TRUE
prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9)
cost=heuristicC(prepped_data$meth_data_imputed)
##############test stuff################################
meth_sel="_mval_scaled_0.999"
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=1)

meth_sel="_mval_scaled_0.999_0.1"
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification")[1],meth_sel=meth_sel,set_scale=set_scale,type=4,cost=0.1)

meth_sel=paste0("_mval_scaled_0.99_",round(cost,3))
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno","histo_segmentation","imaging_segmentation","clinical_annotation","histo_classification"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)

##############use this########################################
meth_sel=paste0("_mval_scaled_sel_0.9_",round(cost,3))
meth_pred_analysis(meth_data_imputed=prepped_data$meth_data_imputed,annotation=prepped_data$annotation,column_annotation=column_annotation,to_analyse=c("histo_segmentation"),set_targets=c("Median AVG Eccentricity Tumor","Median COV Eccentricity Tumor","Mean AVG Eccentricity Tumor","Mean COV Eccentricity Tumor"),meth_sel=meth_sel,set_scale=set_scale,type=4,cost=cost)



##################folowup on regions################################
###functions
get_features=function(pl,prepped_data,factor,rank_cutoff=50){ 
  weights_df=as.data.frame(t(pl[[factor]]$model$W))
  weights_df$region=row.names(weights_df)
  weights_dt=as.data.table(weights_df)
  
  annotation_sel=prepped_data$annotation[!is.na(get(factor))]
  annotation_sel[,bin:=binarize(get(factor),0.2,0.8)]
  
  weights_dt[,rank:=rank(-abs(high),ties.method="random"),]
  
  #now continue the analysis
  sel_weights=weights_dt[rank<=rank_cutoff]
  sel=prepped_data$meth_data_imputed[row.names(prepped_data$meth_data_imputed)%in%annotation_sel[!is.na(bin)]$N_number_seq,sel_weights$region]
  annot_row=data.frame(cont=annotation_sel[!is.na(bin),factor,with=FALSE],bin=annotation_sel[!is.na(bin)]$bin,row.names=annotation_sel[!is.na(bin)]$N_number_seq)
  
  #order samples by annotation
  sel=sel[rownames(annot_row[order(annot_row[,1]),]),]
  sel=sel[,sel_weights[order(high)]$region]
  annot_col=data.frame(weight=weights_dt[region%in%colnames(sel)]$high,row.names=weights_dt[region%in%colnames(sel)]$region)
  
  return(list(annot_col=annot_col,annot_row=annot_row,data=sel))
}

#set up analysis
meth_sel="bval_notscaled_0.9_0.02"
load(paste0("dat_histo_immuno_",meth_sel,".RData"))
####Dont forget to load the matching methylation matrix in the beginning of the script##############
meth_data_dt=meth_data_dt_all[id%in%annotation_all[category=="GBMatch"&IDH=="wt"]$N_number_seq]
set_scale=FALSE
#choose correct min_na_ratio (according to meth_sel)
prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9)
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


for (rank_cutoff in seq(100,200,by=2)){

  features_res=get_features(pl=pl,prepped_data=prepped_data,factor=factor,rank_cutoff=rank_cutoff)

  colors=list(bin=c("high"="#ff9289","low"="#00d8e0"),weight=c("#af8dc3","#f7f7f7","#7fbf7b"))
  ph=pheatmap(t(features_res$data),clustering_distance_rows=dist(t(scale(features_res$data,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(features_res$data,center=TRUE,scale=TRUE)),cluster_rows=TRUE,annotation_col=features_res$annot_row,annotation_row=features_res$annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=features_res$colors,main=paste0(factor," ", rank_cutoff))

  print(ph)
}


#set factor and rank_cutoff
factor="MIB" 
rank_cutoff=186

factor="CD45ro"
rank_cutoff=160

factor="CD163"
rank_cutoff=152

factor="CD3"
rank_cutoff=116

factor="CD68"
rank_cutoff=132

factor="CD8"
rank_cutoff=150

factor="CD34"
rank_cutoff=162

factor="HLA-DR"
rank_cutoff=182


#get features
features_res=get_features(pl=pl,prepped_data=prepped_data,factor=factor,rank_cutoff=rank_cutoff)

sel=features_res$data
annot_row=features_res$annot_row
annot_col=features_res$annot_col


#heatmap analysis (clustering)
colors=list(bin=c("high"="#ff9289","low"="#00d8e0"),weight=c("#af8dc3","#ffffff","#7fbf7b"))
pdf(paste0("feature_analysis/heatmap_",factor,"_",meth_sel,".pdf"),height=5,width=9)
pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=FALSE,annotation_col=annot_row,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=factor)

pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=TRUE,annotation_col=annot_row,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=factor)
dev.off()

png(paste0("feature_analysis/heatmap_",factor,"_",meth_sel,".png"),height=250,width=270)
pheatmap(t(sel),clustering_distance_rows=dist(t(scale(sel,center=TRUE,scale=TRUE))),show_rownames=FALSE,show_colnames=FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,treeheight_col=25,clustering_distance_cols=dist(scale(sel,center=TRUE,scale=TRUE)),cluster_rows=FALSE,annotation_col=annot_row,annotation_row=annot_col,color=colorRampPalette(c("blue" ,"red"))(20),border_color=NA,annotation_colors=colors,main=paste0(factor," top ", rank_cutoff," features"),fontsize=11,legend=FALSE,annotation_legend=FALSE)
dev.off()

#feature distribution
sel_annot=merge(sel,annot_row,by=0)
sel_long=as.data.table(melt(sel_annot,id.vars=c(paste0(gsub("-",".",factor)),"bin","Row.names"),value.name="methyl"))
sel_long=merge(sel_long,data.table(weights=annot_col$weight,variable=row.names(annot_col)),by="variable")
sel_long[,weight_group:=ifelse(weights>0,"weights +","weights -"),]

all_annot=merge(prepped_data$meth_data_imputed,annot_row,by=0,all.y=TRUE)
all_long=as.data.table(melt(all_annot,id.vars=c(gsub("-",".",factor),"bin","Row.names"),value.name="methyl"))
all_long[,weight_group:="core regions"]

sel_all_long=rbindlist(list(sel_long,all_long),use.names=TRUE,fill=TRUE)

spl=unlist(strsplit(sel_long$variable,"_"))
sel_long[,chr:=spl[seq(1,length(spl),3)],]
sel_long[,start:=spl[seq(2,length(spl),3)],]
sel_long[,end:=spl[seq(3,length(spl),3)],]

write.table(unique(sel_long[weight_group=="weights +",c("chr","start","end"),with=FALSE]),paste0("feature_analysis/weight_pos_",factor,"_",meth_sel,".bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(unique(sel_long[weight_group=="weights -",c("chr","start","end"),with=FALSE]),paste0("feature_analysis/weight_neg_",factor,"_",meth_sel,".bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


pdf(paste0("feature_analysis/density_",factor,"_",meth_sel,".pdf"),height=2.5,width=7.5)
ggplot(sel_all_long,aes(x=methyl,col=bin))+ggtitle(factor)+geom_density()+facet_wrap(~weight_group,scales="free_y")
dev.off()

pdf(paste0("feature_analysis/qq_",factor,"_",meth_sel,".pdf"),height=3,width=7.5)
ggplot(sel_all_long) +stat_qq(aes(sample = methyl,col=bin),geom="line")+ggtitle(factor)+xlim(c(-4,4))+ coord_fixed(ratio=8)+facet_wrap(~weight_group) 
dev.off()