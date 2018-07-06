library(project.init)
project.init2("GBMatch")
library(simpleCache)
library("impute")
library(pheatmap)
library(sva)
source(file.path(getOption("PROJECT.DIR"),"src/99-liblinearFunctions.R"))

##functions##########
buildJadd=function(cols,funcs,special){
  r = paste("list(", paste(c(paste0(cols, "=", funcs, "(", cols, ")"),special), collapse=","), ")")
  return(r);
}

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.2-meth_pred/validation")
dir.create(out_dir)
setwd(out_dir)

#create cache including "plate2" and rerun samples (only needed first time)
# eload(loadAnnotation())
# eload(loadUCSCRepeats())
# repeats = SV$repeats
# eload(loadAnnotation())
# BSSamples = SV$psa[file.exists(BSFile) & library=="RRBS", sample_name]; BSSamples
# eload(loadTiles(genomeBuild=genome, tileSize=5000))
# agRegions = get(paste0("tiles5000", genome), env=SV)
# 
# cols=c("hitCount", "readCount", "methyl")
# funcs = c("sum", "sum", "mean")
# special=c("CpGcount=length(na.omit(methyl))")
# jCommand = buildJadd(cols,funcs,special)
# 
# setLapplyAlias(12)
# 
# 
# # 5kb tiles after subtracting repeats
# simpleCache("rrbsTiled5ksub_val", {
#   sampleSummaryList = lapplyAlias(BSSamples, sampleSummaryByRegion,
#                                   regions=agRegions, excludeGR = repeats,
#                                   SV$psa, idColumn = "sample_name", fileColumn="BSFile", 
#                                   cachePrepend="meth.5k.", cacheSubDir="meth/tile5k_sub_val", 
#                                   jCommand=jCommand, 
#                                   readFunction=function(x) {
#                                     ino = BSreadBiSeq(x);
#                                     ino[,methyl:=hitCount/readCount]
#                                   }, mc.preschedule=FALSE)
#   sampleSummaryLong = rbindlist(sampleSummaryList)
#   sampleSummaryLong # Cache this.
# }, recreate=TRUE, noload=TRUE)


#get meth data from caches
cacheName="rrbsTiled5ksub_val"
loadCaches(cacheName)
meth_data_dt_all=get(cacheName)

#load annotation
annotation_all=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined_includeAll.tsv"))
#add rerun samples
rerun_stats=fread(file.path(getOption("PROCESSED.PROJECT"),"results_pipeline/RRBS_stats_summary_rerun.tsv"))
rerun_stats=rerun_stats[grep("rerun",sampleName)]
rerun_stats[,Sex:=annotation_all[N_number_seq==gsub("_rerun","",sampleName)]$Sex,by=1:nrow(rerun_stats)]
rerun_stats[,position_dna:=annotation_all[N_number_seq==gsub("_rerun","",sampleName)]$position_dna,by=1:nrow(rerun_stats)]
setnames(rerun_stats,"sampleName","N_number_seq")
annotation_all=rbindlist(list(annotation_all,rerun_stats),use.names=TRUE,fill=TRUE)
annotation_all[grep("rerun",N_number_seq),c("category","plate_dna_abbrev","IDH"):=list("rerun","rerun","wt"),]

#column name lists
load(file.path(getOption("PROCESSED.PROJECT"),"results_analysis_rev_prelim/01.1-combined_annotation/column_annoation.RData"))
column_annotation=lapply(column_annotation,function(x){grep("N_number|N-number|ID|read_length",x,value=TRUE,invert=TRUE)})

meth_data_dt_all[,region:=paste0(chr,"_",start,"_",end),]

#select samples to use
#beta value
meth_data_dt=meth_data_dt_all[id%in%annotation_all[category%in%c("GBMatch","GBmatch_val","rerun")&IDH=="wt"]$N_number_seq]
prepped_data=prepare_data(meth_data_dt=meth_data_dt,annotation_all=annotation_all,min_na_ratio=0.9,min_na_ratio_samp=0.45)


#create subsets of samples to analyze
prepped_data_prim=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBMatch"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBMatch"])
stopifnot(row.names(prepped_data_prim$meth_data_imputed)==prepped_data_prim$annotation$N_number_seq)

prepped_data_val=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category%in%c("GBmatch_val","rerun")]$N_number_seq,],annotation=prepped_data$annotation[category%in%c("GBmatch_val","rerun")])
stopifnot(row.names(prepped_data_val$meth_data_imputed)==prepped_data_val$annotation$N_number_seq)


prepped_data_val_plate1=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBmatch_val"&plate_dna_abbrev=="EA003_plate1"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBmatch_val"&plate_dna_abbrev=="EA003_plate1"])
stopifnot(row.names(prepped_data_val_plate1$meth_data_imputed)==prepped_data_val_plate1$annotation$N_number_seq)

prepped_data_val_plate2=list(meth_data_imputed=prepped_data$meth_data_imputed[prepped_data$annotation[category=="GBmatch_val"&plate_dna_abbrev=="EA003_plate2"]$N_number_seq,],annotation=prepped_data$annotation[category=="GBmatch_val"&plate_dna_abbrev=="EA003_plate2"])
stopifnot(row.names(prepped_data_val_plate2$meth_data_imputed)==prepped_data_val_plate2$annotation$N_number_seq)



#estimate cost for primary cohort, because this is the one the perdictor is buit on
cost=heuristicC(prepped_data_prim$meth_data_imputed)


#build classifier on primary dataset
classification_prim=check_prediction(data=prepped_data_prim$meth_data_imputed,labels=prepped_data_prim$annotation$Sex,samples=prepped_data_prim$annotation$N_number_seq,cross=10,nReps=3,type=4,cost=cost)

#plot ROC
classification_prim$plot


#use classifier to predict sex for validation cohort

p_val=predict(classification_prim$model,prepped_data_val$meth_data_imputed,proba=FALSE,decisionValues=TRUE)

prediction_val=data.table(sample=prepped_data_val$annotation$N_number_seq,predicted_label=p_val$predictions,true_label=prepped_data_val$annotation$Sex,p_val$decisionValues,prepped_data_val$annotation[,c("plate_dna_abbrev","position_dna","position_rrbs","experiment","pool","adapter"),with=FALSE])

roc_mat=multi_roc(class="m",decisionValues=prediction_val[,"m"],true_labels=prediction_val$true_label)

ggplot(roc_mat,aes(x=fpr,y=tpr))+geom_line(color="blue")+annotate(geom="text",alpha=1,x=0.8,y=0.1,label=paste0("auc=",round(as.numeric(unique(roc_mat$auc)),2)))+annotate(geom="text",x=0,y=1,hjust=0,vjust=1,label="m")

#plot plate positions
prediction_val_long_pre=melt(prediction_val,id.vars=c("sample","predicted_label","true_label","m","f","plate_dna_abbrev","experiment","pool","adapter"),variable.name="position_type",value.name="position")
prediction_val_long=melt(prediction_val_long_pre,id.vars=c("sample","predicted_label","true_label","m","f","position_type", "position","pool","adapter"),variable.name="plate_type",value.name="plate")
prediction_val_long[,plate_type:=ifelse(plate_type=="plate_dna_abbrev","plate_dna",ifelse(plate_type=="experiment","plate_rrbs",NA)),]
prediction_val_long[,position.x:=as.numeric(sub("[A-H]","",position)),]
prediction_val_long[,position.y:=factor(gsub("[0-9]","",position),levels=c("H","G","F","E","D","C","B","A")),]
prediction_val_long[,concordance:=predicted_label==true_label,]

prediction_val_long2=melt(prediction_val_long,id.vars=c("sample","m","f","position_type", "position", "plate_type","plate", "position.x", "position.y","concordance","pool","adapter"),variable.name="label_type",value.name="label")
prediction_val_long2[,plate:=factor(plate,levels=c("EA003_plate1","AK55","EA003_plate2","rerun","AK53","EA003_plate3","AK54")),]

prediction_val_long2=prediction_val_long2[(position_type=="position_rrbs"&plate_type=="plate_rrbs")|(position_type=="position_dna"&plate_type=="plate_dna")]


pdf("validation_cohort_sex_concordance_dna_rrbs.pdf",height=12,width=6)
ggplot(prediction_val_long2[!is.na(plate)],aes(x=factor(position.x),y=position.y,fill=label,col=concordance))+geom_tile(width=0.5,height=0.5,lwd=1)+facet_wrap(plate~label_type,ncol=2)+scale_color_manual(values=c("TRUE"="grey","FALSE"="black"))+xlab("position.x")
dev.off()

write.table( prediction_val[predicted_label!=true_label],"validation_cohort_sex_misclassified.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#only plot DNA plate
pdf("validation_cohort_sex_concordance_dna.pdf",height=7,width=6)
ggplot(prediction_val_long2[plate_type=="plate_dna"&!is.na(plate)],aes(x=factor(position.x),y=position.y,fill=label,col=concordance))+geom_tile(width=0.5,height=0.5,lwd=1)+facet_wrap(plate~label_type,ncol=2)+scale_color_manual(values=c("TRUE"="grey","FALSE"="black"))+xlab("position.x")
dev.off()

#color pool
ggplot(prediction_val_long2,aes(x=factor(position.x),y=position.y,fill=factor(pool),col=concordance))+geom_tile(width=0.5,height=0.5,lwd=1)+facet_wrap(plate~label_type,ncol=2)+scale_color_manual(values=c("TRUE"="grey","FALSE"="black"))+xlab("position.x")+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a", "#33a02c",  "#fb9a99", "#e31a1c",  "#fdbf6f",  "#ff7f00",  "#cab2d6",  "#6a3d9a","#ffff99"))

#color adapter
ggplot(prediction_val_long2[!is.na(plate)],aes(x=factor(position.x),y=position.y,fill=factor(adapter),col=concordance))+geom_tile(width=0.5,height=0.5,lwd=1)+facet_wrap(plate~label_type,ncol=2)+scale_color_manual(values=c("TRUE"="grey","FALSE"="black"))+xlab("position.x")

#write mw sequence for plate 2
predicted_mf_col=paste0(prediction_val_long2[plate=="EA003_plate2"&label_type=="predicted_label"][order(position.x,as.character(position.y))]$label,collapse="")
gregexpr("fmmmmff",predicted_mf_col)

#create random distribution concordance/discordance count
sub_pred=prediction_val[plate_dna_abbrev=="EA003_plate2",c("sample", "predicted_label", "true_label"),with=FALSE]
for (i in 1:1000){
sub_pred[,paste0("rand",i):=sample(true_label,length(true_label),replace=FALSE),]
}
sub_pred_long=melt(sub_pred,id.vars=c("sample", "predicted_label", "true_label"))

sub_pred_long_red=sub_pred_long[,list(concordant=sum(predicted_label==value),discordant=sum(predicted_label!=value),m=sum(value=="m"),f=sum(value=="f")),by="variable"]

sub_pred_long_red_lon=melt(sub_pred_long_red,id.vars="variable",variable.name="count_type")
ggplot(sub_pred_long_red_lon,aes(x=value))+geom_histogram()+facet_wrap(~count_type,scale="free")


#complementary analysis using read coverage on x and y chromosomes

par=data.table(chr=c("chrX","chrX","chrY","chrY"),start=c(10001,155701383,10001,56887903),end=c(2781479,156030895,2781479,57217415),par=c("parX1","parX2","parY1","parY2"))
setkey(par, chr, start, end)
rrbsData_chrX_Y=foverlaps(meth_data_dt_all[chr%in%c("chrX","chrY")],par,type="within",mult="first")
rrbsData_chrX_Y=rrbsData_chrX_Y[is.na(par)]
rrbsData_chrX_Y[,totalReadCount:=sum(readCount),by="id"]
rrbsData_chrX_Y[,maxDetected:=length(unique(regionID)),by="chr"]

rrbsData_normCov=rrbsData_chrX_Y[is.na(par),list(detected=.N/maxDetected[1],normReadCount=sum(readCount)/totalReadCount[1]),by=c("id","chr")]
XY_detected_ratio=rrbsData_normCov[,list(ratio=detected[which(chr=="chrY")]/detected[which(chr=="chrX")]),by=c("id")]

annot_plot=with(annotation_all,data.table(IDH=IDH,Sex=Sex,category=category,plate_dna=plate_dna_abbrev,id=N_number_seq))
annot_plot=merge(annot_plot,setnames(prediction_val[,c("sample","predicted_label"),],"sample","id"),by="id",all=TRUE)


XY_detected_ratio_annot=merge(XY_detected_ratio,annot_plot,by="id")
XY_detected_ratio_annot[,id:=factor(id,levels=unique(id[order(ratio)])),]
XY_detected_ratio_annot[plate_dna=="",plate_dna:=NA,]
XY_detected_ratio_annot[is.na(plate_dna),plate_dna:="other"]


ggplot(XY_detected_ratio_annot[category%in%c("GBmatch_val","GBMatch")],aes(x=id,y=ratio,fill=Sex),)+geom_bar(stat="identity")+geom_point(y=-0.1,aes(col=predicted_label))+theme(axis.text.x  = element_text(angle=90, vjust=0.5,hjust=1))+ylab("YX-ratio")+ylim(c(-0.2,2))+scale_color_manual(values=c("red","blue"),na.value="transparent")

ggplot(XY_detected_ratio_annot[category%in%c("GBmatch_val","GBMatch")],aes(x=id,fill=Sex,group=Sex))+geom_density(alpha=0.5)+geom_point(y=-0.0005,aes(col=predicted_label))+theme(axis.text.x  = element_text(angle=90, vjust=0.5,hjust=1))+scale_color_manual(values=c("red","blue"),na.value="transparent")+ylim(c(-0.001,NA))

pdf("XYratio_density.pdf",height=4,width=7)
ggplot(XY_detected_ratio_annot[category%in%c("GBmatch_val","GBMatch")&plate_dna!="EA003_plate3"],aes(x=log2(ratio),fill=plate_dna,group=plate_dna))+geom_density(alpha=0.5,bw=0.1)+geom_point(y=-0.0005,shape="|",aes(col=Sex),alpha=0.6,size=3)+theme(axis.text.x  = element_text(angle=90, vjust=0.5,hjust=1))+scale_color_manual(values=c("red","blue"),na.value="transparent")+ylim(c(-0.001,NA))+xlab("log2(detectedY/detectedX)")
dev.off()

#########further prediction tests on subsets
#check classification on validation cohort without misclassified samples from primary prediction
#complete val cohort
cost_val=heuristicC(prepped_data_val$meth_data_imputed)
classification_val_complete=check_prediction(data=prepped_data_val$meth_data_imputed,labels=prepped_data_val$annotation$Sex,samples=prepped_data_val$annotation$N_number_seq,cross=10,nReps=3,type=4,cost=cost_val)
#plot ROC
classification_val_complete$plot

#reduced val cohort
correct_samples=prediction_val[predicted_label==true_label]$sample

cost_val_red=heuristicC(prepped_data_val$meth_data_imputed[correct_samples,])
classification_val_red=check_prediction(data=prepped_data_val$meth_data_imputed[correct_samples,],labels=prepped_data_val$annotation[N_number_seq%in%correct_samples,]$Sex,samples=prepped_data_val$annotation[N_number_seq%in%correct_samples,]$N_number_seq,cross=10,nReps=3,type=4,cost=cost_val_red)
#plot ROC
classification_val_red$plot


# separated plate 1 and 2

cost_val_p1=heuristicC(prepped_data_val_plate1$meth_data_imputed)
classification_val_p1=check_prediction(data=prepped_data_val_plate1$meth_data_imputed,labels=prepped_data_val_plate1$annotation$Sex,samples=prepped_data_val_plate1$annotation$N_number_seq,cross=10,nReps=3,type=4,cost=cost_val_p1)
#plot ROC
classification_val_p1$plot

cost_val_p2=heuristicC(prepped_data_val_plate2$meth_data_imputed)
classification_val_p2=check_prediction(data=prepped_data_val_plate2$meth_data_imputed,labels=prepped_data_val_plate2$annotation$Sex,samples=prepped_data_val_plate2$annotation$N_number_seq,cross=10,nReps=3,type=4,cost=cost_val_p2)
#plot ROC
classification_val_p2$plot



#now try reduced validation data set for immune prediction
meth_pred_analysis(meth_data_imputed=prepped_data_val$meth_data_imputed[correct_samples,],annotation=prepped_data_val$annotation[N_number_seq%in%correct_samples,],column_annotation=column_annotation,to_analyse=c("immuno_sex"),set_targets=c("MIB","CD68","CD163","Sex"),meth_sel="_val_red_",set_scale=FALSE,type=4,cost=cost_val_red,cross=10,nReps=3)

meth_pred_analysis(meth_data_imputed=prepped_data_val$meth_data_imputed,annotation=prepped_data_val$annotation,column_annotation=column_annotation,to_analyse=c("immuno_sex"),set_targets=c("MIB","CD68","CD163","Sex"),meth_sel="_val_complete_",set_scale=FALSE,type=4,cost=cost_val,cross=10,nReps=3)

meth_pred_analysis(meth_data_imputed=prepped_data_val_plate1$meth_data_imputed,annotation=prepped_data_val_plate1$annotation,column_annotation=column_annotation,to_analyse=c("immuno_sex"),set_targets=c("MIB","CD68","CD163","Sex"),meth_sel="_val_plate1",set_scale=FALSE,type=4,cost=cost_val_p1,cross=10,nReps=10)

meth_pred_analysis(meth_data_imputed=prepped_data_val_plate2$meth_data_imputed,annotation=prepped_data_val_plate2$annotation,column_annotation=column_annotation,to_analyse=c("immuno_sex"),set_targets=c("MIB","CD68","CD163","Sex"),meth_sel="_val_plate2",set_scale=FALSE,type=4,cost=cost_val_p1,cross=10,nReps=10)



#predict immune infiltration in validation cohort through prim model

meth_pred_analysis(meth_data_imputed=prepped_data_prim$meth_data_imputed,annotation=prepped_data_prim$annotation,column_annotation=column_annotation,to_analyse=c("histo_immuno"),set_targets=c("MIB","CD68","CD163"),meth_sel="_immuno_prim_",set_scale=FALSE,type=4,cost=cost,cross=10,nReps=3)
load("dat_histo_immuno_immuno_prim_.RData")

factor="MIB"
factor="CD68"
factor="CD163"

#reduced val 
p_val=predict(pl[[factor]]$model,prepped_data_val$meth_data_imputed[correct_samples,],proba=FALSE,decisionValues=TRUE)
prediction_val=data.table(sample=prepped_data_val$annotation[N_number_seq%in%correct_samples,]$N_number_seq,predicted_label=p_val$predictions,true_label=binarize(prepped_data_val$annotation[N_number_seq%in%correct_samples,factor,with=FALSE],low=0.2,high=0.8),p_val$decisionValues)

#complete val 
p_val=predict(pl[[factor]]$model,prepped_data_val$meth_data_imputed,proba=FALSE,decisionValues=TRUE)
prediction_val=data.table(sample=prepped_data_val$annotation$N_number_seq,predicted_label=p_val$predictions,true_label=binarize(prepped_data_val$annotation[,factor,with=FALSE],low=0.2,high=0.8),p_val$decisionValues)

#plate1
p_val=predict(pl[[factor]]$model,prepped_data_val_plate1$meth_data_imputed,proba=FALSE,decisionValues=TRUE)
prediction_val=data.table(sample=prepped_data_val_plate1$annotation$N_number_seq,predicted_label=p_val$predictions,true_label=binarize(prepped_data_val_plate1$annotation[,factor,with=FALSE],low=0.2,high=0.8),p_val$decisionValues)

#plate2
p_val=predict(pl[[factor]]$model,prepped_data_val_plate2$meth_data_imputed,proba=FALSE,decisionValues=TRUE)
prediction_val=data.table(sample=prepped_data_val_plate2$annotation$N_number_seq,predicted_label=p_val$predictions,true_label=binarize(prepped_data_val_plate2$annotation[,factor,with=FALSE],low=0.2,high=0.8),p_val$decisionValues)


#plot here
roc_mat=multi_roc(class="high",decisionValues=prediction_val[!is.na(true_label),"high"],true_labels=prediction_val[!is.na(true_label)]$true_label)
ggplot(roc_mat,aes(x=fpr,y=tpr))+geom_line(color="blue")+annotate(geom="text",alpha=1,x=0.8,y=0.1,label=paste0("auc=",round(as.numeric(unique(roc_mat$auc)),2)))+annotate(geom="text",x=0,y=1,hjust=0,vjust=1,label="high")

