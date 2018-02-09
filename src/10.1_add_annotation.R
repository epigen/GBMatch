library(project.init)
project.init2("GBMatch")


##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/")
dir.create(out_dir)
setwd(out_dir)

##get annotation
annotation=fread("annotation_combined.tsv")

#get CNV chromosome arms
CNV_chrom_arms=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/03-CopywriteR/results/CNAprofiles_single_100kb/summary/CNA_Chromosome.tsv"))
CNV_chrom_arms=CNV_chrom_arms[,-c("date","category","surgery","IDH"),with=FALSE]
setnames(CNV_chrom_arms,"sample","N_number_seq")


#get meth heterogeneity
loadCaches("combined_heterogeneity")
setnames(combined_heterogeneity,"id","N_number_seq")


#get mutational load
mut_load=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/04-bissnp/bissnp_var_samp.tsv"))
mut_load=mut_load[,-c("patID","surgery","category"),with=FALSE]
setnames(mut_load,"sample_name","N_number_seq")

#get transcriptional subgroups
transc_subgroups=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.1-GBM_classifier/class_probs_annot_27_noNeural_predRRBS.tsv"))
transc_subgroups=transc_subgroups[!is.na(sub_group)]
transc_subgroups[,sub_group_prob:=get(sub_group),by=sample]
transc_subgroups=transc_subgroups[,-c("sampleName","patID","category","surgery.x", "WHO2016_classification", "Follow-up_years","IDH"),with=FALSE]
setnames(transc_subgroups,"sample","N_number_seq")

#get DipScores
dip_scores_1=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/09-dipScore/dipScores_Astrocyte|ESC_sel.tsv"))
setnames(dip_scores_1,"id","N_number_seq")
dip_scores_2=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/09-dipScore/dipScores_Astrocyte|H1_Cell_Line|H9_Cell_Line_sel_seg.tsv"))
setnames(dip_scores_2,"id","N_number_seq")
dip_scores=merge(dip_scores_1,dip_scores_2,by="N_number_seq")


#get MGMT status
mgmt_status=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/02-meth_overview/mgmt_status_all.tsv"))
setnames(mgmt_status,"sample","N_number_seq")
mgmt_status[,surgery.x:=NULL,]
mgmt_status[,patID:=NULL,]
mgmt_status[,category:=NULL,]
mgmt_status[,IDH:=NULL,]

#column annotation
load("column_annoation.RData")
column_annotation_mol=list(CNV_chrom_arms=names(CNV_chrom_arms), meth_heterogeneity=names(combined_heterogeneity),mutational_load=names(mut_load),transc_subgroups=names(transc_subgroups),dipScores=names(dip_scores),mgmt_status=names(mgmt_status))

column_annotation_combined=c(column_annotation,column_annotation_mol)

#clean from N_numbers
column_annotation_combined_clean=lapply(column_annotation_combined,function(x){grep("N_number",x,invert=TRUE,value=TRUE)})


#column categories
technical=c(column_annotation_combined_clean$random_fragmentation,column_annotation_combined_clean$sequencing_stats,column_annotation_combined_clean$sequencing_annot,"Image width [pix]","Image height [pix]","Block size N","material")
medical=c(column_annotation_combined_clean$psa,column_annotation_combined_clean$clinical_annotation,column_annotation_combined_clean$imaging,column_annotation_combined_clean$imaging_progression,"WHO2016_classification", "WHO2016_classification_comment")
measurement=c(column_annotation_combined_clean$histo_classification,column_annotation_combined_clean$histo_segmentation,column_annotation_combined_clean$imaging_segmentation,column_annotation_combined_clean$CNV_chrom_arms,column_annotation_combined_clean$CNV_GOIs,column_annotation_combined_clean$meth_heterogeneity,column_annotation_combined_clean$mutational_load,column_annotation_combined_clean$transc_subgroups,column_annotation_combined_clean$dipScores,column_annotation_combined_clean$mgmt_status)

column_category=list(technical=technical,medical=medical[!medical%in%technical],measurement=measurement[!measurement%in%c(technical,medical)])




#now merge
merge_objects=c("annotation","CNV_chrom_arms","combined_heterogeneity","mut_load","transc_subgroups","dip_scores","mgmt_status")

combined_annotation=data.table()
for(merge_object in merge_objects){
  print(merge_object)
  if(length(combined_annotation)==0){combined_annotation=get(merge_object)}else{
    combined_annotation=merge(combined_annotation,get(merge_object),by="N_number_seq",all=TRUE)
  }
}

combined_annotation[,flowcell.x:=NULL]
setnames(combined_annotation,"flowcell.y","flowcell")
combined_annotation[,surgery.y:=NULL]
combined_annotation[,surgery:=NULL]
setnames(combined_annotation,"surgery.x","surgery")


#test if all columns are contained in column annoations
sort(names(combined_annotation)[!names(combined_annotation)%in%unlist(column_annotation_combined_clean)])

#when none are missing --> save
save(column_category,file="column_category_combined.RData")
save(column_annotation_combined_clean,file="column_annotation_combined.RData")

write.table(combined_annotation,"annotation_combined_final.tsv",sep="\t",row.names=FALSE,quote=FALSE)
