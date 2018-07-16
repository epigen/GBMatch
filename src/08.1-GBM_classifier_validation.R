#NOTE: This script validates the RRBS based transcriptional subtype classificaton by comparing to classification based on RNA-seq data
library(project.init)
library(simpleCache)
library(pheatmap)
library(ggtern)
library(pROC)
project.init2("GBMatch")

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.1-GBM_classifier")
dir.create(out_dir)
setwd(out_dir)

#set paths and load required files to assemble RNA data
genome="hg38_nc"
results_pipeline="../../results_pipeline_rna/"
sampleAnnotationFile=fread(file.path(getOption("PROJECT.DIR"),"metadata/GBMatch_samples_rna.csv"))
resultsStats=fread(paste0(results_pipeline,"/rnaBitSeq_stats_summary.tsv"))
transcriptAnnotation=fread(tail(system(paste0("ls /data/groups/lab_bock/shared/resources/genomes/",genome,"/transcripts_*"),intern=TRUE),n=1))
setnames(transcriptAnnotation,names(transcriptAnnotation),c("ensG","ensT","transc_start","transc_end","gene","biotype","strand","chr","ensP"))


#collect RNA data function
collectBitSeqData = function(results_pipeline,sampleAnnotationFile,transcriptAnnotation,genome,resultsStats){
  rna_all=data.table()
  sample_no=0
  for (i in 1:nrow(sampleAnnotationFile)){
    print(sampleAnnotationFile[i]$sample_name)
    if (sampleAnnotationFile[i]$run == 1){
      tr_file=paste0(results_pipeline,"/",sampleAnnotationFile[i]$sample_name,"/bowtie1_",genome,"/bitSeq/",sampleAnnotationFile[i]$sample_name,".tr")
      means_file=paste0(results_pipeline,"/",sampleAnnotationFile[i]$sample_name,"/bowtie1_",genome,"/bitSeq/",sampleAnnotationFile[i]$sample_name,".mean")
      counts_file=paste0(results_pipeline,"/",sampleAnnotationFile[i]$sample_name,"/bowtie1_",genome,"/bitSeq/",sampleAnnotationFile[i]$sample_name,".counts")
      if (file.exists(tr_file)&file.exists(means_file)&file.exists(counts_file)){
        sample_no=sample_no+1
        tr=fread(tr_file)
        setnames(tr,names(tr),c("ensG","ensT","length","length_adj"))
        tr[,ensG:=gsub("\\..*","",ensG),]
        tr[,ensT:=gsub("\\..*","",ensT),]        
        tr[,sample_name:=sampleAnnotationFile[i]$sample_name,]
        tr[,sampleID:=paste0("sample_",sample_no),]
        tr[,Raw_reads:=resultsStats[sampleName==sampleAnnotationFile[i]$sample_name,"Raw_reads",with=FALSE],]
        tr[,Aligned_reads:=resultsStats[sampleName==sampleAnnotationFile[i]$sample_name,"Aligned_reads",with=FALSE],]
        if ("Filtered_reads" %in% names(resultsStats)){
          tr[,Filtered_reads:=resultsStats[sampleName==sampleAnnotationFile[i]$sample_name,"Filtered_reads",with=FALSE],]  
        } 
        means=fread(means_file)
        setnames(means,names(means),c("RPKM","RPKM_var"))
        counts=fread(counts_file)
        combined=cbind(tr,means,counts)
        combined[,ensT:=gsub("\\.[0-9]*","",ensT),]
        merge = merge(transcriptAnnotation,combined,by=c("ensT"),all=TRUE)
        rna_all=rbindlist(list(rna_all,merge))
      }
      else{print("file not found")}
    }
    else{print("not selected")}
  }
  return(rna_all)
}

#actually collect and cache the data
simpleCache(recreate=FALSE,"RNA_data",cacheSubDir="RNA",assignToVariable="RNA_data",instruction="collectBitSeqData(results_pipeline,sampleAnnotationFile,transcriptAnnotation,genome,resultsStats)")
RNA_data_uniq=RNA_data[!is.na(sample_name),list(RPKM=max(RPKM)),by=c("sample_name","gene")]

#check out TCGA samples and data
#get classification info from Verhaak 2009
TCGA_genes=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/TCGA_unified_CORE_ClaNC840_transcSubtype_genes_Verhaak2009.txt"))
all_TCGA_genes=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/unifiedScaled_Verhaak2009.txt"))
all_TCGA_genes_long=melt(all_TCGA_genes,id.vars="gene",variable.name="TCGA_sample")

TCGA_samples_annot=as.data.table(t(TCGA_genes[1,-c(1,2),with=FALSE]),keep.rownames="TCGA_sample")
setnames(TCGA_samples_annot,"V1","sample_subtype")
TCGA_genes=TCGA_genes[-1]
setnames(TCGA_genes,c("CLID","Genes highly expressed in each subtype used for Gene Ontology"),c("gene","gene_subtype"))
TCGA_genes[,gene_subtype:=c("","Classical","Mesenchymal","Neural","Proneural")[as.factor(gene_subtype)],]
TCGA_genes_annot=TCGA_genes[,c(1,2),with=FALSE,]

TCGA_genes_long=melt(TCGA_genes,id.vars=c("gene","gene_subtype"),variable.name="TCGA_sample")
TCGA_genes_long[,value:=as.numeric(value),]
TCGA_genes_long=merge(TCGA_genes_long,TCGA_samples_annot,by="TCGA_sample")
TCGA_genes_wide=reshape(TCGA_genes_long[,c("TCGA_sample","gene","value"),],timevar="TCGA_sample",idvar="gene",direction="wide")
setnames(TCGA_genes_wide,names(TCGA_genes_wide),sub("value.","",names(TCGA_genes_wide)))

ggplot(TCGA_genes_long[gene%in%unique(gene)[1:10]],aes(x=value,color=gene))+geom_density()

#plot TCGA heatmap
TCGA_genes_mat=as.matrix(TCGA_genes_wide[,-c(1),with=FALSE,])
row.names(TCGA_genes_mat)=TCGA_genes_wide$gene
col_annot=data.frame(subtype=TCGA_samples_annot$sample_subtype,row.names=TCGA_samples_annot$TCGA_sample)
row_annot=data.frame(subtype=TCGA_genes_annot$gene_subtype,row.names=TCGA_genes_annot$gene)
pheatmap(TCGA_genes_mat, annotation_col=col_annot,annotation_row=row_annot)

#mean expression of each gene in each subgroup
subtype_expression=TCGA_genes_long[,list(mean_expr=mean(value),sd_expr=sd(value)),by=c("sample_subtype","gene","gene_subtype")]
subtype_expression[,gene_subtype_ass:=sample_subtype[which.max(mean_expr)],by="gene"]
concordant_genes=unique(subtype_expression[gene_subtype==gene_subtype_ass]$gene)
pheatmap(TCGA_genes_mat[row.names(TCGA_genes_mat)%in%concordant_genes,], annotation_col=col_annot,annotation_row=row_annot[row.names(row_annot)%in%concordant_genes,,drop=FALSE])

RNA_data[gene%in%concordant_genes]

#some of the genes (51) have different names in the TCGA data set --> match synonyms
missing_genes=concordant_genes[!concordant_genes%in%unique(RNA_data_uniq$gene)]
NCBI_gene_names=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/Homo_sapiens.gene_info_NCBI"),select=c("GeneID","Symbol","Synonyms"))

synonyms_dt=data.table(missing_gene=missing_genes)
synonyms_dt[,c("GeneID","Symbol","Synonyms"):=NCBI_gene_names[grepl(paste0("\\|",missing_gene,"\\|"),paste0("|",Symbol,"|",Synonyms,"|"))],by=1:nrow(synonyms_dt)]

#check genes with multiple hits
synonyms_dt[25]
NCBI_gene_names[grepl(paste0("\\|EDG1\\|"),paste0("|",Symbol,"|",Synonyms,"|"))]
synonyms_dt[49]
NCBI_gene_names[grepl(paste0("\\|ECGF1\\|"),paste0("|",Symbol,"|",Synonyms,"|"))]
#exclude these two genes, because they are not unambiguously assignable
synonyms_dt=synonyms_dt[-c(25,49)]

#exclude genes that don't have matches
synonyms_dt=synonyms_dt[!is.na(GeneID)]

#match synonyms to RNA table
synonyms_dt[,matching_gene:=paste0(c(Symbol, unlist(strsplit(Synonyms,"\\|")))[c(Symbol, unlist(strsplit(Synonyms,"\\|")))%in%RNA_data_uniq$gene],collapse="|"),by=1:nrow(synonyms_dt)]
#exclude genes with multiple matching genes
synonyms_dt=synonyms_dt[!grepl("\\|",matching_gene)]

#exchange matching genes in RNA_data_uniq to match the TCGA data
#first check if matching genes are already there (should be empty)
TCGA_genes[gene%in%synonyms_dt$matching_gene]
RNA_data_uniq=merge(RNA_data_uniq,setnames(synonyms_dt[,c("missing_gene","matching_gene"),with=FALSE],"matching_gene","gene"),all.x=TRUE)
RNA_data_uniq[!is.na(missing_gene),gene:=missing_gene,]

#Now cluster data according to selected and annotated genes
#discard neural genes
concordant_genes_noNeur=unique(subtype_expression[gene_subtype==gene_subtype_ass&gene_subtype!="Neural"&gene%in%RNA_data_uniq$gene]$gene)
RNA_data_uniq[,RPKM_scaled:=scale(log(RPKM),center=TRUE,scale=TRUE),by="gene"]

RNA_data_uniq_wide=reshape(RNA_data_uniq,idvar="gene",timevar="sample_name",drop=c("RPKM","missing_gene"),direction="wide")

RNA_data_uniq_mat=as.matrix(RNA_data_uniq_wide[gene%in%concordant_genes_noNeur,-c("gene"),with=FALSE])
row.names(RNA_data_uniq_mat)=RNA_data_uniq_wide[gene%in%concordant_genes_noNeur]$gene

sel_row_annot=row_annot[row.names(row_annot)%in%row.names(RNA_data_uniq_mat),,drop=FALSE]

RNA_data_uniq_mat=RNA_data_uniq_mat[row.names(sel_row_annot),]

pdf("subtype_validation_samples_heatmap.pdf",height=10,width=14)
pheatmap(RNA_data_uniq_mat,annotation_row=sel_row_annot,cluster_rows=FALSE)
dev.off()
#TCGA data common genes available in both datasets, only neural
pdf("subtype_validation_TCGA_heatmap.pdf",height=10,width=14)
pheatmap(TCGA_genes_mat[row.names(TCGA_genes_mat)%in%concordant_genes_noNeur,!colnames(TCGA_genes_mat)%in%row.names(col_annot[col_annot$subtype=="Neural",,drop=FALSE])],annotation_col=col_annot[col_annot$subtype!="Neural",,drop=FALSE],annotation_row=row_annot[row.names(row_annot)%in%concordant_genes_noNeur,,drop=FALSE])
dev.off()

##Now perform real classification based on mix profiles

#1. create pure expression profiles for each subtype
subtype_profiles=TCGA_genes_long[gene%in%concordant_genes_noNeur&sample_subtype!="Neural",list(mean_expr=mean(value)),by=c("gene","sample_subtype")]
subtype_profiles_wide=reshape(subtype_profiles,idvar="gene",timevar="sample_subtype",direction="wide")


#2. create all mix-profiles in 1% steps
mix_template=data.table()
for (mes in 0:100){
  for (cla in 0:(100-mes)){
    pro = 100-(mes+cla)
    mixID=paste(mes,cla,pro,sep="_")
    print(mixID)
    mix_template=rbindlist(list(mix_template,data.table(mixID=mixID,mes=mes,cla=cla,pro=pro))) 
  }    
}  

generate_mix=function(data,weights){  
  mix=data[,list(mix=weighted.mean(c(mean_expr.Mesenchymal,mean_expr.Classical,mean_expr.Proneural),w=weights)),by="gene"]
  return(mix)
}

mix_profiles=mix_template[,cbind(mix_profile=generate_mix(subtype_profiles_wide,c(mes,cla,pro))),by=c("mixID","mes","cla","pro")]
setnames(mix_profiles,"mix_profile.gene","gene")

#3. calculate correlation of all validation samples with all mix-profiles
mix_profiles_comb=merge(mix_profiles,RNA_data_uniq,by="gene",allow.cartesian=TRUE)
mix_profiles_cor=mix_profiles_comb[,list(cor=cor(mix_profile.mix,RPKM_scaled)),by=c("sample_name","mixID", "mes", "cla", "pro")]
mix_profiles_cor[,majority_sybtype:=c("Mesenchymal","Classical","Proneural")[which.max(c(mes,cla,pro))],by=c("sample_name","mixID")]
mix_profiles_cor[,max_cor_mixID:=mixID[which.max(cor)],by=c("sample_name")]
mix_profiles_cor_max=mix_profiles_cor[mixID==max_cor_mixID][order(sample_name)]

#4. compare to RRBS results
RRBS_classification=fread("class_probs_annot_27_noNeural_predRRBS.tsv")
setnames(RRBS_classification,"sample","sample_name")
compare=merge(RRBS_classification,mix_profiles_cor_max,by="sample_name")
compare[,prob_dist:=as.numeric(dist(rbind(c(Classical,Proneural,Mesenchymal),c(cla/100,pro/100,mes/100)))),by="sample_name"]
#check for same classification
compare[sub_group==majority_sybtype]
write.table(compare,"subtype_validation_compare.tsv",sep="\t",row.names=FALSE,quote=FALSE)

##now analyze the comparison
compare[,plot_id:=paste0(patID,"_",ifelse(is.na(surgery.x),0,surgery.x)),]

pdf("subtype_validation_pred_quality.pdf",height=3,width=5,useDingbats=FALSE)
RRBS_cor=compare[,cor(prob_dist,auc),]
ggplot(compare,aes(x=auc,y=prob_dist))+geom_smooth(method="lm",col="black")+geom_point(size=2.5,stroke=1.5,shape=21,aes(col=majority_sybtype,fill=sub_group))+annotate("text",label=paste0("r=",round(RRBS_cor,3)),hjust=0,x=0.909,y=1.2)+xlab("RRBS based subtype prediction AUC")+ylab("Distance between RNA and RRBS\nbased subtype composition")+scale_fill_discrete(name="RRBS subtype")+scale_color_discrete(name="RNA subtype")
RNA_cor=compare[,cor(prob_dist,cor),]
ggplot(compare,aes(x=cor,y=prob_dist))+geom_smooth(method="lm",col="black")+geom_point(size=2.5,stroke=1.5,pch=21,aes(col=majority_sybtype,fill=sub_group))+annotate("text",label=paste0("r=",round(RNA_cor,3)),hjust=0,x=0.2,y=1.2)+xlab(" Maximum RNA profile correlation")+ylab("Distance between RNA and RRBS\nbased subtype composition")+scale_fill_discrete(name="RRBS subtype")+scale_color_discrete(name="RNA subtype")
dev.off()

#boxplot RNA classification quality vs. concordant/discordant classification 
####Figure S5d
compare[,subtype_classification:=ifelse(majority_sybtype==sub_group,"concordant","discordant"),]
wilcox_diff=compare[,wilcox.test(x=cor[subtype_classification=="concordant"],y=cor[subtype_classification=="discordant"],alternative="g"),]
ttest_diff=compare[,t.test(x=cor[subtype_classification=="concordant"],y=cor[subtype_classification=="discordant"],alternative="g"),]

pdf("subtype_validation_pred_quality_boxpl.pdf",height=3,width=4,useDingbats=FALSE)
ggplot(compare,aes(x=subtype_classification,y=cor))+geom_boxplot()+geom_point(position=position_jitter(width=0.25),alpha=0.8,size=2.5,stroke=1.5,shape=21,aes(col=majority_sybtype,fill=sub_group))+annotate("text",label=paste0("p.wil=",round(wilcox_diff$p.value,3),"\np.ttest=",round(ttest_diff$p.value,3)),hjust=0,x=1.5,y=0.7)+ylab(" Maximum RNA profile correlation")+scale_fill_discrete(name="RRBS subtype")+scale_color_discrete(name="RNA subtype")+stat_summary(fun.data=addNmin, geom="text", vjust=2, col="blue")
dev.off()


#plot correlations for each sample
RRBS_probs=compare[,c("sample_name","plot_id","Classical","Proneural", "Mesenchymal","sub_group","majority_sybtype"),with=FALSE]
setnames(RRBS_probs,c("Classical","Proneural", "Mesenchymal"),c("cla","pro","mes"))
RRBS_probs[,facet_label:=paste0(plot_id,"\nRNA:",majority_sybtype,"\nRRBS:",sub_group)]
mix_profiles_cor_RRBS=merge(mix_profiles_cor,RRBS_probs[,c("sample_name","facet_label"),],by="sample_name")

####Figure S5b
pdf("subtype_validation_correlations.pdf",height=16,width=18,useDingbats=FALSE)
ggtern(data=mix_profiles_cor_RRBS,aes(mes,cla,pro)) + geom_tri_tern(bins=15,fun=mean,aes(value=cor,fill=..stat..))+geom_point(data=RRBS_probs,col="black",size=4,shape=13)+facet_wrap(~facet_label,ncol=7)+scale_fill_gradient2(name="RNA profile\ncorrelation",high="green",mid="lightgrey",low="red") #colors changed to be colorblind friendly
dev.off()

#ROC curves
get_ROC=function(target,values){
  ROC=roc(target,values)
  dt=data.table(TPR=ROC$sensitivities,FPR=1-ROC$specificities,AUC=ROC$auc,run=0,rand=FALSE)
  dt=dt[order(FPR)]
  dt=dt[order(TPR)]
  
  set.seed(234)
  rand_labs=lapply(seq_len(100), function(x) sample(target))
  i=1
  for (rand_lab in rand_labs){
    print(i)
    ROC=roc(rand_lab,values)
    rand_dt=data.table(TPR=ROC$sensitivities,FPR=1-ROC$specificities,AUC=ROC$auc,run=i,rand=TRUE)
    rand_dt=rand_dt[order(FPR)]
    rand_dt=rand_dt[order(TPR)]
    i=i+1
    dt=rbindlist(list(dt,rand_dt))
  }
  dt_int=dt[,approx(FPR,TPR,xout=seq(0,1,0.01),yleft=0,yright=1,method="constant",ties=max),by=c("AUC","run","rand")]
  setnames(dt_int,c("x","y"),c("FPR","TPR"))
  dt_int[FPR==0,TPR:=0,]
  dt_int[FPR==1,TPR:=1,]
  
  return(list(roc_dt=dt,roc_dt_int=dt_int))
}

compare[,c("isMesenchymal","isClassical","isProneural"):=list(majority_sybtype=="Mesenchymal",majority_sybtype=="Classical",majority_sybtype=="Proneural"),]

#relaxed categories
compare[,c("isMesenchymal_rel","isClassical_rel","isProneural_rel"):=list(mes>30,cla>30,pro>30),]

compare_long=melt(compare,id.vars=c("sample_name","sub_group","Classical","Proneural","Mesenchymal"),measure.vars=c("isMesenchymal", "isClassical", "isProneural","isMesenchymal_rel","isClassical_rel","isProneural_rel"),value.name="is.subgroup",variable.name="subgroup_class")
compare_long[,class_prob:=ifelse(subgroup_class%in%c("isMesenchymal","isMesenchymal_rel"),Mesenchymal,ifelse(subgroup_class%in%c("isClassical","isClassical_rel"),Classical,Proneural)),]
compare_long[,Ntrue:=sum(is.subgroup),by="subgroup_class"]
compare_long[,Nfalse:=sum(!is.subgroup),by="subgroup_class"]

ROC_subgroups=compare_long[,get_ROC(is.subgroup,class_prob)$roc_dt,by=c("subgroup_class","Ntrue","Nfalse")]
ROC_subgroups[,facet_label:=paste0(subgroup_class,"\nAUC=",unique(round(mean(AUC[rand==FALSE]),2)),"/",round(mean(AUC[rand==TRUE]),2),"\nTRUE=",Ntrue,"\nFALSE=",Nfalse),by=c("subgroup_class")]

ROC_subgroups_int=compare_long[,get_ROC(is.subgroup,class_prob)$roc_dt_int,by=c("subgroup_class","Ntrue","Nfalse")]
ROC_subgroups_int[,facet_label:=paste0(subgroup_class,"\nAUC=",unique(round(mean(AUC[rand==FALSE]),2)),"/",round(mean(AUC[rand==TRUE]),2),"\nTRUE=",Ntrue,"\nFALSE=",Nfalse),by=c("subgroup_class")]

####Figure S5c
pdf("ROC_subgroups_rand.pdf",height=8,width=6)
ggplot(ROC_subgroups,aes(x=FPR,y=TPR,col=rand,group=run,alpha=rand))+geom_line()+facet_wrap(~facet_label,nrow=3,scale="free")+scale_color_manual(values=c("TRUE"="grey","FALSE"="blue"))+scale_alpha_manual(values=c("TRUE"=0.3,"FALSE"=1))

ggplot(ROC_subgroups_int,aes(x=FPR,y=TPR,group=rand))+stat_summary(geom="ribbon", fun.ymin = function(x) quantile(x, 0.025), fun.ymax = function(x) quantile(x, 0.975), fill="lightgrey",alpha=0.4)+stat_summary(geom="line",aes(col=rand), fun.y=mean)+facet_wrap(~facet_label,nrow=4,scale="free")+scale_color_manual(values=c("TRUE"="grey","FALSE"="blue"))
dev.off()




