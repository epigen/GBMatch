library(project.init)
project.init2("GBMatch")
library(simpleCache)
library(pheatmap)
library(ggtern)

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
row.names(TCGA_genes_mat)=TCGA_genes$gene
col_annot=data.frame(subtype=TCGA_samples_annot$sample_subtype,row.names=TCGA_samples_annot$TCGA_sample)
row_annot=data.frame(subtype=TCGA_genes_annot$gene_subtype,row.names=TCGA_genes_annot$gene)

pheatmap(TCGA_genes_mat, annotation_col=col_annot,annotation_row=row_annot)

#mean expression of each gene in each subgroup
subtype_expression=TCGA_genes_long[,list(mean_expr=mean(value),sd_expr=sd(value)),by=c("sample_subtype","gene","gene_subtype")]
subtype_expression[,gene_subtype_ass:=sample_subtype[which.max(mean_expr)],by="gene"]
concordant_genes=unique(subtype_expression[gene_subtype==gene_subtype_ass]$gene)

pheatmap(TCGA_genes_mat[row.names(TCGA_genes_mat)%in%concordant_genes,], annotation_col=col_annot,annotation_row=row_annot[row.names(row_annot)%in%concordant_genes,,drop=FALSE])

RNA_data[gene%in%concordant_genes]

#some of the genes (~50) have different names in the TCGA data set (for now ignore, add manually later)
concordant_genes[!concordant_genes%in%unique(RNA_data_uniq$gene)]


#Now cluster data according to selected and annotated genes
#discard neural genes
concordant_genes_noNeur=unique(subtype_expression[gene_subtype==gene_subtype_ass&gene_subtype!="Neural"&gene%in%RNA_data_uniq$gene]$gene)

RNA_data_uniq[,RPKM_scaled:=scale(log(RPKM),center=TRUE,scale=TRUE),by="gene"]

RNA_data_uniq_wide=reshape(RNA_data_uniq,idvar="gene",timevar="sample_name",drop="RPKM",direction="wide")

RNA_data_uniq_mat=as.matrix(RNA_data_uniq_wide[gene%in%concordant_genes_noNeur,-c("gene"),with=FALSE])
row.names(RNA_data_uniq_mat)=RNA_data_uniq_wide[gene%in%concordant_genes_noNeur]$gene

pheatmap(RNA_data_uniq_mat,annotation_row=row_annot[row.names(row_annot)%in%row.names(RNA_data_uniq_mat),,drop=FALSE])

#TCGA data common geans available in both datasets, only neural
pheatmap(TCGA_genes_mat[row.names(TCGA_genes_mat)%in%concordant_genes_noNeur,!colnames(TCGA_genes_mat)%in%row.names(col_annot[col_annot$subtype=="Neural",,drop=FALSE])], annotation_col=col_annot[col_annot$subtype!="Neural",,drop=FALSE],annotation_row=row_annot[row.names(row_annot)%in%concordant_genes_noNeur,,drop=FALSE])


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


compare[sub_group==majority_sybtype]

compare[,prob_dist:=dist(rbind(c(Classical,Proneural,Mesenchymal),c(cla/100,pro/100,mes/100))),by="sample_name"]


#plot correlations for each sample
RRBS_probs=compare[,c("sample_name","Classical","Proneural", "Mesenchymal","sub_group","majority_sybtype"),with=FALSE]
setnames(RRBS_probs,c("Classical","Proneural", "Mesenchymal"),c("cla","pro","mes"))
RRBS_probs[,facet_label:=paste0(sample_name,"\nRNA:",majority_sybtype,"\nRRBS:",sub_group)]
mix_profiles_cor_RRBS=merge(mix_profiles_cor,RRBS_probs[,c("sample_name","facet_label"),],by="sample_name")


pdf("subtype_validation_correlations.pdf",height=16,width=18)
ggtern(data=mix_profiles_cor_RRBS,aes(mes,cla,pro)) + geom_tri_tern(bins=15,fun=mean,aes(value=cor,fill=..stat..))+geom_point(data=RRBS_probs,col="black",size=4,shape=13)+facet_wrap(~facet_label,ncol=7)+scale_fill_gradient2(name="RNA profile\ncorrelation",high="green",mid="lightgrey",low="red")
dev.off()

#[sample_name%in%c("N940_16_fv","N1290_12_fv","N2207_13_fv")]


