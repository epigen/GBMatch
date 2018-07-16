#NOTE: This script compares mutational load as determined from RRBS data with mutational load as determined from low coverage WGS data.
#Two different approaches for WGS: based on VCF (commented out) or summary stats file(active).

library(project.init)
project.init2("GBMatch")
library("VariantAnnotation")

annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
annotation_val=annotation[category=="GBmatch_valF"]

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/04-bissnp")
dir.create(out_dir)
setwd(out_dir)


##get WGS and RRBS variant files
#for bcftools
#annotation_val[,vcfs:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_vcf"),pattern=paste0(sub("_fv","",N_number_seq),".*snpeff.vcf$"),full.names=TRUE),by=1:nrow(annotation_val)]
#annotation_val[,vcf_bg:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_bam"),pattern=paste0(sub("_fv","",N_number_seq),".*\\.cov$"),full.names=TRUE),by=1:nrow(annotation_val)]
#wgs stats from bcftools (to compare on a rank basis)
annotation_val[,wgs_stats:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_vcf"),pattern=paste0(sub("_fv","",N_number_seq),".*stats.csv$"),full.names=TRUE),by=1:nrow(annotation_val)]

#for bissnp
annotation_val[,vcfs:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_bissnp/"),pattern=paste0(sub("_fv","",N_number_seq),".*_rg.snp.filtered.sort_annot.vcf$"),full.names=TRUE,recursive=TRUE),by=1:nrow(annotation_val)]
annotation_val[,vcf_bg:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_bissnp"),pattern=paste0(sub("_fv","",N_number_seq),".*cpg.raw.vcf.MethySummarizeList.txt$"),full.names=TRUE,recursive=TRUE),by=1:nrow(annotation_val)]

annotation_val[,vcfs_RRBS:=list.files(file.path(getOption("PROCESSED.PROJECT"),"results_pipeline/"),pattern=paste0(N_number_seq,".*_rg.snp.filtered.sort_annot.vcf$"),full.names=TRUE,recursive=TRUE),by=1:nrow(annotation_val)]
annotation_val[,vcf_bg_RRBS:=list.files(file.path(getOption("PROCESSED.PROJECT"),"results_pipeline/"),pattern=paste0(N_number_seq,".*cpg.raw.vcf.MethySummarizeList.txt$"),full.names=TRUE,recursive=TRUE),by=1:nrow(annotation_val)]

annotation_val=annotation_val[!is.na(vcfs)]

##FUNCTIOONS
read_stats=function(stats_file){
  variants=as.numeric(gsub("[A-za-z ,]","",grep("Number_of_variants_processed , [0-9]",readLines(stats_file),value=TRUE,perl=TRUE)))
  effgenome=as.numeric(gsub("[A-za-z ,]","",grep("Genome_effective_length ,[0-9]",readLines(stats_file),value=TRUE,perl=TRUE)))
  rate=as.numeric(gsub("[A-za-z ,]","",grep("Change_rate , [0-9]",readLines(stats_file),value=TRUE,perl=TRUE)))
  return(list(wgs_var=variants,wgs_effgenome=effgenome,wgs_changerate=rate))
}


condense_vcf=function(sample,patient,vcf_file,vcf_bg_file,covThres=5){
  #load data
  message(sample)
  vcf=readVcf(vcf_file,"hg38")
  
  #for bcftools
  #bg=c()
  #for (i in c(0:20)){
  #bg_temp=as.numeric(system(paste0("grep genome ",vcf_bg_file," | awk '{if($2 > ",i,") sum += $3}END{print sum}'"),intern=TRUE))
  #bg=paste0(bg,"|",bg_temp)
  #}
  bg_conf=as.numeric(gsub("[A-za-z ]","",grep("Confidently called bases[ ]*[0-9]",readLines(vcf_bg_file),value=TRUE,perl=TRUE)))
  bg_all=as.numeric(gsub("[A-za-z ]","",grep("Callable bases[ ]*[0-9]",readLines(vcf_bg_file),value=TRUE,perl=TRUE)))
  bg_reads=as.numeric(gsub("[A-za-z ]","",grep("Average good reads coverage in callable position   [ ]*[0-9]",readLines(vcf_bg_file),value=TRUE,perl=TRUE)))
  
  if (nrow(vcf@fixed)==0){
    message("No SNPs found. Retuning NAs")
    dt=data.table(CHROM=NA,POS=NA,varID=NA,REF=NA,ALT=NA,QUAL=NA,DP=NA,DP4_ref=NA,DP4_alt=NA,DP4=NA,H_count=NA,M_count=NA,H_genes=NA,M_genes=NA,sample_name=sample,patID=patient,bg_calls_all=bg_all,bg_calls_conf=bg_conf,bg_calls_reads=bg_reads)
    return(dt)
  }
  
  #function to parse the annotation
  parse_annot=function(annotations){
    count_H=sum(grepl("HIGH",annotations))
    count_M=sum(grepl("MODERATE",annotations))
    #for bcftools
    #genes_H=paste0(unique(unlist(lapply(grep("HIGH",annotations,value=TRUE),function(x){unlist(strsplit(x,"\\|"))[5]}))),collapse=",")
    #genes_M=paste0(unique(unlist(lapply(grep("MODERATE",annotations,value=TRUE),function(x){unlist(strsplit(x,"\\|"))[5]}))),collapse=",")  
    
    genes_H=paste0(unique(unlist(lapply(grep("HIGH",annotations,value=TRUE),function(x){unlist(strsplit(x,"\\|"))[4]}))),collapse=",")
    genes_M=paste0(unique(unlist(lapply(grep("MODERATE",annotations,value=TRUE),function(x){unlist(strsplit(x,"\\|"))[4]}))),collapse=",")  
    
    ret=data.table(H_count=count_H,M_count=count_M,H_genes=genes_H,M_genes=genes_M)
    ret[ret==""]=NA
    return(ret)
  }
  
  pos=as.data.frame(vcf@rowRanges,row.names=NULL)
  pos$varID=names(vcf@rowRanges)
  row.names(pos)=NULL
  pos$paramRangeID=NULL
  pos$width=NULL
  pos$strand=NULL
  
  #for bcftools
  #DP4_dt=as.data.table(vcf@info$DP4)
  #DP4_dt_red=DP4_dt[,list(ref=sum(value[c(1,2)]),alt=sum(value[c(3,4)]),all=sum(value)),by="group"]
  #df=cbind(pos, vcf@fixed,vcf.info.ANN=unlist(vcf@info$EFF),vcf.info.DP=as.numeric(vcf@info$DP),vcf.info.DP4_ref=DP4_dt_red$ref,vcf.info.DP4_alt=DP4_dt_red$alt,vcf.info.DP4=DP4_dt_red$all)
  #df=df[df$vcf.info.DP>covThres,]
  #dt=data.table(CHROM=as.character(df$seqnames),POS=as.numeric(df$start),varID=as.character(df$varID),REF=as.character(df$REF),ALT=unlist(lapply(df$ALT,paste0,collapse=",")),QUAL=df$QUAL,DP=df$vcf.info.DP,DP4_ref=df$vcf.info.DP4_ref,DP4_alt=df$vcf.info.DP4_alt,DP4=df$vcf.info.DP4,rbindlist(lapply(df$vcf.info.ANN,parse_annot)))
  
  df=cbind(pos, vcf@fixed,vcf.info.ANN=vcf@info$ANN)
  dt=data.table(CHROM=as.character(df$seqnames),POS=as.numeric(df$start),varID=as.character(df$varID),REF=as.character(df$REF),ALT=unlist(lapply(df$ALT,paste0,collapse=",")),QUAL=df$QUAL,rbindlist(lapply(df$vcf.info.ANN,parse_annot)))
  
  dt[,sample_name:=sample,]
  dt[,patID:=patient,]  #enable for complete rerun
  dt[,c("bg_calls_all","bg_calls_conf","bg_calls_reads"):=list(bg_all,bg_conf,bg_reads),]
  return(dt)
}

#read wgs stats
annotation_val[,c("wgs_bcf_var","wgs_bcf_effgenome","wgs_bcf_changerate"):=read_stats(wgs_stats),1:nrow(annotation_val)]
annotation_val[,wgs_bcf_cpm:=wgs_bcf_var/wgs_bcf_effgenome*1000000]


#Perform analysis
#chose mode (rrbs or wgs) run the following once on each to generate both stats (necessary to compare)
modes=c("rrbs","wgs")

for (mode in modes){
  covThres=0 #irrelevant for bissnp
  
  if (mode=="rrbs"){rrbs_var=data.table()}
  if (mode=="wgs"){wgs_var=data.table()}
  
  for(i in c(1:nrow(annotation_val))){
    if (mode=="wgs"){
      cachename=paste0("wgs_",annotation_val$N_number_seq[i])
      message(cachename)
      ##for bcftools
      ##  simpleCache(recreate=FALSE,cachename,{dt=condense_vcf(sample=annotation_val$N_number_seq[i],patient=annotation_val$patID[i],vcf_file=annotation_val$vcfs[i],vcf_bg_file=annotation_val$vcf_bg[i],covThres=covThres);return(dt)},cacheSubDir=paste0("wgs_var_cov",covThres))
      #for bissnp
      simpleCache(recreate=FALSE,cachename,{dt=condense_vcf(sample=annotation_val$N_number_seq[i],patient=annotation_val$patID[i],vcf_file=annotation_val$vcfs[i],vcf_bg_file=annotation_val$vcf_bg[i],covThres=covThres);return(dt)},cacheSubDir=paste0("wgs_var_bissnp_cov",covThres))  
    }
    if (mode=="rrbs"){
      cachename=paste0("rrbs_",annotation_val$N_number_seq[i])
      message(cachename)
      simpleCache(recreate=FALSE,cachename,{dt=condense_vcf(sample=annotation_val$N_number_seq[i],patient=annotation_val$patID[i],vcf_file=annotation_val$vcfs_RRBS[i],vcf_bg_file=annotation_val$vcf_bg_RRBS[i],covThres=covThres);return(dt)},cacheSubDir=paste0("rrbs_var_bissnp_cov",covThres))
    }
    if (!"data.table"%in%is(get(cachename))){next}
    
    var=get(cachename)
    anno_vars=c("patID","surgery.x","category","wgs_bcf_changerate","wgs_bcf_cpm","Unique_CpGs","enrichmentCycles", "Aligned_reads", "Total_CpGs", "meanCoverage", "bisulfiteConversionRate","rand_frag_perc")
    vals=annotation_val[i,anno_vars,with=FALSE]
    var[,(anno_vars):=vals,]
    
    if (mode=="wgs"){wgs_var=rbindlist(list(wgs_var,var))}
    if (mode=="rrbs"){rrbs_var=rbindlist(list(rrbs_var,var))}
  }
  
  #postprocess wgs and rrbs calls
  
  if (mode=="wgs"){
    #for bcftools
    #wgs_var[,DP:=as.numeric(DP),]
    #covThres_sec=5
    #wgs_var_samp=wgs_var[QUAL>=20&DP>covThres_sec,list(all_count=.N,H_count=sum(H_count>0),M_count=sum(M_count>0),bg_calls=as.numeric(unlist(strsplit(bg_calls[1],"\\|"))[covThres_sec+2]),category=category[1]),by=c("patID","sample_name","surgery")]
    
    wgs_var_samp=wgs_var[,list(all_count=.N,H_count=sum(H_count>0),M_count=sum(M_count>0),wgs_bg_calls=as.numeric(bg_calls_conf[1]),wgs_bg_calls_reads=as.numeric(bg_calls_reads[1]),category=category[1]),by=c("patID","sample_name","surgery.x")] #needs to be adapted upon rerun (include bg_calls_reads and bg_calls --> bg_calls_conf)
    wgs_var_samp_long=melt(wgs_var_samp,id.vars=c( "patID","sample_name", "surgery.x", "category", "wgs_bg_calls","wgs_bg_calls_reads"),variable.name="count_type",value.name="wgs_var_count")
    
  }
  
  if (mode=="rrbs"){
    rrbs_var_samp=rrbs_var[,list(all_count=.N,H_count=sum(H_count>0),M_count=sum(M_count>0),rrbs_bg_calls=as.numeric(bg_calls_conf[1]),rrbs_bg_calls_reads=as.numeric(bg_calls_reads[1])),by=c("sample_name",anno_vars)]
    rrbs_var_samp_long=melt(rrbs_var_samp,id.vars=c("sample_name","rrbs_bg_calls","rrbs_bg_calls_reads",anno_vars),variable.name="count_type",value.name="rrbs_var_count")
    
  }
}
#combine and compare wgs and rrbs stats 

var_samp_compare=merge(wgs_var_samp_long,rrbs_var_samp_long,by=c("patID","sample_name","surgery.x","category","count_type"))
var_samp_compare[,c("wgs_var_count_norm","rrbs_var_count_norm"):=list(wgs_var_count/wgs_bg_calls*1000000,rrbs_var_count/rrbs_bg_calls*1000000),]

var_samp_compare[,wgs_bcf_cpm_rank:=rank(wgs_bcf_cpm,ties.method="random"),by=count_type]
var_samp_compare[,rrbs_cpm_rank:=rank(rrbs_var_count_norm,ties.method="random"),by=count_type]
var_samp_compare[,wgs_cpm_rank:=rank(wgs_var_count_norm,ties.method="random"),by=count_type]

cors=unique(var_samp_compare[,cor.test(wgs_var_count_norm,rrbs_var_count_norm,alternative="greater"),by=count_type][,c("count_type","p.value","estimate"),])
pdf("validation_wgs_rrbs_cor.pdf",height=2.5,width=7)
ggplot(var_samp_compare,aes(x=wgs_var_count_norm,y=rrbs_var_count_norm))+geom_point(shape=21,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cors,x=-Inf,y=Inf,aes(label=paste0("r= ",signif(unique(estimate),2),"\np= ",signif(unique(p.value),2))),hjust=-0.2,vjust=1.2)+facet_wrap(~count_type,scale="free")+ylab("Variants per Million (RRBS)")+xlab("Variants per Million (WGS)")
dev.off()

pdf("validation_wgs_rrbs_rank.pdf",height=3,width=8)
ggplot(var_samp_compare,aes(wgs_bcf_cpm_rank,rrbs_cpm_rank))+geom_point()+facet_wrap(~count_type)
ggplot(var_samp_compare,aes(wgs_bcf_cpm_rank,wgs_cpm_rank))+geom_point()+facet_wrap(~count_type)
dev.off()

write.table(var_samp_compare,"validation_compare.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rrbs_var,"validation_rrbs_var.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(wgs_var,"validation_wgs_var.tsv",sep="\t",quote=FALSE,row.names=FALSE)
