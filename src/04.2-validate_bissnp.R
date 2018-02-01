library(project.init)
project.init2("GBMatch")
library("VariantAnnotation")

annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
annotation_val=annotation[category=="GBmatch_valF"]

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/05-bissnp")
dir.create(out_dir)
setwd(out_dir)

#get RRBS varinats
rrbs_var_samp=fread("bissnp_var_samp.tsv")
rrbs_var_samp_long=melt(rrbs_var_samp,id.vars=c( "patID","sample_name", "surgery", "category", "bg_calls"),variable.name="count_type",value.name="rrbs_var_count")
setnames(rrbs_var_samp_long,c( "bg_calls"),c("rrbs_bg_calls"))


##get WGS variant files
annotation_val[,vcfs:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_vcf"),pattern=paste0(sub("_fv","",N_number_seq),".*snpeff.vcf$"),full.names=TRUE),by=1:nrow(annotation_val)]

annotation_val[,vcf_bg:=list.files(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_bam"),pattern=paste0(sub("_fv","",N_number_seq),".*\\.cov$"),full.names=TRUE),by=1:nrow(annotation_val)]

annotation_val=annotation_val[!is.na(vcfs)]


condense_vcf=function(sample,patient,vcf_file,vcf_bg_file,covThres=5){
  #load data
  message(sample)
  vcf=readVcf(vcf_file,"hg38")
  bg=as.numeric(system(paste0("grep genome ",vcf_bg_file," | awk '{if($2 > ",covThres,") sum += $3}END{print sum}'"),intern=TRUE))
  
  if (nrow(vcf@fixed)==0){
    message("No SNPs found. Retuning NAs")
    dt=data.table(CHROM=NA,POS=NA,varID=NA,REF=NA,ALT=NA,QUAL=NA,H_count=NA,M_count=NA,H_genes=NA,M_genes=NA,sample_name=sample,patID=patient,bg_calls=bg)
    return(dt)
  }
  
  #function to farse the annotation
  parse_annot=function(annotations){
    count_H=sum(grepl("HIGH",annotations))
    count_M=sum(grepl("MODERATE",annotations))
    genes_H=paste0(unique(unlist(lapply(grep("HIGH",annotations,value=TRUE),function(x){unlist(strsplit(x,"\\|"))[5]}))),collapse=",")
    genes_M=paste0(unique(unlist(lapply(grep("MODERATE",annotations,value=TRUE),function(x){unlist(strsplit(x,"\\|"))[5]}))),collapse=",")  
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
  
  
  df=cbind(pos, vcf@fixed,vcf.info.ANN=unlist(vcf@info$EFF),vcf.info.DP=as.numeric(vcf@info$DP))
  
  df=df[df$vcf.info.DP>covThres,]
  
  dt=data.table(CHROM=as.character(df$seqnames),POS=as.numeric(df$start),varID=as.character(df$varID),REF=as.character(df$REF),ALT=unlist(lapply(df$ALT,paste0,collapse=",")),QUAL=df$QUAL,rbindlist(lapply(df$vcf.info.ANN,parse_annot)))
  dt[,sample_name:=sample,]
  dt[,patID:=patient,]  #enable for complete rerun
  dt[,bg_calls:=bg,]
  return(dt)
}

covThres=2
wgs_var=data.table()
for(i in c(1:nrow(annotation_val))){
  cachename=paste0("wgs_",annotation_val$N_number_seq[i])
  message(cachename)
  simpleCache(recreate=TRUE,cachename,{dt=condense_vcf(sample=annotation_val$N_number_seq[i],patient=annotation_val$patID[i],vcf_file=annotation_val$vcfs[i],vcf_bg_file=annotation_val$vcf_bg[i],covThres=covThres);return(dt)},cacheSubDir=paste0("wgs_var_cov",covThres))
  
  if (!"data.table"%in%is(get(cachename))){next}
  
  var=get(cachename)
  var[,patID:=annotation_val$patID[i],]
  var[,surgery:=annotation_val$surgery.x[i],]
  var[,category:=annotation_val$category[i],]
  wgs_var=rbindlist(list(wgs_var,var))
  
}

wgs_var_samp=wgs_var[QUAL>=20,list(all_count=.N,H_count=sum(H_count>0),M_count=sum(M_count>0),bg_calls=bg_calls[1],category=category[1]),by=c("patID","sample_name","surgery")]
wgs_var_samp_long=melt(wgs_var_samp,id.vars=c( "patID","sample_name", "surgery", "category", "bg_calls"),variable.name="count_type",value.name="wgs_var_count")
setnames(wgs_var_samp_long,c( "bg_calls"),c("wgs_bg_calls"))

var_samp_compare=merge(wgs_var_samp_long,rrbs_var_samp_long,by=c("patID","sample_name","surgery","category","count_type"))
var_samp_compare[,c("wgs_var_count_norm","rrbs_var_count_norm"):=list(wgs_var_count/wgs_bg_calls*1000000,rrbs_var_count/rrbs_bg_calls*1000000),]

var_samp_compare[,cor(wgs_var_count_norm,rrbs_var_count_norm),by=count_type]
ggplot(var_samp_compare,aes(x=wgs_var_count_norm,y=rrbs_var_count_norm))+geom_point()+geom_smooth(method="lm")+facet_wrap(~count_type,scale="free")






