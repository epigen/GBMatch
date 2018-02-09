library(project.init)
project.init2("GBMatch")
library("VariantAnnotation")

annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/04-bissnp")
dir.create(out_dir)
setwd(out_dir)

annotation[,vcfs:=paste0(getOption("PROCESSED.PROJECT"),"results_pipeline/",N_number_seq,"/bissnp_hg38/",N_number_seq,"_rg.snp.filtered.sort_annot.vcf"),]
annotation[,vcf_bg:=paste0(getOption("PROCESSED.PROJECT"),"results_pipeline/",N_number_seq,"/bissnp_hg38/",N_number_seq,"_rg.cpg.raw.vcf.MethySummarizeList.txt"),]

stopifnot(file.exists(annotation$vcfs)&file.exists(annotation$vcf_bg))


condense_vcf=function(sample,patient,vcf_file,vcf_bg_file){
  #load data
  message(sample)
  vcf=readVcf(vcf_file,"hg38")
  bg=as.numeric(gsub("[A-za-z ]","",grep("Confidently called bases[ ]*[0-9]",readLines(vcf_bg_file),value=TRUE,perl=TRUE)))

  if (nrow(vcf@fixed)==0){
    message("No SNPs found. Retuning NAs")
    dt=data.table(CHROM=NA,POS=NA,varID=NA,REF=NA,ALT=NA,QUAL=NA,H_count=NA,M_count=NA,H_genes=NA,M_genes=NA,sample_name=sample,patID=patient,bg_calls=bg)
    return(dt)
  }
  
  #function to farse the annotation
  parse_annot=function(annotations){
  count_H=sum(grepl("HIGH",annotations))
  count_M=sum(grepl("MODERATE",annotations))
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
  
  
  df=cbind(pos, vcf@fixed,vcf.info.ANN=vcf@info$ANN)
  dt=data.table(CHROM=as.character(df$seqnames),POS=as.numeric(df$start),varID=as.character(df$varID),REF=as.character(df$REF),ALT=unlist(lapply(df$ALT,paste0,collapse=",")),QUAL=df$QUAL,rbindlist(lapply(df$vcf.info.ANN,parse_annot)))
  dt[,sample_name:=sample,]
  dt[,patID:=patient,]  #enable for complete rerun
  dt[,bg_calls:=bg,]
  return(dt)
}


bissnp_var=data.table()
for(i in c(1:nrow(annotation))){
  cachename=paste0("bissnp_",annotation$N_number_seq[i])
  message(cachename)
  simpleCache(recreate=FALSE,cachename,{dt=condense_vcf(sample=annotation$N_number_seq[i],patient=annotation$patID[i],vcf_file=annotation$vcfs[i],vcf_bg_file=annotation$vcf_bg[i]);return(dt)},cacheSubDir="bissnp_var")
  
  if (!"data.table"%in%is(get(cachename))){next}
  
  var=get(cachename)
  var[,patID:=annotation$patID[i],]
  var[,surgery:=annotation$surgery.x[i],]
  var[,category:=annotation$category[i],]
  bissnp_var=rbindlist(list(bissnp_var,var))
  
}


bissnp_var_samp=bissnp_var[,list(all_count=.N,H_count=sum(H_count>0),M_count=sum(M_count>0),bg_calls=bg_calls[1],category=category[1]),by=c("patID","sample_name","surgery")]


bissnp_var_samp_long=melt(bissnp_var_samp,id.vars=c("patID", "sample_name","surgery","bg_calls","category"))
bissnp_var_samp_long[,mut_norm:=value/bg_calls*1000000,]

bissnp_var_pat=bissnp_var_samp[category%in%c("GBMatch","GBmatch_val"),list(all_count=all_count[which.max(bg_calls)], H_count=H_count[which.max(bg_calls)], M_count=M_count[which.max(bg_calls)], bg_calls=bg_calls[which.max(bg_calls)]),by=c("patID","surgery","category")]
bissnp_var_pat_long=melt(bissnp_var_pat,id.vars=c("patID","surgery","category","bg_calls"))
bissnp_var_pat_long[,mut_norm:=value/bg_calls*1000000,]

annot_age=unique(annotation[,list(age=age,surgery=surgery.x,patID=patID,IDH=IDH,sex=sex,timeToSecSurg=timeToSecSurg),])
bissnp_var_pat_annot=merge(bissnp_var_pat_long,annot_age,by=c("patID","surgery"))

bissnp_var_pat_annot[surgery%in%c(1,2),list(cor=cor(mut_norm,age)),by=c("variable","sex","IDH","surgery","category")]
bissnp_var_pat_annot[,cor(bg_calls,value),by=c("variable")]

pdf("norm_mutvsAge.pdf",height=4,width=10.5)
ggplot(bissnp_var_pat_annot[surgery%in%c(1,2)&IDH=="wt"],aes(x=age,y=mut_norm,col=category,group=category))+geom_point(aes(fill=sex),pch=21,alpha=0.6)+geom_smooth(method="lm")+scale_color_manual(values=c("GBMatch"="grey","GBmatch_val"="black"))+facet_wrap(~variable,scale="free")
dev.off()

pdf("norm_mutvsTime.pdf",height=4,width=10.5)
ggplot(bissnp_var_pat_annot[surgery%in%c(1,2)&IDH=="wt"],aes(x=timeToSecSurg,y=mut_norm))+geom_point(aes(fill=sex),pch=21,alpha=0.6)+geom_smooth(method="lm")+facet_wrap(~variable,scale="free")
dev.off()

write.table(bissnp_var_samp,"bissnp_var_samp.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(bissnp_var_pat,"bissnp_var_pat.tsv",sep="\t",quote=FALSE,row.names=FALSE)
#bissnp_var_pat=fread("bissnp_var_pat.tsv")
write.table(bissnp_var,"bissnp_var_combined.tsv",sep="\t",quote=FALSE,row.names=FALSE)

bissnp_var_pat_long[,variable:=factor(variable,levels=c("all_count","M_count","H_count")),]
pdf("norm_mut.pdf",height=3,width=7)
ggplot(bissnp_var_pat_long[surgery<3],aes(y=mut_norm,x=paste0(category,"\nSurgery:",surgery)))+geom_point(col="grey",alpha=1,pch=21,position=position_jitter(width=0.2))+geom_boxplot(outlier.shape=NA,fill="transparent")+facet_wrap(~variable,scale="free_y")+xlab("")+ylab("Mutations per 1M called bases")
dev.off()

bissnp_var_pat_long_wide=reshape(bissnp_var_pat_long[!is.na(surgery)],idvar=c("patID","variable"),timevar="surgery",direction="wide")
bissnp_var_pat_wide=reshape(bissnp_var_pat_long[!is.na(surgery)],idvar=c("patID","bg_calls","surgery"),timevar=c("variable"),drop="value",direction="wide")
bissnp_var_pat_wide=reshape(bissnp_var_pat_wide,idvar="patID",timevar="surgery",direction="wide")


pdf("norm_mut_1vs2.pdf",height=3,width=9)
ggplot(bissnp_var_pat_long_wide,aes(y=mut_norm.2,x=mut_norm.1))+geom_abline(slope=1,intercept=c(0,0),size=2,col="grey",alpha=0.7)+geom_point(alpha=0.5,pch=21)+facet_wrap(~variable,scale="free")+xlab("")+ylab("Mutations per 1M called bases surgery 2")+xlab("Mutations per 1M called bases surgery 1")
dev.off()

pdf("norm_mut_1vs2_col.pdf",height=4,width=6)
ggplot(bissnp_var_pat_wide,aes(y=mut_norm.all_count.2,x=mut_norm.all_count.1,fill=log2((mut_norm.H_count.1+1)/(mut_norm.H_count.2+1))))+geom_abline(slope=1,intercept=c(0,0),size=2,col="grey",alpha=0.7)+geom_point(size=3,alpha=0.8,pch=21)+xlab("")+ylab("Recurring tumor")+xlab("Primary tumor")+scale_fill_gradient2(low="green",high="red",mid="white",name="log2(impact_1/impact_2)")+xlim(c(600,1600))+ylim(c(600,1600))+ggtitle("Mutations per 1M called bases")+theme(aspect.ratio=1)
dev.off()


#check mutated genes
GOIs_sel=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/Known_GB_relevant_Genes_Nik.tsv"))
setnames(GOIs_sel,"x","gene")
OGTS=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/TS-OG_1235122TablesS1-4_vogelsteinetal_2013_red.csv"))
GOIs_sel=merge(GOIs_sel,OGTS,by="gene",all.x=TRUE)


bissnp_var_hit=bissnp_var[!is.na(H_genes)]

#recycling should be ok, because some lines are expanded
bissnp_var_hit_long=bissnp_var_hit[, cbind(gene=unlist(strsplit(x=H_genes,",")), .SD), by=c("sample_name", "bg_calls","patID", "surgery", "category","REF","ALT","QUAL")]

bissnp_var_hit_long=merge(bissnp_var_hit_long,GOIs_sel,by="gene")
bissnp_var_hit_long[QUAL>100,QUAL:=100,]
bissnp_var_hit_long[,recurrent:=ifelse(1%in%surgery&2%in%surgery,TRUE,FALSE),by=c("patID","gene")]

unique(bissnp_var_hit_long$gene)

#needs to be SVG because in pdf the dots don't have the right size
pdf("high_GOI.pdf",height=4.5,width=9,useDingbats=FALSE)
ggplot(bissnp_var_hit_long[surgery%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")],aes(y=gene,x=as.factor(surgery),fill=classification,alpha=QUAL))+geom_dotplot(binaxis="y",stackdir = "center",dotsize=0.8,binwidth=1)+scale_fill_manual(values=c("ND"="black","DRG"="blue","OG"="red","TSG"="green"))+scale_alpha_continuous()+coord_flip()+theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))+theme(legend.position="bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~category,nrow=2,scale="free_y")
dev.off()



