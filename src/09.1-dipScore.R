library(project.init)
project.init2("GBMatch")
source(file.path(getOption("PROJECT.DIR"),"src/99-MIRA.R"))
library(LOLA)

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/09-dipScore")
dir.create(out_dir)
setwd(out_dir)

eload(loadAnnotation())
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
transc_st=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/08.1-GBM_classifier/class_probs_annot_27_noNeural_predRRBS.tsv"))
transc_st=transc_st[!is.na(sub_group)]
transc_st[,sub_group_prob:=get(sub_group),by=sample]
setnames(transc_st, "sample","N_number_seq")
annotation=merge(annotation,transc_st[,-c("patID","IDH","surgery.x","Follow-up_years","WHO2016_classification","category"),],by="N_number_seq")

setkey(annotation, N_number_seq)
setkey(SV$psa,sample_name)
SV$psa[,cell_type:=NULL,]
SV$psa[annotation,cell_type:=ifelse(category=="control","white matter","glioblastoma"),allow=TRUE]
SV$psa[annotation,category:=category,allow=TRUE]


loadCaches("rrbsCg")

#potential new data
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62193

regionDB = loadRegionDB(file.path(Sys.getenv("RESOURCES"),"regions/LOLACore/hg38/"))
cellType_anno=fread(file.path(getOption("PROJECT.DIR"),"metadata/LOLA_annot/CellTypes.tsv"))
regionDB$regionAnno[,GRL_id:=1:nrow(regionDB$regionAnno),]
regionDB$regionAnno=merge(regionDB$regionAnno,cellType_anno,by="cellType",all.x=TRUE)


factor="Astrocyte|ESC"
total_width=5000

#select single factors/cell types
mylist=lapply(c(regionDB$regionAnno[size>1000][grep(factor,antibody)]$GRL_id,regionDB$regionAnno[size>1000][grep(factor,cellType_corr)]$GRL_id),function(x){name=paste0(regionDB$regionAnno[GRL_id==x]$antibody,"__",regionDB$regionAnno[GRL_id==x]$cellType,"_",regionDB$regionAnno[GRL_id==x]$treatment,"_",regionDB$regionAnno[GRL_id==x]$size);dt=as.data.table(as.data.frame(regionDB$regionGRL[[x]]));dt[,start:=ceiling(start-(total_width-width)/2),];dt[,end:=ceiling(end+(total_width-width)/2),];dt[,width:=end-start,];return(list("name"=name,"dt"=dt))})


rangesList=lapply(mylist,function(x){x$dt})
names(rangesList)=as.character(unlist(lapply(mylist,function(x){x$name})))
attach(rangesList)

binCount=21

binResults=binProcess(rangesList, list("rrbsCg"), annotation, binCount)


#plot profiles by subgroup overview (all data)
pdf(paste0("aggregatedMeth_subgroup_",factor,"new.pdf"),height=6,width=8)
for(assay in names(binResults)){
  pl=ggplot(binResults[[assay]]$binnedBSDT[(readCount>1000&(auc>0.8&sub_group_prob>0.8)|category=="control")],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+facet_wrap(~category,scale="free")+ggtitle(assay)+scale_x_discrete(labels=labelBinCenter(binCount))+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))
 print(pl)
}
dev.off()

#plot profiles by subgroup only GBMatch IDH=="wt"
pdf(paste0("aggregatedMeth_subgroup_GBMatch_",factor,".pdf"),height=3.5,width=4)
for(assay in names(binResults)){
  spl=unlist(strsplit(gsub("_\\([0-9]*\\)","",assay),"_"))
  ab=spl[1]
  cellType=ifelse(grepl("embryonic|ESC",spl[3],ignore.case=TRUE),"ESC",ifelse(grepl("Astro|NH-A",spl[3],ignore.case=TRUE),"Astrocytes",spl[3]))
  annot=paste0(ab," (",cellType,")")
  pl=ggplot(binResults[[assay]]$binnedBSDT[readCount>1000&auc>0.8&sub_group_prob>0.8&category=="GBMatch"&IDH=="wt"],aes(x=factor(regionGroupID),y=methyl,col=sub_group))+geom_line(alpha=0.2,aes(group=id))+geom_smooth(aes(group=sub_group), se=FALSE)+annotate(geom="text",x=11,y=0.7,label=annot)+ggtitle(assay)+scale_x_discrete(breaks = c(1,6,11,16,21),labels=labelBinCenter(binCount)[c(1,6,11,16,21)])+xlab(paste0("Genome bins surrounding sites (",total_width/1000,"kb)"))+ylab("DNA methylation")
  print(pl)
}
dev.off()


#collect dip scores
combined_dipScores=binResults[[1]]$dipScores[,c("id","surgery","sub_group","sub_group_prob","category","IDH"),with=FALSE]

for(assay in names(binResults)){
  print(assay)
dipScores=binResults[[assay]]$dipScores[,c("id","dipScore"),with=FALSE]
setnames(dipScores,"dipScore",assay)

if (length(combined_dipScores)==0){combined_dipScores=dipScores}else{
  combined_dipScores=merge(combined_dipScores,dipScores,by="id")
}
}

sort(names(combined_dipScores))
select=c("CTCF__Gliobla_None_62750","CTCF__H1-hESC_None_66533","CTCF__NH-A_None_38465","CHD1_(A301-218A)__H1-hESC_None_7245","NRSF__H1-hESC_None_13281","NRSF__U87_None_11741","NANOG__Embryonic Stem Cell_NA_18532","SOX2__Embryonic Stem Cell_NA_5204","Pol2__Gliobla_None_17405","Pol2__H1-hESC_None_20311","EZH2_(39875)__NH-A_None_6488","EZH2__Embryonic Stem Cell_NA_3204","POU5F1__Embryonic Stem Cell_NA_10374","POU5F1_(SC-9081)__H1-hESC_None_3996","RBBP5_(A300-109A)__H1-hESC_None_16143","POLR2A__Embryonic Stem Cell_NA_14498","KDM4A__Embryonic Stem Cell_NA_20621")

write.table(combined_dipScores[,c("id",select),with=FALSE],paste0("dipScores_",factor,"_sel.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

