library(project.init)
project.init2("GBMatch")

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/02-meth_overview")
dir.create(out_dir)
setwd(out_dir)

annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))


###################MGMT methylation#####################################################################
loadCaches("rrbsCgNoMin")
setnames(rrbsCgNoMin,"id","sample")

#use one of the TCGA beta files to get the coordinates of probes on hg38
array_loc=fread(list.files(file.path(extData,"extData/TCGA/"),pattern=paste0(".*hg38.txt"),recursive = TRUE,full.names=TRUE)[1],select = c(1,3,4,5))
setnames(array_loc,names(array_loc),c("id","chr","start","end"))

#positions:
#cg12981137:chr10  129467311
#cg12434587: chr10  129466945

MGMT_cpgs=c("cg12981137", "cg12434587")
MGMT_loc=array_loc[id%in%MGMT_cpgs]


MGMT_meth_RRBS=rrbsCgNoMin[chr%in%MGMT_loc$chr&start%in%MGMT_loc$start]
sel="_all"

#calculate score
MGMT_meth_RRBS_combined=MGMT_meth_RRBS[,list(methyl=mean(methyl),meth_max=max(methyl),meth_min=min(methyl),readCount=sum(readCount),CpGcount=.N),by="sample"]
MGMT_meth_RRBS_combined=merge(MGMT_meth_RRBS_combined,setnames(annotation[,c("patID","IDH","category","surgery.x","N_number_seq"),with=FALSE],"N_number_seq","sample"),by="sample")


MGMT_meth_RRBS_combined[,mgmt_conf1:=ifelse(readCount>5,ifelse(methyl==1,"meth","unmeth"),NA),]
MGMT_meth_RRBS_combined[,mgmt_conf2:=ifelse(readCount>5,ifelse(methyl>0.5,"meth","unmeth"),NA),]
MGMT_meth_RRBS_combined[,mgmt_conf3:=ifelse(readCount>5,ifelse(methyl>0.36,"meth","unmeth"),NA),]
MGMT_meth_RRBS_combined[,mgmt_conf4:=ifelse(methyl>0.36,"meth","unmeth"),]

setnames(MGMT_meth_RRBS_combined,c("methyl","readCount","CpGcount"),c("mgmt_methyl","mgmt_readCount","mgmt_CpGcount"))

write.table(MGMT_meth_RRBS_combined,file.path(paste0("mgmt_status",sel,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
#MGMT_meth_RRBS_combined=fread(file.path(paste0("mgmt_status",sel,".tsv")))


MGMT_meth_RRBS_combined_red=MGMT_meth_RRBS_combined[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"]

MGMT_meth_RRBS_combined_red[mgmt_readCount>40,mgmt_readCount:=40,]
MGMT_meth_RRBS_combined_red[,id:=paste0(patID,"_",surgery.x),]

MGMT_meth_RRBS_combined_red[,sample:=factor(sample,levels=sample[order(mgmt_methyl,decreasing=TRUE)]),]
MGMT_meth_RRBS_combined_red[,id:=factor(id,levels=unique(id[order(mgmt_methyl,decreasing=TRUE)])),]


pdf(file.path(paste0("mgmt_status_GBMatch_",sel,".pdf")),height=7,width=18)
ggplot(MGMT_meth_RRBS_combined_red,aes(x=id,y=mgmt_methyl,alpha=mgmt_readCount,fill=as.factor(mgmt_CpGcount)))+geom_hline(yintercept = 0.36,lty=20)+geom_bar(stat="identity")+geom_point(pch=22,size=3)+theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust=1))+scale_alpha_continuous(range=c(0.4,1))+ geom_errorbar(aes(ymax=meth_max,ymin=meth_min), width=0.5,alpha=1)+scale_fill_manual(values=c("1"="#a8ddb5","2"="#43a2ca"))+facet_wrap(~category,ncol=1,scales="free_x")
dev.off()

pdf(file.path(paste0("mgmt_status_GBMatch_multi",sel,".pdf")),height=2.5,width=15)
ggplot(MGMT_meth_RRBS_combined_red[patID%in%patID[duplicated(patID)]],aes(x=as.factor(surgery.x),y=mgmt_methyl,alpha=mgmt_readCount,fill=as.factor(mgmt_CpGcount)))+geom_bar(stat="identity")+geom_point(pch=22,size=3)+geom_hline(yintercept = 0.36,lty=20)+theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust=1))+scale_alpha_continuous(range=c(0.4,1))+ geom_errorbar(aes(ymax=meth_max,ymin=meth_min), width=0.5,alpha=1)+scale_fill_manual(values=c("1"="#a8ddb5","2"="#43a2ca"))+facet_grid(~patID,space="free",scale="free")+ylab("MGMT methylation")
dev.off()       
