library(project.init)
project.init2("GBMatch")


##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/02-meth_overview/")
dir.create(out_dir)
setwd(out_dir)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

#only single CpGs (separately, because so big)
simpleCache("rrbsCg")
rrbsCg[,regions:="CG",]
rrbsCg[,regionID:=paste0(chr,"_",start),]
rrbsCg[,end:=start+1,]

#annotate
annot_sub=annotation[,c("N_number_seq","patID","category","material","Center","Siteofsurgery","surgery.x","qualTier","enrichmentCycles","IDH"),with=FALSE]
setnames(annot_sub,"N_number_seq","id")

rrbsCg_annot=merge(annot_sub,rrbsCg,by="id")
rrbsCg_annot[,methyl:=methyl*100,]
rrbsCg_annot[,regionID_ext:=paste0(regions,"_",regionID),]
rm(rrbsCg)

#plot tracks 

#regions to plot
POIs=list(BCL2L11=c("chr2",111113953,111169264),SFRP2=c("chr4",153780540,153794943),MGMT=c("chr10",129445717,129773284))

#now plot
for (POI in POIs){
  print(POI)
  sel_chr=POI[1]
  lower=as.numeric(POI[2])
  upper=as.numeric(POI[3])
  
  sub=rrbsCg_annot[regions=="CG"&surgery.x%in%c(1,2)&category%in%c("GBMatch")&IDH=="wt"&chr==sel_chr&start>lower&end<upper]
  sub[,id:=factor(id,levels=unique(id[order(surgery.x)])),]
  pdf(paste0("methylation_track_CG_",sel_chr,"_",lower,".pdf"),width=5,height=5)
    pl=ggplot(sub)+geom_tile(width=0.005*(upper-lower), aes(fill=methyl,x=start,y=id))+xlab(paste0(sel_chr," (M)"))+scale_fill_gradient(low="blue",high="red")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+theme(legend.position="bottom")+facet_wrap(~surgery.x,ncol=1,scale="free")+xlim(c(lower,upper))
    print(pl)
  dev.off()
}
