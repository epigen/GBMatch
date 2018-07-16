#NOTE: This script produces DNA methykation tracks for selected regions of interest
library(project.init)
project.init2("GBMatch")


##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/02-meth_overview/")
dir.create(out_dir)
setwd(out_dir)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

#get enhancer data to select interesting regions (only needed in discovery phase)
#simpleCache("rrbsEnhancers")
#rrbsEnhancers_regions=rrbsEnhancers[,list(samples_covered=.N,mean_readCount=mean(readCount),mean_CpGcount=mean(CpGcount),mean_methyl=mean(methyl),sd_methyl=sd(methyl)),by=c("regionID","chr","start","end")]
#rrbsEnhancers_regions[samples_covered>500&mean_readCount>200&mean_CpGcount>40][order(sd_methyl)]


#focus on single CpGs (not tiles) 
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
####Figure 1c and S2a
#regions to plot
POIs=list(enhancer51330=c("chr22", 50298971,50304771),enhancer88653=c("chr7", 151114913,151115513),enhancer84097=c("chr7", 1667764,1669964),enhancer49388=c("chr21", 45016885, 45020285),PLXNB2=c("chr22",50274644,50309283),BCL2L11=c("chr2",111113953,111169264),SFRP2=c("chr4",153780540,153794943),MGMT=c("chr10",129445717,129773284),TERT=c("chr5",1248767,1302098))

#now plot
data1c=data.table()
for (POIname in names(POIs)){
  POI=POIs[[POIname]]
  print(POI)
  sel_chr=POI[1]
  lower=as.numeric(POI[2])
  upper=as.numeric(POI[3])
  
  sub=rrbsCg_annot[regions=="CG"&surgery.x%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&chr==sel_chr&start>lower&end<upper]
  sub[,id:=factor(id,levels=unique(id[order(surgery.x)])),]
  pdf(paste0("methylation_track_CG_",POIname,"_",sel_chr,"_",lower,".pdf"),width=5,height=4)
  pl=ggplot(sub)+geom_tile(width=0.005*(upper-lower), aes(fill=methyl,x=start,y=id))+xlab(paste0(sel_chr," (M)"))+scale_fill_gradient(low="blue",high="red")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+geom_vline(xintercept = c(lower,upper),col="lightgrey")+theme(legend.position="bottom")+facet_grid(category+surgery.x~.,scales="free_y",space="free_y")+xlim(c(lower,upper))
  print(pl)
  dev.off()
  data1c=rbindlist(list(data1c,sub[,region:=POIname,]))
}

data1c=data1c[region%in%c("SFRP2","PLXNB2","MGMT"),list(surgery.x,methyl,start,id,category,region)]
write.table(data1c,"Source Data Figure1c.csv",sep=";",quote=FALSE,row.names=FALSE)

#promoter zoom in 
####Figure S2b
POIs=list(MGMT_prom=c("chr10",129466600,129467700),TERT_prom=c("chr5",1294500,1296000),SFRP2_prom=c("chr4",153788500,153793000))

#now plot
for (POIname in names(POIs)){
  POI=POIs[[POIname]]
  print(POI)
  sel_chr=POI[1]
  lower=as.numeric(POI[2])
  upper=as.numeric(POI[3])
  
  sub=rrbsCg_annot[regions=="CG"&surgery.x%in%c(1,2)&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&chr==sel_chr&start>lower&end<upper]
  sub[,id:=factor(id,levels=unique(id[order(surgery.x)])),]
  pdf(paste0("methylation_track_CG_",POIname,"_",sel_chr,"_",lower,".pdf"),width=3,height=4)
  pl=ggplot(sub)+geom_tile(width=0.005*(upper-lower), aes(fill=methyl,x=start,y=id))+xlab(paste0(sel_chr," (M)"))+scale_fill_gradient(low="blue",high="red")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+geom_vline(xintercept = c(lower,upper),col="lightgrey")+theme(legend.position="bottom")+facet_grid(category+surgery.x~.,scales="free_y",space="free_y")+scale_x_continuous(breaks=c(lower,upper))
  print(pl)
  dev.off()
}