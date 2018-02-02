library(project.init)
project.init2("GBMatch")
library(LOLA)
library(RGenomeUtils)
library(MASS)
library(gridExtra)


##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/02-meth_overview/")
dir.create(out_dir)
setwd(out_dir)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))


#overview in tiles
#simpleCache("rrbsTiled1ksub")
#rrbsTiled1ksub[,regions:="tiled_1kb",]
simpleCache("rrbsTiled5ksub")
rrbsTiled5ksub[,regions:="tiled_5kb",]
simpleCache("rrbsProm1kb")
rrbsProm1kb[,regions:="prom_1kb",]
simpleCache("rrbsCGI")
rrbsCGI[,regions:="CGI",]
simpleCache("rrbsEnhancers")
rrbsEnhancers[,regions:="Enhancers",]


#assess batch effects by MDS
rrbsTiled5ksub_wide=as.data.table(dcast(rrbsTiled5ksub,regionID~id,value.var="methyl"))
rrbsTiled5ksub_wide_mat=as.matrix(rrbsTiled5ksub_wide[,-c("regionID"),with=FALSE])
row.names(rrbsTiled5ksub_wide_mat)=rrbsTiled5ksub_wide$regionID

simpleCache(cacheName="rrbsTiled5ksub_dist",instruction="dist(t(rrbsTiled5ksub_wide_mat))",cacheDir=getwd(),assignToVariable="d")
##save(d,file="rrbsTiled5ksub_dist.RData") ##only needed first time when run without simpleCache
mds=isoMDS(d, k=2)
mds_annot_pre=data.table(N_number_seq=row.names(mds$points),MDS1=mds$points[,1],MDS2=mds$points[,2])
mds_annot=merge(mds_annot_pre,annotation,by="N_number_seq")

batch_variables=c("SurgeryDate","material","category","cohort","batch","adapter","flowcell.y")

mds_annot[,SurgeryDate:=as.Date(SurgeryDate,format="%Y-%m-%d"),]
mds_annot[,batch:=as.factor(batch),]


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


pdf("MDS_batch_color.pdf",height=6,width=4)
for (batch_variable in batch_variables){
  pl=ggplot(mds_annot,aes(x=MDS1,y=MDS2,fill=get(batch_variable),col=get(batch_variable)))+geom_point(alpha=0.6)+ggtitle(batch_variable)
  print(pl+ theme(legend.position="none")+ coord_fixed())  
  grid.arrange(g_legend(pl))
}
dev.off()


#actually combine
meth_combined=rbindlist(list(rrbsTiled5ksub,rrbsProm1kb,rrbsCGI,rrbsEnhancers))
rm(list=c("rrbsTiled5ksub","rrbsProm1kb","rrbsCGI","rrbsEnhancers"))


#annotate
annot_sub=annotation[,c("N_number_seq","patID","category","cohort","material","Center","Siteofsurgery","surgery.x","qualTier","enrichmentCycles","IDH"),with=FALSE]
setnames(annot_sub,"N_number_seq","id")

meth_combined_annot=merge(annot_sub,meth_combined,by="id")
meth_combined_annot[,methyl:=methyl*100,]
meth_combined_annot[,regionID_ext:=paste0(regions,"_",regionID),]
rm(meth_combined)

#determine core regions
sample_count=length(unique(meth_combined_annot$id))
core_regions=meth_combined_annot[,list(select=ifelse(sum(!is.na(methyl))>0.75*sample_count,TRUE,FALSE)),by=c( "chr","start", "end", "regionID_ext","regions")]
core_regions[,sum(select),by="regions"]


# check out correlations with enrichment cycles separately for each region
sub=meth_combined_annot[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"]
sub[,cycles_cor:=cor(x=methyl,y=enrichmentCycles),by=c("regionID_ext","regions","category")]
sub[,cycles_cor_class:=cut(cycles_cor,10),]

sub_unique=sub[,list(readCount=mean(readCount), methyl=mean(methyl), CpGcount=mean(CpGcount),cycles_cor=unique(cycles_cor) ,cycles_cor_class=unique(cycles_cor_class) ),by=c("regionID_ext","regions","category")]

pdf("methylation_cycle_correlation.pdf",width=8,height=3)
ggplot(sub_unique[readCount/CpGcount>10],aes(x=cycles_cor,col=regions))+geom_density()+ggtitle(">10 reads per CpG")+facet_wrap(~category)
ggplot(sub_unique[CpGcount>25],aes(x=cycles_cor,col=regions))+geom_density()+ggtitle(">25 CpGs")+facet_wrap(~category)
ggplot(sub_unique[readCount/CpGcount>50],aes(x=cycles_cor,col=regions))+geom_density()+ggtitle(">50 reads per CpG")+facet_wrap(~category)
dev.off()


#now plot overview

#for primary cohort
pdf("methylation_overview_material_prim.pdf",width=5,height=5)
ggplot(meth_combined_annot[(readCount/CpGcount)>10&cohort=="primary"],aes(x=material,y=methyl,group=material))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+xlab("")+ggtitle("> 10 reads per CpG")

ggplot(meth_combined_annot[CpGcount>25&cohort=="primary"],aes(x=material,y=methyl,group=material))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+xlab("")+ggtitle("> 10 reads per CpG")+ggtitle(">25 CpGs per region")
dev.off()

pdf("methylation_overview_qualTier_prim.pdf",width=6,height=5)
ggplot(meth_combined_annot[(readCount/CpGcount)>10&cohort=="primary"],aes(y=methyl,x=as.factor(qualTier),group=as.factor(qualTier)))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+ggtitle("> 10 reads per CpG")

ggplot(meth_combined_annot[CpGcount>25&cohort=="primary"],aes(y=methyl,x=as.factor(qualTier),group=as.factor(qualTier)))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+ggtitle(">25 CpGs per region")
dev.off()


#for validation cohort
pdf("methylation_overview_material_val.pdf",width=4,height=5)
ggplot(meth_combined_annot[(readCount/CpGcount)>10&cohort=="validation"],aes(x=material,y=methyl,group=material))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+xlab("")+ggtitle("> 10 reads per CpG")

ggplot(meth_combined_annot[CpGcount>25&cohort=="validation"],aes(x=material,y=methyl,group=material))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+xlab("")+ggtitle("> 10 reads per CpG")+ggtitle(">25 CpGs per region")
dev.off()

pdf("methylation_overview_qualTier_val.pdf",width=6,height=5)
ggplot(meth_combined_annot[(readCount/CpGcount)>10&cohort=="validation"],aes(y=methyl,x=as.factor(qualTier),group=as.factor(qualTier)))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+ggtitle("> 10 reads per CpG")

ggplot(meth_combined_annot[CpGcount>25&cohort=="validation"],aes(y=methyl,x=as.factor(qualTier),group=as.factor(qualTier)))+geom_violin(fill="black")+geom_boxplot(outlier.size=NA,width=0.2)+facet_wrap(~regions,scale="free")+ylab("% DNA methylation")+ggtitle(">25 CpGs per region")
dev.off()




#scatter plot (first vs. second samples)
sel_samp=meth_combined_annot[regions=="tiled_5kb"&surgery.x%in%c(1,2)&(readCount/CpGcount)>10&CpGcount>25]
sel_samp_wide=as.data.table(dcast(sel_samp,regionID_ext+patID+material~surgery.x,value.var="methyl",fun.aggregate=function(x){x[1]}))
sel_samp_wide[,plotID:=paste0(patID,"_",material),]
cors_samples=sel_samp_wide[!is.na(`1`)&!is.na(`2`),list(corelation=cor(x=`1`,y=`2`)),by=c("patID","material","plotID")]
cors_material=sel_samp_wide[!is.na(`1`)&!is.na(`2`),list(corelation=cor(x=`1`,y=`2`)),by=c("material")]

pdf("methylation_cor1vs2_material.pdf",width=10,height=3)

ggplot(data=sel_samp_wide,aes(x=`1`,y=`2`))+stat_bin_hex(aes(fill=..density..^0.15))+geom_text(data=cors_material,x=15,y=98,aes(label=paste0("r=",round(corelation,3))))+scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+facet_wrap(~material,scale="free")+xlab("Primary")+ylab("Recurring")+theme(aspect.ratio=1,panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

pdf("methylation_cor1vs2_patients.pdf",width=4.3,height=3)
for(patient in unique(sel_samp_wide$plotID)){
  print(patient)
  
  pl=ggplot(data=sel_samp_wide[plotID==patient],aes(x=`1`,y=`2`))+stat_bin_hex(bins=30,aes(fill=..density..^0.01))+annotate("text",x=10,y=98,label=paste0("r=",round(cors_samples[plotID==patient]$corelation,3)))+scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+ggtitle(patient)+xlab("Primary")+ylab("Recurring")+theme(aspect.ratio=1,panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(pl)
}
dev.off()

