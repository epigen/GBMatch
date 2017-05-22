library(project.init)
project.init2("GBMatch")

#annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

annot_ASC=annotation[IDH=="wt"&category=="GBMatch",list(cycles.1=enrichmentCycles[surgery.x==1][1],cycles.2=enrichmentCycles[surgery.x==2][1],Age=min(Age),Sex=unique(Sex), extentOfResection=Extentofresection[surgery.x==1][1],progression_location=progression_location[surgery.x==2][1],timeToSecSurg=unique(timeToSecSurg),category=unique(category),IDH=unique(IDH)),by=patID]
setnames(annot_ASC,"patID","patient")

##set it up
min_reads=20  #20 or 60 (use for furtehr analyses)


out_dir=file.path(getOption("PROCESSED.PROJECT"),paste0("results_analysis/06-methclone/summary_minReads",min_reads))
dir.create(out_dir)
setwd(out_dir)
methclone_files=list.files(path=paste0("../data_minReads",min_reads),pattern="*.output.txt",full.names=TRUE,recursive=TRUE)

##get the data
simpleCache("rbindlist(sapply(methclone_files,function(x){fread(x,select=c(1:48))},simplify=FALSE))",cacheName=paste0("methclone_all_min",min_reads),cacheSubDir="entropy",assignToVariable="methclone")
#for some reason Methclone seldomly gives two result lines for exactly the same location --> only keep location with lower entropy
methclone[,keep:=min(entropy),by=c("chr","start","end","strand","sample")]
methclone=methclone[keep==entropy]
spl=unlist(strsplit(methclone$sample,"__|vs"))
methclone[,patient:=spl[seq(1,to=length(spl),by=5)],]
methclone[,sample_1:=spl[seq(2,to=length(spl),by=5)],]
methclone[,sample_2:=spl[seq(3,to=length(spl),by=5)],]
methclone[,timepoint_1:=spl[seq(4,to=length(spl),by=5)],]
methclone[,timepoint_2:=spl[seq(5,to=length(spl),by=5)],]
methclone[,comparison:=paste0(timepoint_1,"vs",timepoint_2),]
methclone[,comparison_simpl:=gsub("ctr_.","0",comparison),]
methclone_self=methclone[sample_1==sample_2]
methclone=methclone[(sample_1!=sample_2)]

##################################
##Sample centered (self vs self)##
##################################
methclone_red=methclone_self[,c("chr","start","end","entropy1","read1","meth1","sample_1",grep("s0",names(methclone_self),value=TRUE)),with=FALSE]
setnames(methclone_red,c("sample_1",grep("s0",names(methclone_self),value=TRUE)),c("id",gsub("s0:","p",grep("s0",names(methclone_self),value=TRUE))))
simpleCache(paste0("rrbs_entropy_min",min_reads),methclone_red,recreate=FALSE)

#Aggregate entropy over regions
eload(loadTiles(genomeBuild=genome, tileSize=1000))
eload(loadTiles(genomeBuild=genome, tileSize=5000))
eload(loadGencodeGenes("human",versNum=87))
prom1k=promoters(SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]], upstream=1000, downstream=500)
prom1k=prom1k[!duplicated(prom1k)]

tiled_entropy=function(methclone_dt,tiles,exclude=NULL){
  # Build summary J command
  cols=c("entropy1","read1","meth1",grep("p[01]",names(methclone_dt),value=TRUE))
  funcs = c("mean", "sum","mean",rep("mean",16))
  jCommand = buildJ(cols,funcs)
  rrbsAgentropy = BSAggregate(methclone_dt, tiles ,excludeGR=exclude, jCommand=jCommand, splitFactor="id")
  return(rrbsAgentropy)
}

simpleCache(paste0("rrbsAgEntropy_min",min_reads,"_1k"),tiled_entropy(get(paste0("rrbs_entropy_min",min_reads)),SV$tiles1000hg38),cacheSubDir="entropy",recreate=FALSE)
simpleCache(paste0("rrbsAgEntropy_min",min_reads,"_5k"),tiled_entropy(get(paste0("rrbs_entropy_min",min_reads)),SV$tiles5000hg38),cacheSubDir="entropy",recreate=FALSE)
simpleCache(paste0("rrbsAgEntropy_min",min_reads,"_prom1k"),tiled_entropy(get(paste0("rrbs_entropy_min",min_reads)),prom1k),cacheSubDir="entropy",recreate=FALSE)

#now get rid of the frozen samples
methclone=methclone[!(sample_1%in%annotation[category=="GLASS"]$N_number_seq|sample_2%in%annotation[category=="GLASS"]$N_number_seq)]

##########################################
##patient centered eloci/ shift analyses##
##########################################

EPM_tab=methclone[,list(eloci=sum(entropy<(-80)),all=sum(!is.na(entropy)),sd_dentropy=sd(entropy),mean_dentropy=mean(entropy),mean_entropy1=mean(entropy1), mean_entropy2=mean(entropy2),sd_entropy1=sd(entropy1), sd_entropy2=sd(entropy2)),by=c("sample","patient",  "sample_1", "sample_2", "timepoint_1", "timepoint_2","comparison","comparison_simpl")]
EPM_tab[,EPM:=eloci/all*1000000,]

EPM_tab[,decreasing:=median(EPM[grepl("ctr.*vs1",comparison)])>median(EPM[grepl("ctr.*vs2",comparison)]),by="patient"]

EPM_tab=merge(EPM_tab,annot_ASC,by="patient",all.x=TRUE)
#best control: N912-12F = ctr_5 -->remove all others
EPM_tab=EPM_tab[!timepoint_1%in%c("ctr_1","ctr_2","ctr_3","ctr_4")&!timepoint_1%in%c("ctr_1","ctr_2","ctr_3","ctr_4")]

#Entropy per patient 
write.table(EPM_tab[comparison_simpl=="1vs2"],"EPM_1vs2.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(EPM_tab[comparison_simpl=="0vs2"&timepoint_1=="ctr_5"],"EPM_0vs2.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(EPM_tab[comparison_simpl=="0vs1"&timepoint_1=="ctr_5"],"EPM_0vs1.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#check correlation between enrichment cycles and EPM
pdf("EnrichmentCycles_cor.pdf",height=5,width=5)
ggplot(EPM_tab[category=="GBMatch"],aes(y=mean_entropy1,x=cycles.1))+geom_point(alpha=0.6,position=position_jitter(width=0.3,height=0))
ggplot(EPM_tab[category=="GBMatch"],aes(y=mean_entropy2,x=cycles.2))+geom_point(alpha=0.6,position=position_jitter(width=0.3,height=0))
dev.off()


sub=EPM_tab[IDH=="wt"&category=="GBMatch"&timepoint_2<3&timepoint_1!=timepoint_2&cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12]

sub[,comparison:=factor(comparison,levels=c("ctr_5vs1","ctr_5vs2","1vs2")),]
sub=sub[!is.na(comparison)]

sign_1=sub[,wilcox.test(x=EPM[comparison=="ctr_5vs2"],y=EPM[comparison=="ctr_5vs1"])$p.value,]
sign_2=sub[,wilcox.test(x=EPM[comparison=="ctr_5vs1"],y=EPM[comparison=="1vs2"])$p.value,]
sign_3=sub[,wilcox.test(x=EPM[comparison=="ctr_5vs2"],y=EPM[comparison=="1vs2"])$p.value,]

pdf("EPM_boxplot.pdf",height=3.5,width=2)
ggplot(sub,aes(y=log10(EPM+1),x=comparison))+geom_point(size=2,position=position_jitter(width=0.3),alpha=0.6,shape=21,aes(fill=comparison,col=comparison))+geom_boxplot(outlier.shape=NA,fill="transparent",col="black")+theme(legend.position='none')+annotate("text",x=c(1.5,1,2),y=5.9,label=signif(c(sign_1,sign_2,sign_3),3))+scale_fill_manual(values=c("ctr_5vs1"="#bdc9e1","ctr_5vs2"="#9ecae1","1vs2"="#31a354"))+scale_colour_manual(values=c("ctr_5vs1"="#bdc9e1","ctr_5vs2"="#9ecae1","1vs2"="#31a354"))
dev.off()

sub[,patient:=factor(patient,levels=unique(patient[comparison=="1vs2"][order(EPM[comparison=="1vs2"])])),]
pdf("EPM_barplots.pdf",height=5,width=7)
ggplot(sub[!is.na(comparison)],aes(y=log10(EPM+1),x=patient))+geom_bar(stat="identity",aes(fill=comparison))+facet_wrap(~comparison,ncol=1)+scale_fill_manual(values=c("ctr_5vs1"="#bdc9e1","ctr_5vs2"="#9ecae1","1vs2"="#31a354"))+theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1))
dev.off()


#EPM comparisons
sub_cyc=sub
sub_cyc[,log10_EPM:=log10(EPM),]
cors=sub_cyc[,list(cor=cor(timeToSecSurg,log10_EPM),y=max(log10_EPM)),by="comparison"]

pdf("EPM_comps.pdf",height=3.5,width=7)
ggplot(sub_cyc,aes(y=log10_EPM,x=timeToSecSurg/30))+geom_point(shape=21)+geom_smooth(method="lm")+geom_text(data=cors,x=10,aes(y=y,label=paste0("r=",round(cor,3))))+facet_wrap(~comparison,scale="free")+xlab("timeToSecSurg (months)")

ggplot(sub_cyc,aes(y=log10_EPM,x=as.factor(progression_location),group=progression_location))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(fill="transparent",outlier.shape=NA)+facet_wrap(~comparison,scale="free")

ggplot(sub_cyc,aes(y=log10_EPM,x=as.factor(extentOfResection),group=extentOfResection))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(outlier.shape=NA,fill="transparent")+facet_wrap(~comparison)
dev.off()

pdf("EPM_comps_1vs2.pdf",height=3.5,width=3)
ggplot(sub_cyc[comparison=="1vs2"],aes(y=log10_EPM,x=timeToSecSurg/30))+geom_point(shape=21,size=3,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cors[comparison=="1vs2"],x=10,aes(y=y,label=paste0("r=",round(cor,3))))+xlab("timeToSecSurg (months)")

ggplot(sub_cyc[comparison=="1vs2"],aes(y=log10_EPM,x=as.factor(progression_location),group=progression_location))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(fill="transparent",outlier.shape=NA)

ggplot(sub_cyc[comparison=="1vs2"],aes(y=log10_EPM,x=as.factor(extentOfResection),group=extentOfResection))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(outlier.shape=NA,fill="transparent")
dev.off()




