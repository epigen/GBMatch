#NOTE: this script analyzes the methclone (DNA methyaltion entropy) results 
#(3 different read cutoffs and 6 different entropy shift cutoffs)
library(project.init)
project.init2("GBMatch")
library(ChIPpeakAnno)
library(org.Hs.eg.db)
data(TSS.human.GRCh38)
library(LOLA)

#annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

annot_ASC=annotation[IDH=="wt"&category=="GBMatch",list(cycles.1=enrichmentCycles[surgery.x==1][1],cycles.2=enrichmentCycles[surgery.x==2][1],Age=min(Age),Sex=unique(Sex), extentOfResection=Extentofresection[surgery.x==1][1],progression_location=progression_location[surgery.x==2][1],timeToSecSurg=unique(timeToSecSurg),category=unique(category),IDH=unique(IDH)),by=patID]
setnames(annot_ASC,"patID","patient")

#load lola regiondb (needed below)
regionDB_core = loadRegionDB("/data/groups/lab_bock/shared/resources/regions/LOLACore/hg38/")
cellType_conversions=fread(file.path(getOption("PROJECT.DIR"),"metadata/LOLA_annot/CellTypes.tsv"),drop="collection")

#region annotation
eload(loadTiles(genomeBuild=genome, tileSize=1000))
eload(loadTiles(genomeBuild=genome, tileSize=5000))
eload(loadGencodeGenes("human",versNum=87))
prom1k=promoters(SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]], upstream=1000, downstream=500)
prom1k=prom1k[!duplicated(prom1k)]
min_reads_list=c(20,40,60) #40 used for further analysis


for (min_reads in min_reads_list){
  message(paste0("Min reads: ",min_reads))
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
  
  #make sure to only use samples in sample annotation sheet
  methclone=methclone[sample_1%in%annotation$N_number_seq&sample_2%in%annotation$N_number_seq]
  
  methclone_self=methclone[sample_1==sample_2]
  methclone=methclone[(sample_1!=sample_2)]
  
  ##################################
  ##Sample centered (self vs self)##
  ##################################
  methclone_red=methclone_self[,c("chr","start","end","entropy1","read1","meth1","sample_1",grep("s0",names(methclone_self),value=TRUE)),with=FALSE]
  setnames(methclone_red,c("sample_1",grep("s0",names(methclone_self),value=TRUE)),c("id",gsub("s0:","p",grep("s0",names(methclone_self),value=TRUE))))
  simpleCache(paste0("rrbs_entropy_min",min_reads),methclone_red,recreate=FALSE)
  
  #Aggregate entropy over regions
  tiled_entropy=function(methclone_dt,tiles,exclude=NULL){
    # Build summary J command
    cols=c("entropy1","read1","meth1",grep("p[01]",names(methclone_dt),value=TRUE))
    funcs = c("mean", "sum","mean",rep("mean",16))
    jCommand = buildJ(cols,funcs)
    rrbsAgentropy = BSAggregate(methclone_dt, tiles ,excludeGR=exclude, jCommand=jCommand, splitFactor="id")
    return(rrbsAgentropy)
  }
  
  simpleCache(paste0("rrbsAgEntropy_min",min_reads,"_1k"),tiled_entropy(get(paste0("rrbs_entropy_min",min_reads)),SV$tiles1000hg38),cacheSubDir="entropy",recreate=FALSE, noload=TRUE)
  simpleCache(paste0("rrbsAgEntropy_min",min_reads,"_5k"),tiled_entropy(get(paste0("rrbs_entropy_min",min_reads)),SV$tiles5000hg38),cacheSubDir="entropy",recreate=FALSE, noload=TRUE)
  simpleCache(paste0("rrbsAgEntropy_min",min_reads,"_prom1k"),tiled_entropy(get(paste0("rrbs_entropy_min",min_reads)),prom1k),cacheSubDir="entropy",recreate=FALSE, noload=TRUE)
  
  
  ##########################################
  ##patient centered eloci/ shift analyses##
  ##########################################
  #now get rid of the frozen samples
  methclone=methclone[!(sample_1%in%annotation[material=="frozen"]$N_number_seq|sample_2%in%annotation[material=="frozen"]$N_number_seq)]
  
  #annotate the methclone resuls with genetic features
  methclone_gr=with(methclone,GRanges(GRanges(seqnames = Rle(chr), IRanges(start=start, end=end),strand=Rle(strand),sample=sample,entropy=entropy)))
  #gene annotation
  anno=annotatePeakInBatch(methclone_gr, AnnotationData=TSS.human.GRCh38)
  anno=addGeneIDs(annotatedPeak=anno,orgAnn="org.Hs.eg.db", IDs2Add="symbol")
  anno_dt=as.data.table(as.data.frame(anno))
  setnames(anno_dt,"seqnames","chr")
  methclone_anno=merge(anno_dt,methclone,by=c("chr","start","end","entropy","sample","strand"))
  methclone_anno[,distancetoFeature_log10:=ifelse(distancetoFeature>0,log10(distancetoFeature+1),-log10(abs(distancetoFeature)+1)),]
  methclone_anno=merge(methclone_anno,annot_ASC,by="patient",all.x=TRUE)
  
  #create GRanges list for all comparisons (LOLA analysis)
  methclone_grl=GRangesList()
  tests=unique(methclone$sample)
  for (test in tests){
    methclone_grl[[test]]=with(methclone[sample==test],GRanges(seqnames = Rle(chr), IRanges(start=start, end=end),strand=Rle(strand),entropy=entropy)) 
  } 

  #try range of entropy cutoffs according to reviewer comment (-80 used for further analysis)
  entropy_cutoffs=c(-40,-50,-60, -70, -80, -90)
  for (entropy_cutoff in entropy_cutoffs){
    subdir=paste0("EPM_epy_cutoff",entropy_cutoff)
    dir.create(subdir)
    setwd(subdir)  
    
    EPM_tab=methclone[,list(eloci=sum(entropy<(entropy_cutoff)),all=sum(!is.na(entropy)),sd_dentropy=sd(entropy),mean_dentropy=mean(entropy),mean_entropy1=mean(entropy1), mean_entropy2=mean(entropy2),sd_entropy1=sd(entropy1), sd_entropy2=sd(entropy2)),by=c("sample","patient",  "sample_1", "sample_2", "timepoint_1", "timepoint_2","comparison","comparison_simpl")]
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
    print(ggplot(EPM_tab[category=="GBMatch"],aes(y=mean_entropy1,x=cycles.1))+geom_point(alpha=0.6,position=position_jitter(width=0.3,height=0)))
    print(ggplot(EPM_tab[category=="GBMatch"],aes(y=mean_entropy2,x=cycles.2))+geom_point(alpha=0.6,position=position_jitter(width=0.3,height=0)))
    dev.off()
    
    sub=EPM_tab[IDH=="wt"&category=="GBMatch"&timepoint_2<3&timepoint_1!=timepoint_2&cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12]
  
    sub[,comparison:=factor(comparison,levels=c("ctr_5vs1","ctr_5vs2","1vs2")),]
    sub=sub[!is.na(comparison)]
    
    sign_1=sub[,signif(wilcox.test(x=EPM[comparison=="ctr_5vs2"],y=EPM[comparison=="ctr_5vs1"])$p.value,3),]
    sign_2=sub[,signif(wilcox.test(x=EPM[comparison=="ctr_5vs1"],y=EPM[comparison=="1vs2"])$p.value,3),]
    sign_3=sub[,signif(wilcox.test(x=EPM[comparison=="ctr_5vs2"],y=EPM[comparison=="1vs2"])$p.value,3),]
    
    ####Figure S12d
    pdf("EPM_boxplot.pdf",height=3.5,width=2)
    print(ggplot(sub,aes(y=log10(EPM+1),x=comparison))+geom_point(size=2,position=position_jitter(width=0.3),alpha=0.6,shape=21,aes(fill=comparison,col=comparison))+geom_boxplot(outlier.shape=NA,fill="transparent",col="black")+theme(legend.position='none')+annotate("text",x=c(1.5,1,2),y=5.9,label=c(sign_1,sign_2,sign_3))+scale_fill_manual(values=c("ctr_5vs1"="#bdc9e1","ctr_5vs2"="#9ecae1","1vs2"="#31a354"))+scale_colour_manual(values=c("ctr_5vs1"="#bdc9e1","ctr_5vs2"="#9ecae1","1vs2"="#31a354"))+stat_summary(fun.data=addN, geom="text", vjust=-0.5, col="blue"))
    dev.off()
    
    sub[,patient:=factor(patient,levels=unique(patient[comparison=="1vs2"][order(EPM[comparison=="1vs2"])])),]
    pdf("EPM_barplots.pdf",height=5,width=7)
    print(ggplot(sub[!is.na(comparison)],aes(y=log10(EPM+1),x=patient))+geom_bar(stat="identity",aes(fill=comparison))+facet_wrap(~comparison,ncol=1)+scale_fill_manual(values=c("ctr_5vs1"="#bdc9e1","ctr_5vs2"="#9ecae1","1vs2"="#31a354"))+theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1)))
    dev.off()
    
    #EPM comparisons
    sub_cyc=sub
    sub_cyc[,log10_EPM:=log10(EPM),]
    cors=sub_cyc[,list(cor=cor.test(timeToSecSurg,log10_EPM)$estimate,p.value=cor.test(timeToSecSurg,log10_EPM)$p.value,N=length(timeToSecSurg),y=max(log10_EPM)),by="comparison"]
    
    pdf("EPM_comps.pdf",height=3.5,width=7)
    print(ggplot(sub_cyc,aes(y=log10_EPM,x=timeToSecSurg/30))+geom_point(shape=21)+geom_smooth(method="lm")+geom_text(data=cors,x=10,aes(y=y,label=paste0("r=",signif(cor,2)," p=",signif(p.value,2)," N=",N)),size=2.5)+facet_wrap(~comparison,scale="free")+xlab("timeToSecSurg (months)"))
    
    print(ggplot(sub_cyc,aes(y=log10_EPM,x=as.factor(progression_location),group=progression_location))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(fill="transparent",outlier.shape=NA)+facet_wrap(~comparison,scale="free"))
    
    print(ggplot(sub_cyc,aes(y=log10_EPM,x=as.factor(extentOfResection),group=extentOfResection))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(outlier.shape=NA,fill="transparent")+facet_wrap(~comparison))
    dev.off()
    
    ####Figure S12f
    pdf("EPM_comps_1vs2.pdf",height=3.5,width=3)
    print(ggplot(sub_cyc[comparison=="1vs2"],aes(y=log10_EPM,x=timeToSecSurg/30))+geom_point(shape=21,size=3,fill="grey",alpha=0.6)+geom_smooth(method="lm",fill="lightgrey")+geom_text(data=cors[comparison=="1vs2"],x=20,aes(y=y,label=paste0("r=",signif(cor,2)," p=",signif(p.value,2)," N=",N)),size=3)+xlab("timeToSecSurg (months)"))
    
    print(ggplot(sub_cyc[comparison=="1vs2"],aes(y=log10_EPM,x=as.factor(progression_location),group=progression_location))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(fill="transparent",outlier.shape=NA))
    
    print(ggplot(sub_cyc[comparison=="1vs2"],aes(y=log10_EPM,x=as.factor(extentOfResection),group=extentOfResection))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(outlier.shape=NA,fill="transparent"))
    dev.off()
    
    ####Figure S12g
    pdf("distance_loci_tss_1vs2.pdf",height=5,width=6)
    print(ggplot(methclone_anno[comparison_simpl%in%c("1vs2","0vs0")&((cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12)|comparison_simpl%in%c("0vs0"))],aes(x=distancetoFeature_log10,linetype=entropy<(entropy_cutoff),col=comparison_simpl=="0vs0"))+geom_density()+scale_color_discrete(name="Control")+scale_linetype_discrete(name="Eloci"))
    dev.off()
    
    ##run LOLA
    #first run LOLA on eloci in any 1vs2 comparison
    uset_all=unique(subset(methclone_gr,entropy<(entropy_cutoff)&grepl("__1vs2",sample)&sample%in%methclone_anno[(cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12)]$sample))
    univ_all=unique(subset(methclone_gr,grepl("__1vs2",sample)&sample%in%methclone_anno[(cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12)]$sample))
    N_samples=length(unique(uset_all$sample))
    N_regions=length(uset_all)
    
    simpleCache(cacheName=paste0("allEloci_",min_reads,"epy",abs(entropy_cutoff)),{df=runLOLA(userSets=GRangesList(allEloci__all1vsall2__1vs2all=uset_all),userUniverse=univ_all,regionDB=regionDB_core);return(df)},cacheSubDir="methcloneLOLA",recreate=FALSE)
    
    locResults=get(paste0("allEloci_",min_reads,"epy",abs(entropy_cutoff)))
    locResults[,sample:="allEloci__all1vsall2__1vs2all",]
    
    #then run for all patients/contols separately
    for(list in c(grep("__ctr_[1-5]vsctr_[1-5]",names(methclone_grl),value=TRUE),grep("__1vs2",names(methclone_grl),value=TRUE))){
      uset=subset(methclone_grl[[list]],entropy<(entropy_cutoff))
      univ=methclone_grl[[list]]
      if (length(uset)<5|length(univ)<10){
        print("Too fiew loci. Skipping!")
        next
      }
      simpleCache(cacheName=paste0(list,"_",min_reads,"epy",abs(entropy_cutoff)),{df=runLOLA(userSets=uset,userUniverse=univ,regionDB=regionDB_core);return(df)},cacheSubDir="methcloneLOLA")
      locResults=rbindlist(list(locResults,get(paste0(list,"_",min_reads,"epy",abs(entropy_cutoff)))[,sample:=list]))
      rm(list=c(paste0(list,"_",min_reads,"epy",abs(entropy_cutoff))))
    }
    
    collections=c("codex","encode_tfbs")
    locResults=locResults[collection%in%collections]
    locResults[,p.adjust:=p.adjust(10^(-pValueLog),method="BY")]
    locResults[,mlog10p.adjust:=-log10(p.adjust),]
    spl=unlist(strsplit(locResults$sample,"__|vs"))
    locResults[,patient:=spl[seq(1,to=length(spl),by=5)],]
    locResults[,sample_1:=spl[seq(2,to=length(spl),by=5)],]
    locResults[,sample_2:=spl[seq(3,to=length(spl),by=5)],]
    locResults[,timepoint_1:=spl[seq(4,to=length(spl),by=5)],]
    locResults[,timepoint_2:=spl[seq(5,to=length(spl),by=5)],]
    locResults[,comparison:=paste0(timepoint_1,"vs",timepoint_2),]
    locResults[,comparison_simpl:=gsub("ctr_.","0",comparison),]
    
    locResults[cellType==""|is.na(cellType),cellType:="Not defined",]
    locResults=merge(locResults,cellType_conversions,by="cellType",all=TRUE)
    locResults[description=="T-cell acute lymphoblastic leukaemia (T-ALL) cell line.",c("Lineage1","Lineage","cellType_corr"):=list(Lineage1="Lymphoid",Lineage="Lymphoid",cellType_corr="T lymphocyte"),]
    locResults[,target:=toupper(sub("-","",unlist(lapply(antibody,function(x){spl=unlist(strsplit(x,"_|eGFP-"));spl[spl!=""][1]})))),]
    locResults=locResults[!is.na(userSet)]
    locResults_red=locResults[,list(p.adjust=min(p.adjust),logOddsRatio=logOddsRatio[which.min(p.adjust)]),by=c("sample","patient","cellType_corr","target","comparison_simpl")]
    locResults_red[,dbSet:=paste0(target,":",cellType_corr),]
    locResults_red_annot=merge(locResults_red,annot_ASC,by="patient",all.x=TRUE)
    write.table(locResults_red_annot,"LOLA_locResults_condensed.tsv",quote=FALSE,sep="\t",row.names=FALSE)
    
    sub=locResults_red_annot[dbSet%in%locResults_red_annot[p.adjust<0.001]$dbSet&((cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12)|comparison_simpl%in%c("0vs0","1vs2all"))]
    sub[,p.adjust:=ifelse(p.adjust==0,min(p.adjust[p.adjust>0]),p.adjust)]
    sub[,N_comparisons:=length(unique(sample)),by="comparison_simpl"]
    sub[,facet_label:=paste0(comparison_simpl,"\n",N_comparisons," comparisons"),]
    
    pdf("LOLA_signif_all.pdf",height=7,width=9)
    print(ggplot(sub,aes(x=dbSet,y=-log10(p.adjust),col=target))+geom_boxplot()+geom_abline(intercept=-log10(0.001),slope=0,lty=2)+facet_wrap(~facet_label,scale = "free_x")+coord_flip())
    dev.off()
    
    ####Figure S12h
    pdf("LOLA_signif_combined.pdf",height=3,width=4)
    print(ggplot(sub[sample=="allEloci__all1vsall2__1vs2all"&p.adjust<0.001],aes(x=dbSet,y=-log10(p.adjust),col=target,size=logOddsRatio))+geom_point()+geom_abline(intercept=-log10(0.001),slope=0,lty=2)+facet_wrap(~facet_label,scale = "free_x")+ylim(0,NA)+coord_flip()+ggtitle(paste0("N=", N_regions,"regions/",N_samples," samples")))
    print(ggplot(sub[sample=="allEloci__all1vsall2__1vs2all"&p.adjust<0.001],aes(x=dbSet,y=logOddsRatio,size=-log10(p.adjust),col=target))+geom_point()+facet_wrap(~facet_label,scale = "free_x")+coord_flip()+ggtitle(paste0("N=", N_regions,"regions/",N_samples," samples")))
    dev.off()
    setwd("../")
  }
  rm(list=grep("allEloci|rrbs_entropy_|rrbsAgEntropy_|methclone",objects(),value=TRUE))
}

