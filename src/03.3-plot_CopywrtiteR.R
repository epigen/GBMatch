library(project.init)
project.init2("GBMatch")
library(gsubfn)
#library(GenomicRanges)
library(gtools)
#library(rtracklayer)
#library(pheatmap)
library(ggrepel)


CNA_dir="CNAprofiles_single_100kb"

setwd(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/03-CopywriteR/results",CNA_dir))
controls=c("controls")

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

##get gene annotation
eload(loadGencodeGenes("human",versNum=87))
SV$gencodeContainer$genesGR$gene_name=SV$gencodeContainer$genes$external_gene_name
SV$gencodeContainer$genesGR$ensG=SV$gencodeContainer$genes$ensembl_gene_id
genes_gr=SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]]
genes_gr=genes_gr[!duplicated(genes_gr$gene_name)]

chroms_rle=Rle(as.factor(gsub("chr","",as.character(seqnames(genes_gr)))))
seqlevels(genes_gr)=levels(chroms_rle)
seqnames(genes_gr)=chroms_rle

genes_red=as.data.table(as.data.frame(genes_gr))
setnames(genes_red,names(genes_red),c("chrom","transc_start","transc_end","width","strand","gene","ensG"))


#only of single
control_samples=c("N1471_11A" ,"N1783_14" , "N321_15C"  ,"N61_15D"   ,"N912_12F")


#load segmentation data
load("segment.Rdata")
segments_dt=as.data.table(segment.CNA.object$output)
segments_dt[,sample:=unlist(lapply(ID,function(x){unlist(strsplit(gsub("log2.|.bam","",x),"\\.vs\\."))[1]})),]
segments_dt[,control:=unlist(lapply(ID,function(x){unlist(strsplit(gsub("log2.|.bam","",x),"\\.vs\\."))[2]})),]
segments_dt[,normalized:=ifelse(control=="none",FALSE,TRUE)]
segments_dt[,is_control:=ifelse(grepl(paste0(control_samples,collapse="|"),sample),TRUE,FALSE)]
segments_dt[,chrom:=gsubfn("23|24",list("23"="X","24"="Y"),as.character(chrom)),]
segments_dt=segments_dt[!sample%in%c(controls)]

#--------------------------

#annotate
annotation_sub=unique(annotation[,c("patID","N_number_seq","surgery.x","SurgeryDate","category","IDH"),with=FALSE])
setnames(annotation_sub,c("SurgeryDate","N_number_seq","surgery.x"),c("date","sample","surgery"))
segments_dt=merge(segments_dt,annotation_sub,by=c("sample"),all=TRUE)
segments_dt[,sample_short:=sample,]
#----------------------------------------------------

segments_ctr_dt=segments_dt[normalized==TRUE&num.mark>=8]
segments_ctr_dt[,Feature:=paste0("seg_",1:nrow(segments_ctr_dt)),]
segments_ctr_gr=with(segments_ctr_dt, GRanges(seqnames = Rle(chrom), IRanges(start=loc.start, end=loc.end),strand=Rle("*"),ID=Feature))

#load tiles
load("log2_counts_summary.Rdata")
counts_stats=log2_counts_summary
counts_stats[,Chromosome:=gsub("chr","",Chromosome),]
counts_stats=counts_stats[abs(set_mean)<10]  #need to do this because of weirdly high values (up to 10^9)
counts_stats[,dataSet:="GB_match"]
counts_stats_gr=with(counts_stats, GRanges(seqnames = Rle(Chromosome), IRanges(start=Start, end=End),strand=Rle("*"),ID=Feature))


#overlap segments with tiles
overlap=as.data.table(findOverlaps(counts_stats_gr,segments_ctr_gr,type="within"))
segments_ctr_annot=cbind(segments_ctr_dt[overlap$subjectHits],setnames(counts_stats[overlap$queryHits],c("Feature"),c("Feature.ct")))
segments_ctr_annot_mean=segments_ctr_annot[,list(set_median=mean(set_median,na.rm=TRUE),set_mean=mean(set_mean,na.rm=TRUE),set_sd=sd(set_sd,na.rm=TRUE),Ntiles=.N),by=c("sample","ID","chrom","loc.start","loc.end","normalized","is_control","surgery","date","sample_short","num.mark","seg.mean","patID","category","IDH")]
#-----------------------------------------

segments_ctr_annot_mean[,set_sig:=ifelse((abs(seg.mean)>set_sd),ifelse(seg.mean<0,"deletion","amplification"),"nc"),]
segments_ctr_annot_mean[,chrom:=factor(chrom,levels=mixedsort(unique(chrom))),]
segments_ctr_annot_mean[,sdChange:=abs(seg.mean/set_sd),]
segments_ctr_annot_mean[,sdChange_top:=ifelse(sdChange>5,5,sdChange),]

#now plot
dir.create("summary")

#distribution of change over sd
pdf(file.path("summary/sdChange_distribution.pdf"),height=4,width=4)
ggplot(segments_ctr_annot_mean[!chrom%in%c("X","Y")],aes(x=sdChange_top,col=is_control))+geom_density()+xlim(c(0,5))
dev.off()

#distribution of number of marks (num.mark)
pdf(file.path("summary/numMarks_distribution.pdf"),height=4,width=8)
ggplot(segments_dt,aes(x=num.mark,col=is_control))+geom_density()+facet_wrap(~normalized,labeller = "label_both")
dev.off()

#segment distribution: controls, not normalized
pdf(file.path("summary/segMean_distribution.pdf"),height=4,width=8)
ggplot(segments_dt,aes(x=seg.mean,col=is_control))+geom_density()+facet_wrap(~normalized,labeller = "label_both")
dev.off()


#CNV profiles for all samples
height=nrow(segments_ctr_annot_mean[,.N,by=c("sample_short","date")])+4

pdf(file.path("summary/CNAprofiles.pdf"),height=height,width=16)
ggplot(segments_ctr_annot_mean,aes(y=seg.mean,yend=seg.mean,col=set_sig))+geom_segment(aes(x = loc.start, xend = loc.end))+facet_grid(patID+surgery+sample_short~chrom,scales="free_x")+ylab("log2 fc")+xlab("")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+scale_color_manual(values=c("amplification"="red","deletion"="green","nc"="black"))+ylim(c(-2,2))
dev.off()


#CNV frequency chromosome plot simple 
sig_CNA=segments_ctr_annot_mean[set_sig!="nc"]
sig_CNA[,frag_len:=loc.end-loc.start,]
sig_CNA[,length_rank:=rank(-frag_len,ties.method="random"),by=c("chrom","set_sig")]
sig_CNA[,length_rank:=ifelse(set_sig=="deletion",-length_rank,length_rank),]
sig_CNA_gr=with(sig_CNA, GRanges(seqnames = Rle(chrom), IRanges(start=loc.start, end=loc.end),strand=Rle("*")))


height=max(sig_CNA[,.N,by="chrom"]$N)*0.035+5
pdf(file.path("summary/CNA_combiProfile.pdf"),height=height,width=20)
ggplot(sig_CNA,aes(y=length_rank,yend=length_rank,col=set_sig))+geom_segment(aes(x = loc.start, xend = loc.end))+facet_grid(~chrom,scales="free_x")+xlab("")+ylab("")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+scale_color_manual(values=c("amplification"="red","deletion"="green","nc"="black"))
dev.off()


#create matrix for chromosome wide amplification/deletion

mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
cytoband <- as.data.table(getTable( ucscTableQuery(mySession, track="cytoBand", table="cytoBand")))
cytoband[,chrom:=gsub("chr","",chrom),]
cytoband[,chromArm:=substr(name,1,1),]
cytoband_gr=with(cytoband,GRanges(seqnames = Rle(chrom), IRanges(start=chromStart, end=chromEnd),strand=Rle("*"),name=name))
ov=as.data.frame(findOverlaps(sig_CNA_gr,cytoband_gr,type="any",select="all"))
sig_CNA_cb=unique(cbind(sig_CNA[ov$queryHits],cytoband[ov$subjectHits,c("chromArm"),with=FALSE]))

amp_del_dt=sig_CNA_cb[,list(var_bases=sum(frag_len)),by=c("sample","date","chrom","chromArm","set_sig","surgery","category","IDH")]
amp_del_wide=dcast(amp_del_dt,sample+date+category+IDH+surgery~chrom+chromArm+set_sig,value.var="var_bases")
amp_del_wide[is.na(amp_del_wide)]=0

write.table(amp_del_wide,file.path("summary/CNA_Chromosome.tsv"),sep="\t",row.names=FALSE,quote=FALSE)


#also annoatete complete table with chromosome arm (1p 19q plot)
segments_ctr_annot_mean_gr=with(segments_ctr_annot_mean, GRanges(seqnames = Rle(chrom), IRanges(start=loc.start, end=loc.end),strand=Rle("*")))
ov_all=as.data.frame(findOverlaps(segments_ctr_annot_mean_gr,cytoband_gr,type="any",select="all"))
segments_ctr_annot_mean_cb=unique(cbind(segments_ctr_annot_mean[ov_all$queryHits],cytoband[ov_all$subjectHits,c("chromArm"),with=FALSE]))
annota_sub2=annotation[,c("N_number_seq","WHO2016_classification"),with=FALSE]
setnames(annota_sub2,"N_number_seq","sample_short")
segments_ctr_annot_mean_cb=merge(segments_ctr_annot_mean_cb,annota_sub2,by="sample_short")
segments_ctr_annot_mean_cb[,loc.len:=loc.end-loc.start,]
#mark bad samples (deletion in all segments) --> probably due to super low coverage
segments_ctr_annot_mean_cb[,status:=ifelse((sum(loc.len[set_sig=="deletion"])/2)>(sum(loc.len[set_sig!="deletion"])),"fail","pass"),by="sample_short"]

segments_ctr_annot_mean_cb_1p19q=segments_ctr_annot_mean_cb[(chrom==1&chromArm=="p")|(chrom==19&chromArm=="q")][status=="pass"]
segments_ctr_annot_mean_cb_1p19q[,region:=paste0(chrom,"_",chromArm),]
segments_ctr_annot_mean_cb_1p19q_red=segments_ctr_annot_mean_cb_1p19q[,list(loc.len=sum(loc.len)),by=c("sample_short","WHO2016_classification","category","IDH","region","set_sig")]

segments_ctr_annot_mean_cb_1p19q_red=segments_ctr_annot_mean_cb_1p19q_red[,list(dom_set=ifelse((max(c(0,loc.len[set_sig=="deletion"]))/sum(loc.len)>0.1),"deletion",set_sig[which.max(loc.len)]),percent=ifelse((max(c(0,loc.len[set_sig=="deletion"]))/sum(loc.len)>0.1),loc.len[set_sig=="deletion"]/sum(loc.len),max(loc.len)/sum(loc.len))),by=c("sample_short","WHO2016_classification","category","IDH","region")]

segments_ctr_annot_mean_cb_1p19q_wide=reshape(segments_ctr_annot_mean_cb_1p19q_red,idvar=c("sample_short","WHO2016_classification","category","IDH"),timevar="region",direction="wide")
segments_ctr_annot_mean_cb_1p19q_wide[,mean_len_perc:=mean(c(percent.1_p,percent.19_q)),by="sample_short"]

pdf(file.path("summary/1p19q.pdf"),height=3,width=6)
ggplot(segments_ctr_annot_mean_cb_1p19q_wide,aes(x=dom_set.1_p,y=dom_set.19_q,fill=WHO2016_classification,size=mean_len_perc))+geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0.3, dodge.width = NULL),alpha=0.5)+xlab(label="Chromosome 1p")+ylab(label="Chromosome 19q")+scale_fill_manual(values=c("red","grey","grey","grey","grey","red","grey"))+scale_size_continuous(range=c(1,4))
dev.off()


#overlap significant CNVs with genes
overlap=as.data.table(findOverlaps(genes_gr,sig_CNA_gr,type="within"))
sig_CNA_annot=cbind(sig_CNA[overlap$subjectHits],genes_red[overlap$queryHits])


#Select most frequent aberrations/genes and annotate all samples accordingly
genes_frequency=sig_CNA_annot[category=="GBMatch"&IDH=="wt",list(N=.N,min_frag_len_sep=min(frag_len),med_frag_len_sep=median(frag_len)),by=c("gene","set_sig","chrom","transc_start","surgery")]
genes_frequency[,rank_freq:=rank(-N,ties.method="min"),by=c("set_sig","surgery")]
genes_frequency[,rank_len:=rank(med_frag_len_sep,ties.method="min"),by=c("set_sig","surgery")]
genes_frequency[,max_rank:=pmax(rank_freq,rank_len),]
genes_frequency[,combined_rank:=rank(max_rank,ties.method="random"),by=c("set_sig","surgery")]

#make wide 
genes_frequency[,min_frag_len:=min_frag_len_sep[which.max(N)],by=c("gene")]
genes_frequency_wide=as.data.table(dcast(genes_frequency,gene+chrom+min_frag_len+surgery~set_sig,value.var="N"))
genes_frequency_wide[is.na(amplification),amplification:=0,]
genes_frequency_wide[is.na(deletion),deletion:=0,]
genes_frequency_wide[,ad_ratio:=log2(amplification/deletion),]
genes_frequency_wide[,type:=ifelse(ad_ratio>0,"amplification",ifelse(ad_ratio<0,"deletion","tie")),]
genes_frequency_wide[,frequency:=pmax(amplification,deletion),]
genes_frequency_wide[,rank_freq:=rank(-frequency,ties.method="min"),by="type"]
genes_frequency_wide[,rank_rat:=rank(-abs(ad_ratio),ties.method="min"),by="type"]
genes_frequency_wide[,rank_len:=rank(min_frag_len,ties.method="min"),by="type"]
genes_frequency_wide[,max_rank:=pmax(rank_freq,rank_rat,rank_len),]
genes_frequency_wide[,combined_rank:=rank(max_rank,ties.method="random"),by="type"]

write.table(genes_frequency_wide[order(combined_rank)],"summary/CNV_genes.tsv",sep="\t",quote=FALSE,row.names=FALSE)


cat(genes_frequency_wide[order(combined_rank)][type=="amplification"]$gene[1:500],sep="\n",file="summary/top_amplified_genes.txt")
cat(genes_frequency_wide[order(combined_rank)][type=="deletion"]$gene[1:500],"\n",file="summary/top_deleted_genes.txt")


#selcet genes of interest
GOIs_sel=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/Known_GB_relevant_Genes_Nik.tsv"))
setnames(GOIs_sel,"x","gene")
OGTS=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/TS-OG_1235122TablesS1-4_vogelsteinetal_2013_red.csv"))

GOIs_sel=merge(GOIs_sel,genes_red,by="gene")
GOIs_sel=merge(GOIs_sel,OGTS,by="gene",all.x=TRUE)
GOIs_sel[,sel:=TRUE,]
GOIs_sel=GOIs_sel[,c("gene","chrom","transc_start","sel","classification"),with=FALSE]
GOIs_sel_ext=rbindlist(list(copy(GOIs_sel[,surgery:=1,]),copy(GOIs_sel[,surgery:=2,]),copy(GOIs_sel[,surgery:=3,]),copy(GOIs_sel[,surgery:=4,]),copy(GOIs_sel[,surgery:=5,]),copy(GOIs_sel[,surgery:=6,])))


genes_frequency[,keep:=set_sig[which.min(min_frag_len_sep)],by=c("gene","surgery")]
GOIs_diff=genes_frequency[set_sig!="tie"&!is.na(surgery)&keep==set_sig,c("gene","set_sig","chrom","transc_start","surgery"),with=FALSE]
GOIs_comb=merge(GOIs_diff,GOIs_sel_ext[surgery%in%c(1,2)],by=c("gene","chrom","transc_start","surgery"),all.y=TRUE)
GOIs_comb[is.na(classification),classification:="n.d.",]

GOIs_comb=GOIs_comb[chrom%in%c(as.character(1:22),"X","Y")&sel==TRUE]
GOIs_comb[is.na(set_sig),set_sig:="sel"]
GOIs_comb[,ypos:=ifelse(set_sig=="deletion",-200,ifelse(set_sig=="amplification",200,0)),]

sig_CNA[,length_rank_sep:=rank(abs(length_rank),ties.method="random"),by=c("surgery","set_sig","chrom")]
sig_CNA[,length_rank_sep:=ifelse(set_sig=="deletion",-length_rank_sep,length_rank_sep),]

GOIs_comb[,chrom:=factor(chrom,levels=c(as.character(1:22),"X","Y")),]
sig_CNA[,chrom:=factor(chrom,levels=c(as.character(1:22),"X","Y")),]

pdf(file.path("summary/CNA_combiProfile_genes.pdf"),height=height*0.5,width=30)
ggplot(sig_CNA[surgery%in%c(1,2)&category=="GBMatch"&IDH=="wt"])+geom_segment(aes(col=set_sig,y=length_rank_sep,yend=length_rank_sep,x = loc.start, xend = loc.end))+geom_point(data=GOIs_comb,aes(x=transc_start,y=ypos))+geom_text_repel(data=GOIs_comb,force=20,aes(x=transc_start,y=ypos,label=gene,col=classification))+facet_wrap(~surgery+chrom,scales="free_x",ncol=24)+xlab("")+ylab("")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+scale_color_manual(values=c("ND"="black","DRG"="blue","OG"="red","TSG"="green","amplification"="red","deletion"="green","nc"="black"))+ylim(c(-250,250))
dev.off()


