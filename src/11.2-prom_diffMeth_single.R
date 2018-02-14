library(project.init)
project.init2("GBMatch")
library(LOLA)
library(RGenomeUtils)
library(enrichR)



##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/11.2-diffMeth_single/")
dir.create(out_dir)
setwd(out_dir)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

annot_ASC=annotation[IDH=="wt",list(sample.1=N_number_seq[surgery.x==1][1],cycles.1=enrichmentCycles[surgery.x==1][1],cycles.2=enrichmentCycles[surgery.x==2][1],Age=min(Age),Sex=unique(Sex), extentOfResection=Extentofresection[surgery.x==1][1],progression_location=progression_location[surgery.x==2][1],timeToSecSurg=unique(timeToSecSurg)),by=patID]



#chose here
##list all available diff meths
diffMeth_files=list.files(file.path(getOption("RCACHE.DIR"),paste0("diffMeth/Prom1k")),full.names=TRUE)

#patient wise combi for GBMatch
combinations=annotation[category=="GBMatch"&IDH=="wt",list(N1=N_number_seq[surgery.x==1],N2=N_number_seq[surgery.x==2]),by=c("patientID","patID")]
combinations[,cache:=unlist(lapply(paste0(N1,"vs",N2,".RData"),function(x){ret=grep(x,diffMeth_files,value=TRUE);return(paste0(ret,""))})),]
stopifnot(nrow(combinations[cache==""])==0)


#set significance criteria
diffmeth_plim=0.001
diffmeth_diff=0.75
diffmeth_readlim=20
all_signif=data.table()
all=data.table()

for (i in 1:nrow(combinations)){

  print(combinations[i]$cache)
  rm(ret)
  load(combinations[i]$cache)
  ret[,signif:=p.value.adj<diffmeth_plim&abs(methyl.1-methyl.2)>diffmeth_diff&readCount.1>diffmeth_readlim&readCount.2>diffmeth_readlim]
  
  all_signif=rbindlist(list(all_signif,ret[signif==TRUE]))
  all=rbindlist(list(all,ret))

}

#add patID
all_signif[,hypo_meth:=ifelse(methyl.1<methyl.2,1,ifelse(methyl.1>methyl.2,2,"tie")),]
all_signif[,sample.1:=sub("RRBS_cpgMethylation_","",sampleName.1),]
all_signif[,patID:=sub("_[1-9]$","",id,perl=TRUE),]

all[,sample.1:=sub("RRBS_cpgMethylation_","",sampleName.1),]
all[,patID:=sub("_[1-9]$","",id,perl=TRUE),]


sub=all_signif[,.N,by=c("gene_name","hypo_meth")][gene_name%in%gene_name[N>4]]
sub[,N:=ifelse(hypo_meth==1,N,-N),]
sub[,plot_order:=-N[which.max(abs(N))],by="gene_name"]
#fix order for DCN (equel + and -)
sub[gene_name=="DCN",plot_order:=-4,]
sub[,gene_name:=factor(gene_name,levels=unique(gene_name[order(plot_order,decreasing=TRUE)])),]


sub_expand <- rbind(sub, cbind(expand.grid(gene_name=levels(sub$gene_name), hypo_meth=levels(as.factor(sub$hypo_meth))), N=NA,plot_order=NA))

pdf("promoter_diff_meth_recurrent.pdf",width=5,height=5)
ggplot(sub_expand,aes(x=gene_name,y=N,fill=factor(hypo_meth,levels=c("2","1")),group=hypo_meth,alpha=ifelse(abs(N)<5,"off","on")))+geom_bar(col="black",stat="identity",position=position_dodge(width=0))+theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))+scale_alpha_discrete(range=c(0.25,1))+theme(legend.position="bottom",legend.box = "vertical")
dev.off()

sub_signif=all_signif[gene_name%in%sub$gene_name]
sub_all=all[readCount.1>diffmeth_readlim&readCount.2>diffmeth_readlim]
mycor=sub_all[,cor(methyl.1,methyl.2)]

pdf("promoter_diff_meth_scatter.pdf",width=6.5,height=4)
ggplot(data=sub_all,aes(x=methyl.1*100,y=methyl.2*100))+stat_bin_hex(bins=30,aes(fill=..density..^0.1))+geom_point(data=sub_signif,shape=21,size=3,aes(col=factor(hypo_meth,levels=c("2","1"))))+annotate("text",x=15,y=98,label=paste0("r=",round(mycor,3)))+scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+xlab("% DNA methylation primary tumor")+ylab("% DNA methylation recurring tumor")+theme(aspect.ratio=1)
dev.off()

#follow-up on recurrently differentially methylated genes
sub_rec=sub[abs(N)>=5]
sub_rec=sub_rec[,list(hypo_meth=as.character(sum(as.numeric(hypo_meth)))),by="gene_name"]
sub_rec[hypo_meth=="3",hypo_meth:="both",]

sel_recurring=merge(sub_all,sub_rec,by="gene_name")

sel_recurring[,gene_name:=factor(gene_name,levels=unique(gene_name[order(hypo_meth)])),]
pdf("rec_promoter_all_boxpl.pdf",width=7,height=4)
ggplot(sel_recurring,aes(x=gene_name,y=methyl.2-methyl.1,fill=hypo_meth,col=hypo_meth))+geom_point(shape=21,alpha=0.5,position=position_jitter(height=0,width=0.2))+geom_boxplot(col="black",fill="transparent",outlier.shape=NA)+geom_hline(yintercept=c(0.8,-0.8),lty=20,col="grey")+theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))
dev.off()

#now group patients according to methylation trend of recurrently differentially methylated genes
sel_recurring_trend=sel_recurring[hypo_meth!="both"]
sel_recurring_trend[,trend:=ifelse(hypo_meth=="1",1,-1),]
sel_recurring_trend[,diff_meth:=methyl.2-methyl.1,]
sel_recurring_trend_pat=sel_recurring_trend[,list(diff_trend_cor=cor(trend,diff_meth),diff_trend_dist=as.numeric(dist(rbind(trend,diff_meth),method="manhattan")),Ngenes=.N),by=c("patID")]
sel_recurring_trend_pat[,diff_trend_dist_norm:=diff_trend_dist/Ngenes,]

pdf("sel_recurring_trend.pdf",height=2.5,width=3.5)
ggplot(sel_recurring_trend_pat,aes(x=Ngenes,y=diff_trend_dist_norm))+geom_point(shape=21)
ggplot(sel_recurring_trend_pat,aes(x=diff_trend_dist_norm))+geom_histogram(col="darkgrey",fill="grey",binwidth=0.05)+ylim(c(0,20))+scale_x_continuous(breaks=seq(from=0,to=2,by=0.2),limits=c(0,2))+ggtitle("all patients")+geom_vline(xintercept=c(0.7,1.15),lty=20)
ggplot(sel_recurring_trend_pat[Ngenes>4],aes(x=diff_trend_dist_norm))+geom_histogram(col="darkgrey",fill="grey",binwidth=0.05)+ylim(c(0,20))+scale_x_continuous(breaks=seq(from=0,to=2,by=0.2),limits=c(0,2))+ggtitle(">4 genes")+geom_vline(xintercept=c(0.7,1.15),lty=20)
dev.off()

sel_recurring_trend=merge(sel_recurring,sel_recurring_trend_pat,by="patID")
sel_recurring_trend[,annot:=ifelse(diff_trend_dist_norm>=1,"anti-trend",ifelse(diff_trend_dist_norm<0.8,"trend",NA)),]
sel_recurring_trend[,gene_name:=factor(gene_name,levels=levels(sub_expand$gene_name)),]

pdf("promoter_diff_meth_recurrent_trend.pdf",width=5.15,height=3)
ggplot(sel_recurring_trend[!is.na(annot)],aes(x=gene_name,y=methyl.2-methyl.1))+geom_point(position=position_jitter(height=0,width=0.3),aes(group=patID,col=annot),alpha=0.2)+geom_smooth(se=FALSE,aes(group=annot,col=annot),span=0.3)+geom_hline(yintercept=c(0,0.75,-0.75),lty=20,col="black")+theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))+theme(legend.position="bottom",legend.box = "vertical")+scale_color_manual(values=c("navy","firebrick"))
dev.off()

write.table(sel_recurring_trend_pat,"sel_recurring_trend_pat.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#enrichr analysis of diffmeth genes
##run enrichr
sub_uniq=sub[abs(N)>=5]
enrichr_res_hypo_1=try(as.data.table(enrichGeneList(unique(sub_uniq[hypo_meth==1]$gene_name),databases = c("Panther_2016"),fdr.cutoff=0.05)),silent = TRUE)
enrichr_res_hypo_1[,hypo_meth:="1",]
enrichr_res_hypo_2=try(as.data.table(enrichGeneList(unique(sub_uniq[hypo_meth==2]$gene_name),databases = c("Panther_2016"),fdr.cutoff=0.05)),silent = TRUE)
enrichr_res_hypo_2[,hypo_meth:="2",]


enrichr_res_sub=rbindlist(list(enrichr_res_hypo_1[pval<0.05][order(genes)],enrichr_res_hypo_2[pval<0.05][order(genes)]))
enrichr_res_sub[,category_simpl:=gsub("\\(.*\\)|P0....|Homo sapiens","",gsub("_"," ",category)),]
enrichr_res_sub[,annot:=ifelse(grepl("devel|gastru|guidance",category),"Development",ifelse(grepl("DNA|apopto|Apopto",category),"Apoptosis",ifelse(grepl("Wnt|Cadherin",category),"Wnt signalling",ifelse(grepl("Pentose|Fructose|Glyco",category),"Glycolysis","Immune response")))),]

enrichr_res_sub[,category_simpl:=factor(category_simpl,levels=category_simpl[order(-pval)]),]

pdf("promoter_diff_meth_enrichr.pdf",width=6,height=3)
ggplot(enrichr_res_sub,aes(x=category_simpl,y=-log10(pval)))+geom_bar(stat="identity",aes(fill=annot),alpha=0.6)+coord_flip()+geom_text(aes(label=genes),hjust=1)+scale_fill_manual(values=c("Development"="#7fc97f","Apoptosis"="#beaed4","Immune response"="#fb9a99","Wnt signalling"="#fdc086"))+geom_hline(yintercept=-log10(0.05),lty=20,col="grey")+facet_wrap(~hypo_meth,labeller="label_both",scale="free_y",ncol=1)+xlab("")
dev.off()

#now group samples according to methylation status in developmental/ DNA damage response genes (from enrichr)
enrichr_genes=fread(file.path(getOption("PROJECT.DIR"),"metadata/gene_lists/Panther_2016_from_enrichr"),col.names=c("category","genes"))

term_list=c("Development","Apoptosis","Wnt signalling","Immune response")
genes_enrichr=data.table()
for (term in term_list){
  genes=unique(unlist(strsplit(enrichr_genes[category%in%enrichr_res_sub[annot==term]$category]$genes,"\t")))
  genes_dt=data.table(term=term,gene_name=genes)
  genes_enrichr=rbindlist(list(genes_enrichr,genes_dt))
}

genes_enrichr[,N_terms:=.N,by="gene_name"]

sub_all_term=merge(sub_all,genes_enrichr[N_terms<2],by="gene_name")
sub_all_term[,diff_meth:=methyl.2-methyl.1,]
sub_all_term_pat=sub_all_term[,list(mean_diffmeth=mean(diff_meth,na.rm=TRUE),sd_diffmeth=sd(diff_meth,na.rm=TRUE),median_diffmeth=median(diff_meth,na.rm=TRUE),Ngenes=.N),by=c("patID","term")]

sub_all_term_pat_wide=reshape(sub_all_term_pat[Ngenes>4],idvar="patID",timevar="term",drop=c("Ngenes","median_diffmeth","sd_diffmeth"),direction="wide")
write.table(sub_all_term_pat_wide,"promoter_diff_meth_enrichr_term.tsv",sep="\t",quote=FALSE,row.names=FALSE)


#relative diffmeth comparisons
min_diff=0
sub_all_agg=sub_all[,list(signif=sum(signif),all=.N),by=c("sampleName.1","sampleName.2")]
sub_all_agg[,sample.1:=sub("RRBS_cpgMethylation_","",sampleName.1),]
sub_all_agg[,DPM:=signif/all*1000000,]
sub_all_agg=merge(annot_ASC,sub_all_agg,by="sample.1")

write.table(sub_all_agg,"DPM.tsv",sep="\t",quote=FALSE,row.names=FALSE)

sub_cyc=sub_all_agg[cycles.1<16&cycles.2<16&cycles.1>12&cycles.2>12]
sub_cyc[,log10_DPM:=log10(DPM),]
cors=sub_cyc[,cor.test(timeToSecSurg,log10_DPM,alternative="less"),]

pdf("DPT_comps.pdf",height=3.5,width=3)
ggplot(sub_cyc,aes(y=log10_DPM,x=timeToSecSurg/30))+geom_point(shape=21,size=3,fill="grey")+geom_smooth(method="lm",fill="lightgrey")+annotate("text",x=-Inf,y=Inf,label=paste0("r=",round(cors$estimate,2),"\np.value=",round(cors$p.value,2)),hjust=-0.2,vjust=1.2)+xlab("timeToSecSurg (months)")

ggplot(sub_cyc,aes(y=log10_DPM,x=as.factor(progression_location),group=progression_location))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(fill="transparent",outlier.shape=NA)

ggplot(sub_cyc,aes(y=log10_DPM,x=as.factor(extentOfResection),group=extentOfResection))+geom_point(shape=21,position=position_jitter(width=0.2))+geom_boxplot(outlier.shape=NA,fill="transparent")

dev.off()

#follow-up on Wnt genes (specifically SFRP2)

#get genes
eload(loadGencodeGenes("human",versNum=87))
prom1k=promoters(SV$gencodeContainer$genesGR[SV$gencodeContainer$genes[,which(gene_biotype=="protein_coding")]], upstream=1000, downstream=500)
genes=SV$gencodeContainer$genes[gene_biotype=="protein_coding"]
genes_annot=cbind(as.data.table(as.data.frame(prom1k)),genes[,c("external_gene_name","ensembl_gene_id","chromosome_name"),])
stopifnot(genes_annot$seqnames==genes_annot$chromosome_name)
setnames(genes_annot,"seqnames","chr")
#for some reasons 20 genes have different ENSGs but same position and gene name --> remove
genes_annot=genes_annot[!duplicated(cbind(chr,start,end)),]


#get measures of promoter DNA Meth and heterogeneity

simpleCache("pdrPromoters1ksub", assignToVariable="pdr_ag")
simpleCache("entropy/rrbsAgEntropy_min40_prom1k",assignToVariable="entropy_ag")
simpleCache("rrbsProm1kb",assignToVariable="meth_ag")
setnames(meth_ag,"CpGcount","CpGcount_meth")

merged_heterogeneity_pre=merge(pdr_ag,entropy_ag,by=c("id","chr","start","end","regionID"),all=TRUE)
merged_heterogeneity=merge(merged_heterogeneity_pre,meth_ag,by=c("id","chr","start","end","regionID"),all=TRUE)
setnames(merged_heterogeneity,"id","N_number_seq")

merged_heterogeneity_annot=merge(merged_heterogeneity,annotation[,c("N_number_seq","rand_frag_perc","Unique_CpGs","enrichmentCycles","Ct_postLigation","Ct_postConversion","Follow-up_years","Age","surgery.x","patID","category","IDH"),with=FALSE],by="N_number_seq")
merged_heterogeneity_annot_genes=merge(genes_annot,merged_heterogeneity_annot,by=c("chr","start","end"))


#store meth values for other analysis
SFRP2_meth=unique(merged_heterogeneity_annot_genes[external_gene_name%in%c("SFRP2"),c("N_number_seq","methyl","readCount","CpGcount_meth"),with=FALSE])
setnames(SFRP2_meth,c("methyl","readCount","CpGcount_meth"),c("SFRP2_meth","SFRP2_readCount","SFRP2_CpG_count"))
write.table(SFRP2_meth,"SFRP2_meth.tsv",quote=FALSE,sep="\t",row.names=FALSE)


#some plots

sub=merged_heterogeneity_annot_genes[external_gene_name%in%c("SFRP2")&category%in%c("GBMatch","GBmatch_val")&surgery.x%in%c(1,2)&IDH=="wt"&(readCount/CpGcount_meth)>10&CpGcount_meth>20]
sub[,NsurgAvail:=.N,by=patID]
sub_long=melt(sub,id.vars=c("patID","N_number_seq","category","external_gene_name","surgery.x","NsurgAvail"),measure.vars=c("methyl","PDRa"))
#remark: too little coverage in entropy...

pdf("SFRP2_followup_boxpl.pdf",height=3,width=5)
ggplot(sub_long,aes(y=value,x=paste0(category,"\nSurgery ",surgery.x)))+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=0.3),alpha=0.5)+facet_wrap(variable~external_gene_name,ncol=3)+xlab("")
dev.off()

pdf("SFRP2_followup_cor.pdf",height=4,width=5.5)
sub[,meth_pdr_cor:=cor(methyl,PDRa),by=c("surgery.x","category","external_gene_name")]
ggplot(sub,aes(y=PDRa,x=methyl,col=paste0(category," ",external_gene_name,"\nSurgery ",surgery.x," r=",round(meth_pdr_cor,3))))+geom_point(alpha=0.6)+geom_smooth(method="lm",alpha=0.2)+scale_color_discrete(name="Chort")+facet_wrap(~external_gene_name)
dev.off()

pdf("SFRP2_followup_trend.pdf",height=6,width=5)
ggplot(sub_long[NsurgAvail>1],aes(y=value, x=paste0("Surgery ",surgery.x),group=paste0(patID,"_",variable)))+geom_line(aes(lty=variable))+geom_point(aes(shape=variable),position=position_jitter(width=0.01),alpha=0.6)+facet_wrap(external_gene_name~patID,ncol=2)+xlab("")
dev.off()

