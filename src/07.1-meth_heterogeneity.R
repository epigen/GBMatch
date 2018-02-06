library(project.init)
project.init2("GBMatch")
library(LOLA)
library(RGenomeUtils)
library(hexbin)
library("RColorBrewer")
library("gridExtra")

##set directories
out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/07-meth_heterogeneity/")
dir.create(out_dir)
setwd(out_dir)

##get annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
setnames(annotation,"N_number_seq","id")


#get measures of DNA Meth heterogeneity

simpleCache("pdrPromoters1ksub", assignToVariable="pdr_ag")
simpleCache("entropy/rrbsAgEntropy_min40_prom1k",assignToVariable="entropy_ag")
simpleCache("rrbsProm1kb",assignToVariable="meth_ag")
setnames(meth_ag,"CpGcount","CpGcount_meth")

merged_heterogeneity_pre=merge(pdr_ag,entropy_ag,by=c("id","chr","start","end","regionID"),all=TRUE)
merged_heterogeneity=merge(merged_heterogeneity_pre,meth_ag,by=c("id","chr","start","end","regionID"),all=TRUE)

merged_heterogeneity_annot=merge(merged_heterogeneity,annotation[,c("id","rand_frag_perc","Unique_CpGs","enrichmentCycles","Ct_postLigation","Ct_postConversion","Follow-up_years","Age","surgery.x","patID","category","IDH"),with=FALSE],by="id")


#define regions to include (core regions)
sample_count=length(unique(merged_heterogeneity$id))
core_regions=merged_heterogeneity[,list(pdr=ifelse(sum(!is.na(PDRa))>0.75*sample_count,TRUE,FALSE),entropy=ifelse(sum(!is.na(entropy1))>0.10*sample_count,TRUE,FALSE),meth=ifelse(sum(!is.na(methyl))>0.75*sample_count,TRUE,FALSE)),by=c( "chr","start", "end", "regionID")]
regions_count=core_regions[,list(pdr=sum(pdr),entropy=sum(entropy),meth=sum(meth)),]
regions_count

all_means=merged_heterogeneity[!chr%in%c("chrY","chrX"),list(mean_entropy=mean(entropy1[regionID%in%core_regions[entropy==TRUE]$regionID],na.rm=TRUE),mean_pdr=mean(PDRa[(readCount/CpGcount_meth)>20&regionID%in%core_regions[pdr==TRUE]$regionID],na.rm=TRUE),mean_meth=mean(methyl[(readCount/CpGcount_meth)>20&regionID%in%core_regions[meth==TRUE]$regionID],na.rm=TRUE)),by="id"]

#chach this to include in general annotation and use in other analyses (e.g. survival)
simpleCache("combined_heterogeneity",all_means,recreate=FALSE,assignToVariable="all_means")


all_means_annot=merge(all_means,annotation[,c("id","rand_frag_perc","Unique_CpGs","enrichmentCycles","Ct_postLigation","Ct_postConversion","Follow-up_years","Age","surgery.x","patID","category","IDH"),with=FALSE],by="id")


#plot heterogeneity per surgery
sub_cycles=all_means_annot[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2)&enrichmentCycles<16&enrichmentCycles>12]
sub_cycles[,annot_entropy:=ifelse(mean_entropy>quantile(mean_entropy,0.8),"high",ifelse(mean_entropy<quantile(mean_entropy,0.2),"low",NA)),by=c("surgery.x","category")]
sub_cycles[,annot_pdr:=ifelse(mean_pdr>quantile(mean_pdr,0.8),"high",ifelse(mean_pdr<quantile(mean_pdr,0.2),"low",NA)),by=c("surgery.x","category")]

#calculate significance of difference


sign_1_pdr=sub_cycles[,wilcox.test(x=mean_pdr[surgery.x==1&category=="GBMatch"],y=mean_pdr[surgery.x==2&category=="GBMatch"])$p.value,]
sign_2_pdr=sub_cycles[,wilcox.test(x=mean_pdr[surgery.x==1&category=="GBMatch"],y=mean_pdr[surgery.x==1&category=="GBmatch_val"])$p.value,]
sign_3_pdr=sub_cycles[,wilcox.test(x=mean_pdr[surgery.x==2&category=="GBMatch"],y=mean_pdr[surgery.x==1&category=="GBmatch_val"])$p.value,]

sign_1_epy=sub_cycles[,wilcox.test(x=mean_entropy[surgery.x==1&category=="GBMatch"],y=mean_entropy[surgery.x==2&category=="GBMatch"])$p.value,]
sign_2_epy=sub_cycles[,wilcox.test(x=mean_entropy[surgery.x==1&category=="GBMatch"],y=mean_entropy[surgery.x==1&category=="GBmatch_val"])$p.value,]
sign_3_epy=sub_cycles[,wilcox.test(x=mean_entropy[surgery.x==2&category=="GBMatch"],y=mean_entropy[surgery.x==1&category=="GBmatch_val"])$p.value,]


pdf("means_comp.pdf",height=3.5,width=3)
ggplot(sub_cycles,aes(y=mean_entropy,x=paste0(category,"\nSurg: ",surgery.x),group=paste0(category,"\nSurg: ",surgery.x)))+geom_point(shape=21,size=3,col="grey",position=position_jitter(width=0.2,height=0))+geom_boxplot(fill="transparent",outlier.shape=NA)+guides(fill=FALSE)

ggplot(sub_cycles,aes(y=mean_entropy,x=paste0(category,"\nSurg: ",surgery.x),group=paste0(category,"\nSurg: ",surgery.x)))+geom_point(aes(fill=annot_entropy),shape=21,size=3,col="grey",alpha=0.7,position=position_jitter(width=0.2,height=0))+geom_boxplot(fill="transparent",outlier.shape=NA)+annotate("text",x=0.5,y=30,label=paste0(round(c(sign_1_epy,sign_2_epy,sign_3_epy),6),collapse="  "),hjust=0)+guides(fill=FALSE)+scale_fill_brewer(palette="Set2")

ggplot(sub_cycles,aes(y=mean_pdr,x=paste0(category,"\nSurg: ",surgery.x),group=paste0(category,"\nSurg: ",surgery.x)))+geom_point(shape=21,size=3,col="grey",position=position_jitter(width=0.2,height=0))+geom_boxplot(fill="transparent",outlier.shape=NA)+guides(fill=FALSE)

ggplot(sub_cycles,aes(y=mean_pdr,x=paste0(category,"\nSurg: ",surgery.x),group=paste0(category,"\nSurg: ",surgery.x)))+geom_point(aes(fill=annot_pdr),shape=21,size=3,col="grey",alpha=0.7,position=position_jitter(width=0.2,height=0))+geom_boxplot(fill="transparent",outlier.shape=NA)+guides(fill=FALSE)+annotate("text",x=0.5,y=0.27,label=paste0(round(c(sign_1_pdr,sign_2_pdr,sign_3_pdr),6),collapse="  "),hjust=0)+guides(fill=FALSE)+scale_fill_brewer(palette="Set2")
dev.off()


#compare single regions
merged_heterogeneity_annot_red=na.omit(merged_heterogeneity_annot)[enrichmentCycles<16&enrichmentCycles>12&category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&CpGcount>15]

#scatter plots
entropy_pdr_cor=merged_heterogeneity_annot_red[,cor(entropy1,PDRa),by=c("surgery.x","category")]

pdf("compare_heterogeneity.pdf",height=3,width=9)
ggplot(merged_heterogeneity_annot_red,aes(x=PDRa,y=entropy1))+stat_binhex(aes(fill=log10(..count..)))+scale_fill_gradient(low="white",high="blue")+geom_smooth(se=FALSE,col="black")+geom_text(data=entropy_pdr_cor,x=0.1,y=190,aes(label=paste0("r=",round(V1,3))))+facet_wrap(category~surgery.x)
ggplot(merged_heterogeneity_annot_red,aes(x=methyl,y=entropy1))+stat_binhex(aes(fill=log10(..count..)))+scale_fill_gradient(low="white",high="blue")+geom_smooth(se=FALSE,col="black")+facet_wrap(category~surgery.x)
ggplot(merged_heterogeneity_annot_red,aes(x=methyl,y=PDRa))+stat_binhex(aes(fill=log10(..count..)))+scale_fill_gradient(low="white",high="blue")+geom_smooth(se=FALSE,col="black")+facet_wrap(category~surgery.x)
dev.off()


#assess epi-allel distribution in high and low heterogeneity samples
merged_heterogeneity_annot_red=merge(merged_heterogeneity_annot_red,sub_cycles[,c("id","annot_entropy","annot_pdr"),with=FALSE],by="id")
merged_heterogeneity_annot_red_long=melt(merged_heterogeneity_annot_red,id.vars=c("id","patID","category","Follow-up_years","surgery.x","annot_entropy", "annot_pdr","regionID","entropy1","PDRa"),measure.vars=grep("^p[10].*",names(merged_heterogeneity_annot_red),value=TRUE,perl=TRUE))

categories=c("GBMatch","GBmatch_val")
for (sel_category in categories){

condensed=merged_heterogeneity_annot_red_long[category==sel_category,list(sum=sum(value),Nregions=.N),by=c("id","patID","category","Follow-up_years","variable","annot_entropy", "annot_pdr","surgery.x")]
condensed[,rel_freq:=sum/(Nregions),]
condensed[,plotID:=paste0(patID,"_",surgery.x),]
condensed[,variable:=factor(variable,levels=c("p0000","p0001","p0010","p0100","p1000","p0011","p0101","p1001","p0110","p1010","p1100","p0111","p1011","p1101","p1110","p1111"))]
condensed[,plotID:=factor(plotID,levels=unique(plotID[order(rel_freq[variable=="p0000"])]))]

pl1=ggplot(condensed[!(is.na(annot_entropy)&is.na(annot_pdr))&variable!="p0000"],aes(x=plotID,y=rel_freq,fill=variable))+geom_point(y=-2,aes(col=annot_entropy))+geom_point(y=-5,aes(col=annot_pdr))+geom_bar(stat="identity")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5), legend.box = "vertical",legend.position="bottom")+scale_fill_manual(values=colorRampPalette(brewer.pal(name="Blues", n = 9)[2:9])(15))+ylim(c(-6,22))+facet_wrap(~surgery.x,scale="free")+scale_color_brewer(palette="Set2")+ guides(fill = guide_legend(nrow = 2))

pl2=ggplot(unique(condensed[,c("plotID","category","surgery.x","annot_entropy","annot_pdr"),with=FALSE])[!(is.na(annot_entropy)&is.na(annot_pdr))],aes(x=plotID))+geom_density(aes(group=annot_pdr,col=annot_pdr,fill=annot_pdr),alpha=0.2)+facet_wrap(~surgery.x,scale="free")+annotate("text",label="PDR class",x=1,hjust=0,y=0)+guides(col=FALSE,fill=FALSE)+theme(axis.title.x=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")

pl3=ggplot(unique(condensed[,c("plotID","category","surgery.x","annot_entropy","annot_pdr"),with=FALSE])[!(is.na(annot_entropy)&is.na(annot_pdr))],aes(x=plotID))+geom_density(aes(group=annot_entropy,col=annot_entropy,fill=annot_entropy),alpha=0.2)+facet_wrap(~surgery.x,scale="free")+annotate("text",label="Entropy class",x=1,hjust=0,y=0)+guides(col=FALSE,fill=FALSE)+theme(axis.title.x=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")

pdf(paste0("epiallele-freqs_",sel_category,".pdf"),height=6,width=ifelse(sel_category=="GBMatch",8,7.5))
grid.arrange(pl3, pl2, pl1, ncol=1, nrow=3, heights=c(1,1,5))
dev.off()
}

##investigate effect of enrichment cycles
#make correlation plot of mean heterogeneity measures and multiple quality measures
my.panel.cor=function(x, y, digits=3, prefix="", cex.cor=2, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  cex.cor=7;
  text(0.5, 0.5, txt, cex = max(c(cex.cor * abs(r/2),1)))
}

pdf("means_cor_explain.pdf",height=20,width=20)
pairs(all_means_annot[category=="GBMatch"&IDH=="wt",-c("id","surgery.x","patID","category","IDH"),with=FALSE],upper.panel= my.panel.cor,main="GBMatch all selected regions")
pairs(all_means_annot[category=="GBmatch_val"&IDH=="wt",-c("id","surgery.x","patID","category","IDH"),with=FALSE],upper.panel= my.panel.cor,main="GBmatch_val all selected regions")
dev.off()


#compare surgery 1 and 2

sub=all_means_annot[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&surgery.x%in%c(1,2),c("mean_entropy",  "mean_pdr", "mean_meth","surgery.x","patID","enrichmentCycles","category"),with=FALSE]
all_means_annot_long=melt(sub,id.vars=c("patID","surgery.x","enrichmentCycles","category"))
all_means_annot_wide=reshape(all_means_annot_long, idvar=c("patID","variable","category"),timevar="surgery.x",direction="wide")

coords=all_means_annot_wide[,list(max=max(c(unlist(value.1),unlist(value.2)),na.rm=TRUE)),by=c("variable","category")]

pdf(paste0("heterogeneity_cycles.pdf"),height=3,width=11)
ggplot()+geom_segment(aes(x = 0, y = 0, xend = max, yend = max),lty=20,lwd=2,colour = "grey",data=coords)+geom_point(data=all_means_annot_wide[category=="GBMatch"],aes(x=value.1,y=value.2,fill=log2(enrichmentCycles.1/enrichmentCycles.2)),shape=21,size=3)+facet_wrap(~variable, scales="free")+scale_fill_gradient2(high="red",mid="white",low="green",name="cy_rat")+xlab("Primary tumor")+ylab("Recurring tumor")+ theme(aspect.ratio=1)

ggplot(data=all_means_annot_wide)+geom_point(aes(x=enrichmentCycles.1+0.2,y=value.1,fill="1",col=category),shape=21,size=2.5,alpha=0.6,position=position_jitter(height=0,width=0.1))+geom_point(aes(x=enrichmentCycles.2-0.2,y=value.2,fill="2",col=category),shape=21,size=2.5,alpha=0.6,position=position_jitter(height=0,width=0.1))+facet_wrap(~variable, scales="free_y")+xlab("enrichment Cycles") + ylab("Value")+scale_color_manual(values=c("GBMatch"="grey","GBmatch_val"="black"))

ggplot(all_means_annot_long,aes(x=value,group=surgery.x,col=factor(surgery.x)))+geom_density()+facet_wrap(category~variable,scale="free")

dev.off()
