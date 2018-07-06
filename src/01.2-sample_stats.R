library(project.init)
project.init2("GBMatch")
library(ggmap)
library(maps)
library(ggrepel)

out_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.2-sample_stats/")
dir.create(out_dir)
setwd(out_dir)

annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))


###############################
##plot sequencing/ RRBS stats##
###############################

annotation[,sample_group:=ifelse(N_number_st%in%annotation[material=="RCL"]$N_number_st,N_number_st,NA),]

pdf("rand_frag_RCL.pdf",height=5,width=6)
ggplot(annotation[!is.na(sample_group)&!material=="frozen"],aes(x=material,y=rand_frag_perc,fill=sample_group,group=sample_group))+geom_line(alpha=0.6)+geom_point(alpha=0.5,shape=21,size=3)+ylim(c(0,100))
dev.off()

annotation[,material_cohort:=paste0(cohort,"\n",material),]
annotation[,time_to_prep:=as.numeric(as.Date(flowcell_date)-as.Date(SurgeryDate)),]
cor_time=annotation[,list(date=min(SurgeryDate),y=rand_frag_perc[which.min(as.Date(SurgeryDate))],time_cor=cor(rand_frag_perc,-time_to_prep),N=nrow(na.omit(cbind(rand_frag_perc,-time_to_prep)))),by=c("material_cohort")]

####Figure S1C start #################################
pdf("rand_fragVStime.pdf",height=3.5,width=4.5)
ggplot(annotation,aes(x=as.Date(SurgeryDate),y=rand_frag_perc,fill=material_cohort,group=material_cohort,col=material_cohort))+geom_point(alpha=0.5,shape=21,size=3)+geom_text(hjust=0,data=cor_time,aes(x=as.Date(date),y=y,label=paste0(material_cohort,": r=",round(time_cor,3)," (N=",N,")")))+ylim(c(0,100))+geom_smooth(method="lm") +ylab("% randomly fragmented reads")+xlab("Surgery date")
dev.off()

#use FFPE necrosis values because they are not available for frozen
annotation[,`Relative share necrosis`:=ifelse(is.na(`Relative share necrosis`),`Relative share necrosis`[material=="FFPE"],`Relative share necrosis`),by="N_number_st"]

cor_nec=annotation[,list(`Relative share necrosis`=max(`Relative share necrosis`,na.rm=TRUE),rand_frag_perc=rand_frag_perc[which.max(`Relative share necrosis`)],cor=cor(`Relative share necrosis`,rand_frag_perc,use="pairwise.complete"),N=nrow(na.omit(cbind(`Relative share necrosis`,rand_frag_perc)))),by="material_cohort"]
pdf("rand_fragVSnecrosis.pdf",height=3.5,width=4.5)
ggplot(annotation,aes(x=`Relative share necrosis`,y=rand_frag_perc,fill=material_cohort,group=material_cohort,col=material_cohort))+geom_point(alpha=0.5,shape=21,size=3)+geom_text(data=cor_nec,hjust=1,aes(label=paste0(material_cohort,": r=",round(cor,3)," (N=",N,")")))+ylim(c(0,100))+geom_smooth(method="lm")+ylab("% randomly fragmented reads")+xlab("Relative share necrosis")
dev.off()

cor_CpG=annotation[,list(Raw_reads=max(Raw_reads),Unique_CpGs=Unique_CpGs[which.max(Raw_reads)],cor=cor(Raw_reads,Unique_CpGs),N=nrow(na.omit(cbind(Raw_reads,Unique_CpGs)))),by="material_cohort"]
pdf("readsVSCpG.pdf",height=3.5,width=4.5)
ggplot(annotation,aes(x=Raw_reads/1000000,y=Unique_CpGs/1000000,fill=material_cohort,group=material_cohort,col=material_cohort))+geom_point(alpha=0.5,shape=21,size=3)+geom_text(data=cor_CpG,hjust=1,aes(label=paste0(material_cohort,": r=",round(cor,3)," (N=",N,")")))+geom_smooth(method="lm")+ylab("Unique CpGs (Million)")+xlab("Raw reads (Million)")
dev.off()
####Figure S1C start #################################


#bisylfite conversion
####Figure S1D
pdf("bisulfite_conversion.pdf",width=4.5,height=2.5)
ggplot(annotation)+geom_point(aes(x=material,y=K1_unmethylated_meth*100,col="Unmethylated",fill="Unmethylated"),position=position_jitter(width=0.3,height=0),alpha=0.5,shape=21)+geom_point(aes(x=material,y=K3_methylated_meth*100,col="Methylated",fill="Methylated"),position=position_jitter(width=0.3,height=0),alpha=0.5,shape=21)+geom_hline(yintercept=5,color="grey",lty=20)+geom_hline(yintercept=95,color="grey",lty=20)+facet_grid(~cohort,space="free_x",scale="free_x")+ylab("Control methylation (%)")
dev.off()

########################
##Supplementary tables##
########################

#Supplementary Table 1 (patients)
pat_stats=annotation[,list(N_number_1st=N_number_seq[order(category,surgery.x)][1],AgeAtDiag=min(Age),categories=paste0(unique(category),collapse=","),cohorts=paste0(unique(cohort),collapse=","),materials=paste0(unique(material),collapse=","),WHO2016_classification=paste0(unique(WHO2016_classification[order(surgery.x)]),collapse=","),Nsamples=.N),by=c("patID","Center","Sex","IDH","NoOfSurgeries","timeToFirstProg","Follow-up_years", "VitalStatus","StuppComplete")]
setcolorder(pat_stats,c("N_number_1st","patID","Center","Sex","AgeAtDiag","IDH","WHO2016_classification","NoOfSurgeries","Nsamples","timeToFirstProg","Follow-up_years", "VitalStatus","StuppComplete","cohorts","categories","materials"))
setnames(pat_stats,names(pat_stats),c("Sample ID","Patient ID","Center","Sex","Age (Diagnosis)","IDH status","WHO2016 Classifications","Number of Surgeries","Number of Samples","Time to Progression (m)","Follow-up Time (m)","Status","Stupp Complete","Cohorts","Categories","Materials"))
pat_stats[,`Time to Progression (m)`:=round(`Time to Progression (m)`/30.4,1),]
pat_stats[,`Follow-up Time (m)`:=round(`Follow-up Time (m)`*12,1),]
write.table(pat_stats,file="Pat_stats.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#calculate basic pat stats
#GBMatch (primary cohort)
pat_stats_summary_GBMatch=pat_stats[grepl("GBMatch",Categories),list(mean_ttp=mean(`Time to Progression (m)`),mean_followup=mean(`Follow-up Time (m)`),sd_ttp=mean(`Time to Progression (m)`),sd_followup=mean(`Follow-up Time (m)`),median_ttp=median(`Time to Progression (m)`),median_followup=median(`Follow-up Time (m)`),male=sum(Sex=="m"),female=sum(Sex=="f"),median_age=as.double(median(`Age (Diagnosis)`)),mean_age=mean(`Age (Diagnosis)`),sd_age=sd(`Age (Diagnosis)`)),by="IDH status"]
write.table(pat_stats_summary_GBMatch,file="Pat_stats_summary_GBMatch.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#GBmatch_val (vaidation cohort)
pat_stats_summary_GBmatch_val=pat_stats[grepl("GBmatch_val([^F]|$)",Categories,perl=TRUE),list(mean_ttp=mean(`Time to Progression (m)`,na.rm=TRUE),mean_followup=mean(`Follow-up Time (m)`,na.rm=TRUE),sd_ttp=mean(`Time to Progression (m)`,na.rm=TRUE),sd_followup=mean(`Follow-up Time (m)`,na.rm=TRUE),median_ttp=median(`Time to Progression (m)`,na.rm=TRUE),median_followup=median(`Follow-up Time (m)`,na.rm=TRUE),male=sum(Sex=="m"),female=sum(Sex=="f"),median_age=as.double(median(`Age (Diagnosis)`)),mean_age=mean(`Age (Diagnosis)`),sd_age=sd(`Age (Diagnosis)`)),by="IDH status"]
write.table(pat_stats_summary_GBmatch_val,file="Pat_stats_summary_GBMatch_val.tsv",sep="\t",quote=FALSE,row.names=FALSE)


#Supplementary Table 2 (RRBS stats)
rrbs_stats=annotation[,c("N_number_seq","patID","Sex","Center","tissue","IDH","surgery.x","material","category","inputDNA","adapter","enrichmentCycles","instrument_model","read_length","read_type", "library","Raw_reads","Aligned_reads","Alignment_rate","rand_frag_perc","Unique_CpGs","meanCoverage","bisulfiteConversionRate","K1_unmethylated_meth","K3_methylated_meth","qualTier"),with=FALSE]
setnames(rrbs_stats,names(rrbs_stats),c("Sample ID","Patient ID","Sex","Center","Sample Type","IDH status","Surgery Number","Material","Category","Input DNA [ng]","Adapter","Enrichment Cycles","Sequencing Instrument","Read Length","Read Type","Library","Raw Reads","Aligned Reads","Alignment Rate [%]","Fragmented Reads [%]","Unique GpGs","Mean Coverage","Non-CpG Conversion Rate","DNA methylation (unmethylated spike-in)","DNA methylation (methylated spike-in)","Quality Tier"))
write.table(rrbs_stats,file="RRBS_stats.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#some summary stats for paper
write.table(rrbs_stats[,list(`Median CpGs`=median(`Unique GpGs`)),by="Material"],file="stats.txt",sep="\t",quote=FALSE,row.names=FALSE)
cat(paste0("\n\nMin conv rate:\n",1-rrbs_stats[,max(`DNA methylation (unmethylated spike-in)`)]),file="stats.txt",append=TRUE)
cat(paste0("\n\nFFPE samples <500000 GpGs:\n",nrow(rrbs_stats[Material=="FFPE"&`Unique GpGs`<500000])," (",nrow(rrbs_stats[Material=="FFPE"&`Unique GpGs`<500000])/nrow(rrbs_stats[Material=="FFPE"]),")"),file="stats.txt",append=TRUE)
sub=pat_stats[grepl("FFPE",Materials)]
cat(paste0("\n\nFFPE patients\nIDHwt patients: ",nrow(sub[`IDH status`=="wt"]),"\nIDHmut patients: ",nrow(sub[`IDH status`=="mut"]),"\nControl patients: ",nrow(sub[Categories=="control"]),"\nAll patients: ",nrow(sub)),file="stats.txt",append=TRUE)
sub=rrbs_stats[Material=="FFPE"]
cat(paste0("\n\nFFPE samples\nIDHwt patients: ",nrow(sub[`Patient ID`%in%pat_stats[`IDH status`=="wt"]$`Patient ID`]),"\nIDHmut patients: ",nrow(sub[`Patient ID`%in%pat_stats[`IDH status`=="mut"]$`Patient ID`]),"\nControl patients: ",nrow(sub[Category=="control"]),"\nAll patients: ",nrow(sub)),file="stats.txt",append=TRUE)


#########################
##Cohort overview plots##
#########################

#Cohort over time
####Figure 1B and S1B
clinical_annot_long_date=melt(annotation[category%in%c("GBMatch","GBmatch_val")&IDH=="wt"&!is.na(VitalStatus)],measure.vars =c("DateOfBirth","DateOfDeath_LastFollow-up","ProgressionDate","SurgeryDate","TreatmentDate"),variable.name="Event",value.name = "Date_form")
clinical_annot_long_date[Event!="DateOfDeath_LastFollow-up",VitalStatus:="alive",]
clinical_annot_long_date[Event=="SurgeryDate",Event:=paste0("Surgery ",surgery.x),]
clinical_annot_long_date[Event=="ProgressionDate",Event:=paste0("Progression ",surgery.x-1),]
clinical_annot_long_date[,Date_form:=as.Date(Date_form),]

eventOrder=eventOrder=c("DateOfBirth",as.character(unique(clinical_annot_long_date[,Event[order(surgery)],])[-c(1:2)]),"DateOfDeath_LastFollow-up")
clinical_annot_long_date[,Event:=factor(Event,levels=eventOrder),]

myOrder=as.character(clinical_annot_long_date[Event=="Surgery 1"][order(Date_form)]$patID)
clinical_annot_long_date[,patID:=factor(patID,levels=unique(myOrder)),]

pl=ggplot(clinical_annot_long_date[Event!="DateOfBirth"],aes(x=Date_form,y=patID,col=factor(surgery.x)))+
  geom_line()+
  geom_point(data=clinical_annot_long_date[grepl("DateOfDeath_LastFollow-up",Event)],aes(x=Date_form,y=patID,fill=factor(VitalStatus)),col="black",shape=22,size=3)+
  geom_point(data=clinical_annot_long_date[grepl("Progression",Event)],aes(x=Date_form,y=patID),fill="white",shape=21,size=4)+
  geom_point(data=clinical_annot_long_date[grepl("Surgery",Event)],aes(x=Date_form,y=patID),shape=3,size=2)+
  scale_fill_manual(values=c("white","black","grey","grey"))+
  scale_color_manual(values=c("#e41a1c", "#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf","black"))+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))+
  coord_flip()+facet_grid(~category,scales="free",space="free")

pdf("timeline_prog_surg.pdf",height=5,width=28)
print(pl)
dev.off()

pl2=pl+geom_point(data=clinical_annot_long_date[grepl("Surgery",Event)&!is.na(N_number_seq)],aes(x=Date_form,y=patID),shape=16,size=4)

pdf("timeline_Nno.pdf",height=8,width=28)
print(pl2)
dev.off()

dataFig1b=clinical_annot_long_date[Event!="DateOfBirth"&!is.na(Date_form)&category=="GBMatch"][,list(Date_form,patID,surgery.x,VitalStatus,Event,category)]
write.table(dataFig1b,"Source Data Figure 1b.csv",sep=";",quote=FALSE,row.names=FALSE)


#single patient
####Figure 1A
patient="pat_014"
sub=clinical_annot_long_date[patID==patient&Event!="DateOfBirth"]

pl=ggplot(sub,aes(x=Date_form,y=patID,col=factor(surgery.x)))+
  geom_line()+
  geom_point(data=sub[grepl("DateOfDeath_LastFollow-up",Event)],aes(x=Date_form,y=patID,fill=factor(VitalStatus)),col="black",shape=22,size=3)+
  geom_point(data=sub[grepl("Progression",Event)],aes(x=Date_form,y=patID),fill="white",shape=21,size=4)+
  geom_point(data=sub[grepl("Surgery",Event)],aes(x=Date_form,y=patID),shape=3,size=2)+
  geom_text(data=sub[grepl("Treatment",Event)],aes(x=Date_form,y=1.1,label=CTXTreatment),col="black")+
  geom_text(data=sub[grepl("Treatment",Event)],aes(x=Date_form,y=0.9,label=RTXTreatment),col="black")+
  scale_fill_manual(values=c(alive="white",dead="black",lost="grey"))+
  scale_color_manual(values=c("#e41a1c", "#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf"))+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))
pdf(paste0("timeline_",patient,".pdf"),height=3,width=7)
print(pl)
dev.off()


#plot locations
####Figure S1A
mapWorld <- borders("world","austria", colour="gray50", fill="white") # create a layer of borders

centers=pat_stats[grepl("GBMatch|GBmatch_add|GBmatch_val([^F]|$)",Categories),list(N=.N,IDH_count=paste0(length(`IDH status`[`IDH status`=="wt"]),"/",IDHwt=length(`IDH status`[`IDH status`=="mut"]))),by="Center"]
centers[Center=="Rudolfstiftung",Center:="Rudolfstiftung Vienna",]
centers[,Center_simpl:=gsub("MedUni |Rudolfstiftung ","",Center),]
centers=cbind(centers,as.data.table(geocode(centers$Center_simpl)))
while(any(is.na(c(centers$lon,centers$lat)))){
try(centers[is.na(lon),c("lon","lat"):=try(geocode(Center_simpl),silent=TRUE),],silent=TRUE)
}
centers


pdf(paste0("centers_locations.pdf"),height=3,width=7)
ggplot(centers)+mapWorld+geom_point(aes(x=as.numeric(lon),y=as.numeric(lat),size=N,col=N),alpha=0.6)+geom_text(aes(x=as.numeric(lon),y=as.numeric(lat),label=Center))+geom_text(aes(x=as.numeric(lon),y=as.numeric(lat),label=IDH_count))+scale_size(range=c(6,20))+scale_color_gradient(low="blue",high="red")+theme_bw()
dev.off()

#plot succession of treatments
sub=annotation[category=="GBMatch"&treatment>0]
sub[,plot_CTX:=ifelse(grepl("Stupp|TMZ",CTXTreatment),"Stupp/TMZ",ifelse(grepl("No therapy|none|Watchful",CTXTreatment),"No therapy",ifelse(grepl("Unknown",CTXTreatment),"Unknown",CTXTreatment))),]

pdf("ctx_treatments.pdf",height=5,width=8)
ggplot(sub,aes(x=as.factor(treatment),fill=plot_CTX,group=plot_CTX))+geom_bar(col="black",position="fill")+scale_fill_manual(values=c(terrain.colors(9, alpha = 1)[1:8],"light grey"))
dev.off()

