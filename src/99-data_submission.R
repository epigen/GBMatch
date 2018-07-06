library(project.init)
project.init2("GBMatch")

#For progression cohort
submission_dir="/data/groups/lab_bock/jklughammer/projects/Glioblastoma_match/forDataSubmission"
annotation=fread(file.path(submission_dir,"Suppl. Table 2 - RRBS summary v5.txt"))

#For validation cohort
submission_dir="/data/groups/lab_bock/jklughammer/projects/Glioblastoma_match/forDataSubmission_val"
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.2-sample_stats/RRBS_stats.tsv"))
annotation=annotation[Category%in%c("GBmatch_val","GBmatch_valF")]
annotation[,Multisector:="No",]


out_dir=file.path(submission_dir)
dir.create(file.path(out_dir,"christoph_bock/processed_data"), recursive=TRUE)
dir.create(file.path(submission_dir,"raw_data"), recursive=TRUE)
setwd(out_dir)

annotation_sub=annotation[,list(raw_data_path=file.path(getOption("PROCESSED.PROJECT"),"results_pipeline",`Sample ID`,"raw",paste0(`Sample ID`,".bam")),processed_data_path=file.path(getOption("PROCESSED.PROJECT"),"results_pipeline",`Sample ID`,"biseq_hg38",paste0("RRBS_cpgMethylation_",`Sample ID`,".bed")),`Sample ID`=`Sample ID`,Sample_name=paste0(`Patient ID`,"_",`Sample ID`),title=paste0("RRBS ",`Sample Type`," ",`Patient ID`," Surg.",`Surgery Number`," ",`Material`),source_name=ifelse(!`Sample Type`%in%c("Normal white matter control","White matter"),"brain (tumor)","brain (white matter)"),organism="homo sapiens",IDH=`IDH status`,Sex=Sex,Center=`Center`,`Sample Type`=`Sample Type`,`Surgery Number`=`Surgery Number`,Multisector=Multisector,Material=Material,molecule="genomic DNA",`Sequencing Instrument`=paste0("Illumina ",`Sequencing Instrument`),`Read Length`=`Read Length`),]

annotation_sub[,raw_data_file:=paste0(Sample_name,"_unmapped.bam"),]
annotation_sub[,processed_data_file:=paste0(Sample_name,"_cpgmeth.bed"),]

#sort by Sample_name
annotation_sub=annotation_sub[order(Sample_name)]

#calculate processed checksum
annotation_sub[,processed_checksum:=unlist(strsplit(system(paste0(" md5sum ",processed_data_path),intern=TRUE)," "))[1],by=1:nrow(annotation_sub)]
#calculate raw checksum
annotation_sub[,raw_checksum:=unlist(strsplit(system(paste0(" md5sum ",raw_data_path),intern=TRUE)," "))[1],by=1:nrow(annotation_sub)]


#copy processed data to destination
annotation_sub[,system(paste0("cp ",processed_data_path," ",getwd(),"/christoph_bock/processed_data/",processed_data_file)),by=1:nrow(annotation_sub)]

#copy raw data to destination
annotation_sub[,system(paste0("cp ",raw_data_path," ",getwd(),"/raw_data/",raw_data_file)),by=1:nrow(annotation_sub)]

#run EGA preparation of raw data (takes a long time for many files. Run directly on commandline)
system(paste0("cd ", getwd(),"/christoph_bock/processed_data/;java -jar /data/groups/lab_bock/jklughammer/resources/tools_general/EgaCryptor/EgaCryptor.jar -file *.bam"))

#collect encryped file names, and md5 checksums for processed data (EGA)
annotation_sub[,raw_data_file_encr:=paste0(raw_data_file,".gpg"),]
annotation_sub[,raw_checksum_ega:=system(paste0("cat raw_data/",raw_data_file,".md5"),intern=TRUE),by=1:nrow(annotation_sub)]
annotation_sub[,raw_encr_checksum_ega:=system(paste0("cat raw_data/",raw_data_file_encr,".md5"),intern=TRUE),by=1:nrow(annotation_sub)]

write.table(annotation_sub,"annotation_compiled.tsv",quote=FALSE,row.names=FALSE,sep="\t")

#other data types (RNA and WGS)

#WGS
annotation_sub_WGS=annotation[grep("_fv",`Sample ID`),list(raw_data_path=list.files(path=file.path(file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/gb_unmapped_bam/")),pattern=paste0(gsub("_fv","",`Sample ID`),".*[0-9].bam$"),full.names=TRUE,recursive=TRUE),processed_data_path=list.files(path=file.path(getOption("PROCESSED.PROJECT"),"WGS_CNV/cnvkit_vs_pga_hg38/"),pattern=paste0(gsub("_fv","",`Sample ID`),".*call.cns$"),full.names=TRUE),Sample_name=paste0(`Patient ID`,"_",`Sample ID`),title=paste0("WGS ",`Sample Type`," ",`Patient ID`," Surg.",`Surgery Number`," ",`Material`),source_name=ifelse(!`Sample Type`%in%c("Normal white matter control","White matter"),"brain (tumor)","brain (white matter)"),organism="homo sapiens",IDH=`IDH status`,Sex=Sex,Center=`Center`,`Sample Type`=`Sample Type`,`Surgery Number`=`Surgery Number`,Multisector=Multisector,Material=Material,molecule="genomic DNA",`Sequencing Instrument`="Illumina HiSeq_3000",`Read Length`=50),by=`Sample ID`]
annotation_sub_WGS=annotation_sub_WGS[!is.na(raw_data_path)]

annotation_sub_WGS[,raw_data_file:=paste0(Sample_name,"_wgs_unmapped.bam"),]
annotation_sub_WGS[,processed_data_file:=paste0(Sample_name,"_wgs_call.cns"),]

#sort by Sample_name
annotation_sub_WGS=annotation_sub_WGS[order(Sample_name)]

#calculate processed checksum
annotation_sub_WGS[,processed_checksum:=unlist(strsplit(system(paste0(" md5sum ",processed_data_path),intern=TRUE)," "))[1],by=1:nrow(annotation_sub_WGS)]
#calculate raw checksum
annotation_sub_WGS[,raw_checksum:=unlist(strsplit(system(paste0(" md5sum ",raw_data_path),intern=TRUE)," "))[1],by=1:nrow(annotation_sub_WGS)]

#copy processed data to destination
dir.create(file.path(out_dir,"christoph_bock/processed_data_WGS"), recursive=TRUE)
annotation_sub_WGS[,system(paste0("cp ",processed_data_path," ",getwd(),"/christoph_bock/processed_data_WGS/",processed_data_file)),by=1:nrow(annotation_sub_WGS)]

#copy raw data to destination
annotation_sub_WGS[,system(paste0("cp ",raw_data_path," ",getwd(),"/raw_data/",raw_data_file)),by=1:nrow(annotation_sub_WGS)]

#run EGA preparation of raw data (takes a long time for many files. Run directly on commandline)
system(paste0("cd ", getwd(),"/christoph_bock/processed_data/;java -jar /data/groups/lab_bock/jklughammer/resources/tools_general/EgaCryptor/EgaCryptor.jar -file *_wgs_unmapped.bam"))

#collect encryped file names, and md5 checksums for processed data (EGA)
annotation_sub_WGS[,raw_data_file_encr:=paste0(raw_data_file,".gpg"),]
annotation_sub_WGS[,raw_checksum_ega:=system(paste0("cat raw_data/",raw_data_file,".md5"),intern=TRUE),by=1:nrow(annotation_sub_WGS)]
annotation_sub_WGS[,raw_encr_checksum_ega:=system(paste0("cat raw_data/",raw_data_file_encr,".md5"),intern=TRUE),by=1:nrow(annotation_sub_WGS)]

write.table(annotation_sub_WGS,"annotation_compiled_WGS.tsv",quote=FALSE,row.names=FALSE,sep="\t")


#RNA
processed_data_path_RNA=file.path(out_dir,"christoph_bock/processed_data_RNA")
dir.create(processed_data_path_RNA, recursive=TRUE)
annotation_sub_RNA=annotation[grep("_fv",`Sample ID`),list(raw_data_path=list.files(path=file.path(getOption("PROCESSED.PROJECT"),"results_pipeline_rna/"),pattern=paste0(`Sample ID`,".bam$"),recursive=TRUE,full.names=TRUE),processed_data_path=file.path(processed_data_path_RNA,paste0(`Patient ID`,"_",`Sample ID`,"_rna.RPKM")),Sample_name=paste0(`Patient ID`,"_",`Sample ID`),title=paste0("RNA ",`Sample Type`," ",`Patient ID`," Surg.",`Surgery Number`," ",`Material`),source_name=ifelse(!`Sample Type`%in%c("Normal white matter control","White matter"),"brain (tumor)","brain (white matter)"),organism="homo sapiens",IDH=`IDH status`,Sex=Sex,Center=`Center`,`Sample Type`=`Sample Type`,`Surgery Number`=`Surgery Number`,Multisector=Multisector,Material=Material,molecule="RNA",`Sequencing Instrument`="Illumina HiSeq_3000",`Read Length`=50),by=`Sample ID`]
annotation_sub_RNA=annotation_sub_RNA[!is.na(raw_data_path)]

annotation_sub_RNA[,raw_data_file:=paste0(Sample_name,"_rna_unmapped.bam"),]
annotation_sub_RNA[,processed_data_file:=paste0(Sample_name,"_rna.RPKM"),]

#sort by Sample_name
annotation_sub_RNA=annotation_sub_RNA[order(Sample_name)]

#prepare processed data and copy raw data
prep_processed_data=function(sampleAnnot_merged,results_dir,prefix){

for (i in 1:nrow(sampleAnnot_merged)){
  sample_name=sampleAnnot_merged[i]$`Sample ID`
  sampleID=sampleAnnot_merged[i]$Sample_name
  
  #prep processed
  tr=fread(paste0(results_dir,sample_name,"/",prefix,"/bitSeq/",sample_name,".tr"))
  setnames(tr,names(tr),c("ensG","ensT","transcript_length","transcript_length_adj"))
  mean=fread(paste0(results_dir,sample_name,"/",prefix,"/bitSeq/",sample_name,".mean"))
  setnames(mean,names(mean),c("RPKM","RPKM_variance"))
  mean=cbind(tr,mean)
  write.table(mean,sampleAnnot_merged[i]$processed_data_path,sep="\t",quote=FALSE,row.names=FALSE)
}
}
prep_processed_data(sampleAnnot_merged=annotation_sub_RNA,results_dir=file.path(getOption("PROCESSED.PROJECT"),"results_pipeline_rna/"),prefix="bowtie1_hg38_nc")

#copy raw data to destination
annotation_sub_RNA[,system(paste0("cp ",raw_data_path," ",getwd(),"/raw_data/",raw_data_file)),by=1:nrow(annotation_sub_RNA)]

#calculate processed checksum
annotation_sub_RNA[,processed_checksum:=unlist(strsplit(system(paste0(" md5sum ",processed_data_path),intern=TRUE)," "))[1],by=1:nrow(annotation_sub_RNA)]
#calculate raw checksum
annotation_sub_RNA[,raw_checksum:=unlist(strsplit(system(paste0(" md5sum ",raw_data_path),intern=TRUE)," "))[1],by=1:nrow(annotation_sub_RNA)]

#run EGA preparation of raw data (takes a long time for many files. Run directly on commandline)
system(paste0("cd ", getwd(),"/christoph_bock/processed_data/;java -jar /data/groups/lab_bock/jklughammer/resources/tools_general/EgaCryptor/EgaCryptor.jar -file *_rna_unmapped.bam"))

#collect encryped file names, and md5 checksums for processed data (EGA)
annotation_sub_RNA[,raw_data_file_encr:=paste0(raw_data_file,".gpg"),]
annotation_sub_RNA[,raw_checksum_ega:=system(paste0("cat raw_data/",raw_data_file,".md5"),intern=TRUE),by=1:nrow(annotation_sub_RNA)]
annotation_sub_RNA[,raw_encr_checksum_ega:=system(paste0("cat raw_data/",raw_data_file_encr,".md5"),intern=TRUE),by=1:nrow(annotation_sub_RNA)]

write.table(annotation_sub_RNA,"annotation_compiled_RNA.tsv",quote=FALSE,row.names=FALSE,sep="\t")



