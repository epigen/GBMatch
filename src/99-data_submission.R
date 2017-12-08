library(project.init)
project.init2("GBMatch")

submission_dir="/data/groups/lab_bock/jklughammer/projects/Glioblastoma_match/forDataSubmission"


out_dir=file.path(submission_dir)
dir.create(file.path(out_dir,"christoph_bock/processed_data"), recursive=TRUE)
dir.create(file.path(submission_dir,"raw_data"), recursive=TRUE)
setwd(out_dir)

annotation=fread(file.path(submission_dir,"Suppl. Table 2 - RRBS summary v5.txt"))

annotation_sub=annotation[,list(raw_data_path=file.path(getOption("PROCESSED.PROJECT"),"results_pipeline",`Sample ID`,"raw",paste0(`Sample ID`,".bam")),processed_data_path=file.path(getOption("PROCESSED.PROJECT"),"results_pipeline",`Sample ID`,"biseq_hg38",paste0("RRBS_cpgMethylation_",`Sample ID`,".bed")),Sample_name=paste0(`Patient ID`,"_",`Sample ID`),title=paste0("RRBS ",`Sample Type`," ",`Patient ID`," Surg.",`Surgery Number`," ",`Material`),source_name=ifelse(`Sample Type`!="Normal white matter control","brain (tumor)","brain (white matter)"),organism="homo sapiens",gender=Sex,IDH=`IDH status`,Center=`Center`,`Sample Type`=`Sample Type`,`Surgery Number`=`Surgery Number`,Multisector=Multisector,Material=Material,molecule="genomic DNA",`Sequencing Instrument`=paste0("Illumina ",`Sequencing Instrument`),`Read Length`=`Read Length`),]


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



