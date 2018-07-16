#NOTE: This script runs methclone for the analysis of DNA methylation entropy and entropy shifts between primary tumor and recurrency
library(project.init)
project.init2("GBMatch")
library(data.table)
library("Rsamtools")

##run separately for different read thresholds = minimum number of reads to consider a region
#min_reads=40
#min_reads=60
min_reads=20

out_dir=file.path(getOption("PROCESSED.PROJECT"),paste0("results_analysis/06-methclone/data_minReads",min_reads))
dir.create(out_dir,recursive=TRUE)
logdir=file.path(out_dir,"log")
dir.create(logdir,recursive=TRUE)
setwd(getOption("PROCESSED.PROJECT"))
run_methclone=paste0(getOption("PROJECT.DIR"),"/src/99-run_methclone.sh ")

#get the annotation
annotation=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))
annotation[,aligned_bam:=paste0(resolveLink("results_pipeline_bismark"),"/",N_number_seq,"/bismark_hg38/",N_number_seq,".bam"),]
annotation[category=="control",patID:=paste0("ctr_",1:length(category)),]
annotation[,surgery:=ifelse(!is.na(surgery.x),as.character(surgery.x),as.character(resectionNo)),]
annotation[category=="control",surgery:=patID,]
annotation[category=="GBmatch_valF"&is.na(surgery),surgery:=patID,]

#function that actually submits the comparisons
submit=function(job,logdir,run_methclone,bam1,bam2,out_dir,sleep){
  tempdir=file.path(out_dir,job)
  dir.create(tempdir,showWarnings = FALSE)
  bam1_exp=normalizePath(bam1)
  bam2_exp=normalizePath(bam2)
  if (is.na(bam1_exp)|is.na(bam2_exp)){next}
  
  sbatch_command=paste0("sbatch --export=ALL --get-user-env --job-name=methclone_",job," --workdir=",tempdir," --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6000 --partition=shortq --time=12:00:00 -o ",logdir,"/",job,"_%j.log ", paste0(c(run_methclone,bam1_exp,bam2_exp,out_dir,job,min_reads),collapse=" "))
  
  if (!file.exists(file.path(tempdir,paste0(job,".output.txt")))){
    system(sbatch_command)
  }else{
    message(paste0(job," already complete. Not submitting!"))
  }
  return(sbatch_command)
  Sys.sleep(sleep)
}


##Create the test combinations

#compare all samples from one patient (1 vs 2)
test_combinations=annotation[order(surgery)][!is.na(patID)&category%in%c("GBMatch","GBmatch_add","GLASS")&surgery%in%c(1:2),if(length(aligned_bam)>1){combi_files=combn(aligned_bam,2);combi_names=combn(N_number_seq,2);combi_surgery=combn(surgery,2);list(f1=combi_files[1,],f2=combi_files[2,],s1=combi_names[1,],s2=combi_names[2,],surg1=combi_surgery[1,],surg2=combi_surgery[2,])}else{message(paste0("INFO: Only 1 sample for ",patID,".Omitting."))},by="patID"]

#samples vs. controls
controls=annotation[category=="control"&`Unique_CpGs`>1000000]
test_combinations_ctr=annotation[order(surgery)][!is.na(patID)&category%in%c("GBMatch","GBmatch_add","GBmatch_val","GBmatch_valF","GLASS","control")&surgery%in%c(1:2),data.table(patID=patID,f2=aligned_bam,f1=controls$aligned_bam,s2=N_number_seq,s1=controls$N_number_seq,surg1=controls$patID,surg2=surgery),by=1:length(N_number_seq)]
test_combinations_ctr=test_combinations_ctr[surg1!=surg2]    

#samples vs. self to calculate entropy (not entropy shifts)
test_combinations_self=annotation[order(surgery)][!is.na(patID)&surgery%in%c(1:2),data.table(patID=patID,f2=aligned_bam,f1=aligned_bam,s2=N_number_seq,s1=N_number_seq,surg1=surgery,surg2=surgery),by=1:length(N_number_seq)]


#controls vs. contols
test_combinations_ctrVSctr=annotation[order(surgery)][!is.na(patID)&category%in%c("control"),data.table(patID=patID,f2=aligned_bam,f1=controls$aligned_bam,s2=N_number_seq,s1=controls$N_number_seq,surg1=controls$patID,surg2=surgery),by=1:length(N_number_seq)]
test_combinations_ctrVSctr=test_combinations_ctrVSctr[surg1!=surg2]    


##check if s1 is always before s2 (should be empty)
test_combinations[surg1>surg2]

##create job name
test_combinations[,comparisonID:=paste0(patID,"__",s1,"vs",s2,"__",surg1,"vs",surg2),by="patID"]
test_combinations_ctr[,comparisonID:=paste0(patID,"__",s1,"vs",s2,"__",surg1,"vs",surg2)]
test_combinations_self[,comparisonID:=paste0(patID,"__",s1,"vs",s2,"__",surg1,"vs",surg2)]
test_combinations_ctrVSctr[,comparisonID:=paste0(patID,"__",s1,"vs",s2,"__",surg1,"vs",surg2)]


##create commands and submit
### NOTE:if indices have to be created for samples that are used by multiple comparisons indexing will be called multiple times. Only submit one job that creates the index. 

commands=test_combinations[,submit(job=comparisonID,logdir,run_methclone,bam1=f1,bam2=f2,out_dir),by=1:nrow(test_combinations)]

commands_ctr=test_combinations_ctr[,submit(job=comparisonID,logdir,run_methclone,bam1=f1,bam2=f2,out_dir,sleep=30),by=1:nrow(test_combinations_ctr)]
commands_self=test_combinations_self[,submit(job=comparisonID,logdir,run_methclone,bam1=f1,bam2=f2,out_dir,sleep=30),by=1:nrow(test_combinations_self)]

commands_ctrVSctr=test_combinations_ctrVSctr[,submit(job=comparisonID,logdir,run_methclone,bam1=f1,bam2=f2,out_dir,sleep=30),by=1:nrow(test_combinations_ctrVSctr)]

