#NOTE: Script to handle supplementary data for the website
library(project.init)
project.init2("GBMatch")

annotation_all=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

wd="/data/groups/lab_bock/jklughammer/projects/Glioblastoma_match/forPublication/"
setwd(wd)

copyFiles=function(inDir,outDir,annotation,filter="*"){
  dir.create(outDir,recursive=TRUE)
  files=data.table(file_paths=list.files(path=inDir,full.names=TRUE))
  files[,patID:=substring(file_paths,first=regexpr("pat_[0-9]{3}_1",file_paths),last=regexpr("pat_[0-9]{3}_1",file_paths)+6),]
  toCopy=files[patID%in%annotation$patID&grepl(filter,file_paths)]
  toCopy[,system(command=paste0("cp ", file_paths," ",outDir)),1:nrow(toCopy)]
  system(command=paste0("zip -r ",outDir,".zip ",outDir))
}


#histo raw 
copyFiles("GBMatch_Revision_Slide_Scans/SlideScans_Overview/","GBMatch_val_histo_raw",annotation_all)

#histo seg 
copyFiles("GBMatch_Revision_Slide_Scans/SlideScans_Overview_annot/","GBMatch_val_histo_seg",annotation_all)

#histo lib 
copyFiles("GBMatch_Revision_Slide_Scans/HE_Slide_scan_library/","GBMatch_val_histo_lib",annotation_all)

#MRI raw
copyFiles("segmentation_images_second_batch/","GBMatch_val_MRI_raw",annotation_all,"_T1c")

#MRI seg
copyFiles("segmentation_images_second_batch/","GBMatch_val_MRI_seg",annotation_all,"segmentation.png")




