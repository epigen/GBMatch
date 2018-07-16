#NOTE: this script was used to create track hubs for the UCSC genome browser
library(project.init)
project.init2("GBMatch")

annotation_all=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))

out_dir="/data/groups/lab_bock/public_html/papers/GBMatch/data/tracks/"
genome="hg38"

dir.create(paste0(out_dir,genome),recursive=TRUE)

#get bigWig files
annotation_all[,bw_file:=dirdata("results_pipeline/", N_number_seq, "/bigwig_", genome, "/RRBS_", N_number_seq, ".bw"),]

#set surgery to 0 when surgery is na
annotation_all[is.na(surgery.x),surgery.x:=0,]

#create track name
annotation_all[,track_name:=paste0(patID,"_",surgery.x),]


##hub.txt
hub="hub GBMatch
shortLabel GBMatch RRBS
longLabel GBMatch RRBS
genomesFile genomes.txt
email jklughammer@cemm.oeaw.ac.at"
cat(hub,file=paste0(out_dir,"hub.txt"))

##genomes.txt
geneomes="genome hg38
position chr11:2167853-2182439
trackDb hg38/RRBS_trackDB.txt"
cat(geneomes,file=paste0(out_dir,"genomes.txt"))


##create trackdb
db_header="track GBMatch
compositeTrack on
subGroup1 category category GBMatch=GBMatch GBmatch_val=GBmatch_val GBmatch_valF=GBmatch_valF GBmatch_add=GBmatch_add GLASS=GLASS GBmatch_rcl=GBmatch_rcl multiselector=multiselector control=control 
subGroup2 surgery surgery 1=1 2=2 3=3 4=4 5=5 6=6 nd=nd
subGroup3 IDH IDH wt=wt mut=mut
shortLabel GBMatch_RRBS
longLabel GBMatch_RRBS
type bigWig
spectrum on
maxHeightPixels 200:50:8
viewLimits 0:1
visibility hide
dimensions dimA=category dimA=surgery dimA=IDH\n\n"

cat(db_header,file=paste0(out_dir,genome,"/RRBS_trackDB.txt"))

for (i in c(1:nrow(annotation_all))){  
  annot_sub=annotation_all[i]
  track=paste0("RRBS_",annot_sub$N_number_seq)
  sg1= annot_sub$category
  sg2= annot_sub$surgery.x
  sg3= annot_sub$IDH
  sl= annot_sub$track_name
  ll= paste0(annot_sub$track_name,"_",annot_sub$N_number_seq)
  bdu= paste0("RRBS_",annot_sub$N_number_seq,".bw")
  
  track=paste0(
    "\ttrack ",track,"\n",
    "\tparent GBMatch on\n",
    "\ttype bigWig\n",
    "\tsubGroups category=",sg1," surgery=",sg2," IDH=",sg3,"\n",
    "\tcolor 115,191,113\n",
    "\tshortLabel ",sl,"\n",
    "\tlongLabel ",ll,"\n",
    "\tbigDataUrl ",bdu,"\n\n"
  )
  cat(track,file=paste0(out_dir,genome,"/RRBS_trackDB.txt"),append=TRUE)
  file.copy(from=annot_sub$bw_file,to=paste0(out_dir,genome),overwrite=FALSE)
}

