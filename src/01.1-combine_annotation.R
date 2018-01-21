library(project.init)
project.init2("GBMatch")
source(file.path(getOption("PROJECT.DIR"),"src/99-prep_annotations.R"))

wd=file.path(getOption("PROJECT.DIR"),"metadata/sample_annotations")
setwd(wd)

our_dir=file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation")
dir.create(our_dir,recursive=TRUE)

#psa=fread("ProjectSampleAnnotation_add.txt")
psa=fread("GBMatch_projectSampleAnnotation_rev.csv")

#histo_all=prep_histo("GBMatch_histo part 3_20022017.csv")
histo_all=prep_histo("GBMatch_histo_immuno_rev.csv","GBMatch_histo_immuno_EZH2_rev.csv")
histo=histo_all$histo_annotation
histo_description=histo_all$histo_description

sequencing_stats=fread(file.path(getOption("PROCESSED.PROJECT"),"results_pipeline/RRBS_stats_summary.tsv"))
setnames(sequencing_stats,"sampleName","N_number")
sequencing_annot=fread(file.path(getOption("PROJECT.DIR"),"metadata/GBMatch_samples.csv"))
setnames(sequencing_annot,"sample_name","N_number")
sequencing=merge(sequencing_annot,sequencing_stats,by="N_number")

#clinical=prep_clin_annot("GBMatch_clinical data_20032017_mod.csv");clinical
clinical=prep_clin_annot("GBMatch_clinical data_2018_rev.csv");clinical

imaging=fread("GBMatch_imaging_04012016.tsv")
imaging_prog=fread("GBMatch_MR_progressiontype.csv")

#imaging_auto=fread("TumorFeaturesFixed.tsv")
imaging_auto=fread("GBMatch_MR_auto_rev.csv")
setnames(imaging_auto,"N-Nummer","N_number")

#prep=fread("GBMatch_libraryPrep.csv")
prep=fread("GBMatch_libraryPrep_rev.csv")

segmentation=fread("GBMatch_histo_part4_02142017.tsv")

#include fragmentation check
############bash######################################
#Run this command in results_pipeline
#for bam in `ls */raw/*.bam`; do printf "$bam\t"; printf "`samtools view $bam |grep -r -P "\t[CT]GG" |wc -l`\t";printf "`samtools view $bam |grep -r -P "\t[CT]GA" |wc -l`\t";printf "`samtools view -c $bam`\n"; done > fragmentation_check.tsv
######################################################
frag_check=fread(file.path(getOption("PROCESSED.PROJECT"),"results_pipeline/fragmentation_check.tsv"))
setnames(frag_check,names(frag_check),c("bam","mspI","taqI","all"))
frag_check[,N_number:=gsub(".bam","",gsub(".*/","",bam)),]
frag_check[,rand_frag:=all-(mspI+taqI),]
frag_check[,rand_frag_perc:=rand_frag/all*100,]
sequencing=merge(sequencing,frag_check,by="N_number",all.x=TRUE)


psa[,N_number_st:=gsub("\\.|_|-| |/|A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
sequencing[,N_number_st:=gsub("\\.|_|-| |rr|A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
sequencing[,N_number_st:=gsub("^M","N",N_number_st, perl=TRUE),]
sequencing[,N_number_st:=gsub("14v","14",N_number_st, perl=TRUE),]
prep[,N_number_st:=gsub("\\.|_|-|/| |rr|A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
prep[,N_number_st:=gsub("^M","N",N_number_st, perl=TRUE),]
prep[,N_number_st:=gsub("14v","14",N_number_st, perl=TRUE),]
histo[,N_number_st:=gsub("\\.|_|-| |A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
histo_description[,N_number_st:=gsub("\\.|_|-| |A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
clinical[,N_number_st:=gsub("\\.|_|-| |A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
imaging[,N_number_st:=gsub("_|-| |A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
imaging_auto[,N_number_st:=gsub("\\.|_|-| |A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number, perl=TRUE),]
imaging_prog[,patientID:=gsub(" ","",patientID, perl=TRUE),]
#for now remove entries without N_number (progression without surgery)
clinical=clinical[!is.na(N_number_st)]
segmentation[,N_number_seg:=gsub(" .*","",FileName)]
segmentation[,N_number_st:=gsub("\\.|_|-| |A|B|C|D|E|F|G1|G2|G3|Z2|Z","",N_number_seg, perl=TRUE)]


#remove GPO samples
sequencing=sequencing[!grepl("GPO",N_number_st)]

#change TCGA to GLASS
psa[category=="TCGA",category:="GLASS"]

#add cohort category
psa[,cohort:=ifelse(category%in%c("GBmatch_val","GBmatch_valF"),"validation","primary"),]

#remove N669-06A because not sequenced
histo=histo[N_number!="N669-06A"]
#histo_description=histo_description[N_number!="N669-06A"] #needed for multiselector annoatation
segmentation=segmentation[N_number_seg!="N669-06A"]

#remove N1267-07 from prep table because only repeat is seqeunced
prep=prep[N_number!="N1267-07"]

#imaging translations
imaging[,border:=c(NA,"sharp","diffuse")[factor(border)],]
imaging[,contrast:=c(NA,"necrotic","cystic","solid","mixed solid/necrotic","mixed necrotic/solid")[factor(contrast)],]

#segmentation: -1 -->NA
segmentation[segmentation==-1]=NA
#sort out duplicate samples in segmentation
seg_match=merge(segmentation[N_number_st%in%N_number_st[duplicated(N_number_st)],c("N_number_st","N_number_seg"),with=FALSE],sequencing[,c("N_number_st","N_number"),with=FALSE],by="N_number_st")
seg_match[,N_number_seg_st:=gsub("_|-| |","",N_number_seg, perl=TRUE),]
seg_match[,N_number_seq_st:=gsub("_|-| |","",N_number, perl=TRUE),]
seg_match[,keep:=ifelse(N_number_seg_st==N_number_seq_st,TRUE,FALSE),]
segmentation=rbindlist(list(segmentation[N_number_seg%in%seg_match[keep==TRUE]$N_number_seg],segmentation[!N_number_st%in%N_number_st[duplicated(N_number_st)]]))

#change N_number_st in segmentation because N_number is off by 1 (2 blocks from same operation)
segmentation[N_number_seg=="N370-08",N_number_st:="N36908",]
segmentation[N_number_seg=="N694-15",N_number_st:="N69315",]

#keep 986 samples separate 
histo[grep("N986",N_number),N_number_st:=gsub("_|-| ","",N_number)]
histo_description[grep("N986",N_number),N_number_st:=gsub("_|-| ","",N_number)]
psa[grep("N986",N_number),N_number_st:=gsub("_|-| ","",N_number)]
prep[grep("N986",N_number),N_number_st:=gsub("_|-| ","",N_number)]
sequencing[grep("N986",N_number),N_number_st:=gsub("_|-| ","",N_number)]
segmentation[grep("N986",N_number_seg),N_number_st:=gsub("_|-| ","",N_number_seg)]

#match rk samples #
psa[N_number=="N173-14_rk",N_number_st:="N17414rk",]
prep[N_number=="N173-14_rk",N_number_st:="N17414rk",]
sequencing[N_number=="N173_14_rk",N_number_st:="N17414rk",]
segmentation[N_number_seg=="N173-14_rk",N_number_st:="N17414rk",]
psa[N_number=="N656-15_rk",N_number_st:="N65815rk",]
prep[N_number=="N656-15_rk",N_number_st:="N65815rk",]
sequencing[N_number=="N656_15_rk",N_number_st:="N65815rk",]
segmentation[N_number_seg=="N656-15_rk",N_number_st:="N65815rk",]


setnames(psa,"N_number","N_number_psa")
setnames(clinical,"N_number","N_number_clin")
setnames(imaging,"N_number","N_number_img")
setnames(imaging_auto,"N_number","N_number_imgA")
setnames(histo,"N_number","N_number_hist")
setnames(histo_description,"N_number","N_number_hDesc")
setnames(sequencing,"N_number","N_number_seq")
setnames(prep,"N_number","N_number_prep")


#check for duplicated names (if they are actually the same samples it's ok for them to have dupl. names)
stopifnot(nrow(psa[N_number_st%in%psa[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(clinical[N_number_st%in%clinical[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(imaging[N_number_st%in%imaging[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(imaging_auto[N_number_st%in%imaging[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(histo[N_number_st%in%histo[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(histo_description[N_number_st%in%histo_description[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(sequencing[N_number_st%in%sequencing[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(prep[N_number_st%in%prep[duplicated(N_number_st)]$N_number_st])==0)
stopifnot(nrow(segmentation[N_number_st%in%segmentation[duplicated(N_number_st)]$N_number_st])==0)


column_annotation=list(psa=names(psa),histo_immuno=names(histo),histo_classification=names(histo_description),histo_segmentation=names(segmentation),sequencing_stats=c(names(sequencing_stats),"qualTier"),sequencing_annot=names(sequencing_annot),library_prep=names(prep),clinical_annotation=c(names(clinical),"timeToFirstProg","timeToSecSurg"),imaging=names(imaging),imaging_progression=names(imaging_prog),imaging_segmentation=names(imaging_auto),random_fragmentation=names(frag_check))
save(column_annotation,file = paste0(our_dir,"/column_annoation.RData"))


merge1=merge(prep,sequencing , all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
merge1[,N_number_st:=gsub("tr3","",N_number_st, perl=TRUE),]
merge1.1=merge(psa,merge1, all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
merge1.1=merge(merge1.1,segmentation,all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
#adjust names of glass samples in psa sothat they find their ffpe counterparts
merge1.1[N_number_psa=="N1131-12_f",N_number_st:="N113212f",]
merge1.1[N_number_psa=="N1154/55-13_f",N_number_st:="N115413f",]
merge1.1[N_number_psa=="N1592-95-12_f",N_number_st:="N159312f",]
merge1.1[N_number_psa=="N1575-78-12_f",N_number_st:="N157612f",]
#adjust names of frozen revision samples in psa sothat they find their ffpe counterparts
merge1.1[N_number_psa=="N1164-14_fv",N_number_st:="N116514fv",]
merge1.1[N_number_psa=="N1196-13_fv",N_number_st:="N119513fv",]
merge1.1[N_number_psa=="N1917-12_fv",N_number_st:="N191612fv",]
merge1.1[N_number_psa=="N1934-14_fv",N_number_st:="N193514fv",]
merge1.1[N_number_psa=="N486-13_fv",N_number_st:="N49013fv",]
merge1.1[N_number_psa=="N724-16_fv",N_number_st:="N72716fv",]
merge1.1[N_number_psa=="N838-13_fv",N_number_st:="N84513fv",]
merge1.1[N_number_psa=="N879-12_fv",N_number_st:="N88612fv",]

merge1.1[,N_number_st:=gsub("rk","",N_number_st, perl=TRUE)]
merge1.1[,N_number_st:=gsub("fv","",N_number_st, perl=TRUE)]
merge1.1[,N_number_st:=gsub("f","",N_number_st, perl=TRUE)]

merge2=merge(merge1.1,histo , all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
merge2[,N_number_st:=gsub("I|II|III","",N_number_st, perl=TRUE)]
merge2.1=merge(merge2,histo_description , all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
merge2.1[,N_number_st:=gsub("A|B|C|D|E","",N_number_st, perl=TRUE)] #because of N986 samples
merge2.1[,N_number_st:=gsub("^N66906$","N66806",N_number_st, perl=TRUE)]
merge3=merge(clinical,imaging , all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
merge3.1=merge(merge3,imaging_auto , all=TRUE,by="N_number_st",suffixes = c(".x",".y"))
merge3.2=merge(merge2.1,merge3.1 , all=TRUE,by="N_number_st",suffixes = c(".x",".y"))

#add imaging_prog via patientID
merge4=merge(merge3.2,imaging_prog,by="patientID",all=TRUE)


#quality annotate samples based on the number or unique CpGs
merge4[,qualTier:=ifelse(Unique_CpGs>3000000,1,ifelse(Unique_CpGs>2000000,2,ifelse(Unique_CpGs>1000000,3,4))),]

#missing control
#for now exclude GBmatch_add because these samples have no annotation (yet)
merge4[is.na(category),category:="Missing"]
sub=merge4[category=="GBmatch_add"]
sub=merge4[category=="GBmatch_rcl"]
sub=merge4[category=="GBMatch"]
sub=merge4[category=="GBmatch_valF"]
sub=merge4[category=="GBmatch_val"]
sub=merge4[category=="multiselector"]
sub=merge4[category=="control"]
sub=merge4[category=="Missing"]

sub[is.na(N_number_psa)|is.na(N_number_seq)|is.na(N_number_seg)|is.na(N_number_hist)|is.na(N_number_hDesc)|is.na(N_number_clin)|is.na(N_number_img)|is.na(N_number_prep),c(grep("N_number_",names(merge4),value=TRUE),"category"),with=FALSE]

write.table(merge4[,c(grep("N_number_",names(merge4),value=TRUE),"category"),with=FALSE],paste0(our_dir,"/combine_check_table.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

#check N_number_imgA separately because it is very sparse
#merge4[category!="GBmatch_add"][!is.na(N_number_imgA),c(grep("N_number_",names(merge4),value=TRUE),"category"),with=FALSE]


#remove all entries without N_number_seq (samples included in the annotation but not prepped because sample was used up/ missing)
merge4=merge4[!is.na(N_number_seq)]

#remove sample N591_10 because it is a recurrency in the validation cohort (primary tumor was oparated on in Turkey)
merge4=merge4[N_number_seq!="N591_10"]


#preliminarily use "individual" as patientID and patID
message(paste0("Missing patID:\n",paste0(merge4[is.na(patID)]$N_number_seq,collapse=",")))
merge4[is.na(patID),patID:=paste0("ind_",individual),]
merge4[is.na(patientID),patientID:=paste0("ind_",individual),]


#preliminarily use "date" as "SurgeryDate"
message(paste0("Missing SurgeryDate:\n",paste0(merge4[is.na(SurgeryDate)]$N_number_seq,collapse=",")))
merge4[is.na(SurgeryDate),SurgeryDate:=as.Date(date,'%d.%m.%Y'),]
#preliminarily use "resectionNo" as surgery.x
message(paste0("Missing Surgery number:\n",paste0(merge4[is.na(surgery.x)]$N_number_seq,collapse=",")))
merge4[is.na(surgery.x),surgery.x:=as.double(resectionNo),]

#fix center naming (st. poelten)
merge4[,Center:=gsub("\\?","oe",Center),]

#translate site of surgery and side
translations=list(Siteofsurgery=c("0"=NA, "1"="frontal","2"="temporal", "3"="parietal", "4"="occipital", "5"="central","6"="limbic"),Side=c("1"="right", "2"="left", "3"="median"))

for (i in 1:length(translations)){
  name=names(translations[i])
  translation=translations[[i]]
  print(name)
  if (!(name%in%names(merge4))){message("missing");next}
  for (j in 1:length(translation)){
    target=names(translation[j])
    value=  translation[j]
    merge4[,eval(name):=as.character(get(name)),]
    merge4[get(name)==target,eval(name):=value,]
  } 
}
merge4

#calculate time to first progression
merge4[,timeToFirstProg:=as.numeric(unique(as.Date(na.omit(c(ProgressionDate[surgery.x==2],FirstProgressionDate))))-as.Date(unique(SurgeryDate[surgery.x==1]))),by=patID]
merge4[,timeToSecSurg:=as.numeric(as.Date(unique(SurgeryDate[surgery.x==2]))-as.Date(unique(SurgeryDate[surgery.x==1]))),by=patID]

#now make sure that there is only one entry per sample and category (some samples were sequenced more thatn once)
dubl=merge4[category!="multiselector",list(.N,N_number_seq),by=c("category","patID","surgery.x")][N>1]
merge4[N_number_seq%in%dubl$N_number_seq]$N_number_seq
#remove "tr3" samples
merge4=merge4[grep("_tr3$",N_number_seq,invert=TRUE)]

#Change IDH status from samples that show G-CIMP methylation pattern (previous analysis)
#candidate samples: "N706_15","N536_07","N554_03B","N1814_09A","N127_05A","N307_08C"
#validated samples: "N536_07","N554_03B"
#exclude nevertheless??? N706_15, N307_08C

IDH_pat=merge4[N_number_seq%in%c("N536_07","N554_03B","N706_15","N307_08C")]$patID
merge4[patID%in%IDH_pat,list(patID,N_number_seq,surgery.x)]
merge4[patID%in%IDH_pat,genotype:="mtIDH1"]
merge4[patID%in%IDH_pat,IDH:="mut"]

#Kick out secundary GBM
merge4=merge4[phenotype!="Secondary GBM"]

#remove patients with no relapse
surg_check= merge4[,list(max=max(surgery.x),min=min(surgery.x),categories=paste0(category,collapse=",")),by=patID]
surg_check[min!=1|max<2]
merge4=merge4[!patID%in%surg_check[(min!=1|max<2)&categories=="GBMatch"]$patID]

#complement complement Sex with sex
merge4[is.na(Sex),Sex:=sex]


write.table(merge4,paste0(our_dir,"/annotation_combined.tsv"),row.names = FALSE,sep="\t",quote=FALSE)


