#NOTE: helper script containing functions to process and make available/usable some types of annotations.
library(data.table)
library(ggplot2)


#prepare immunohistochemistry annotations
prep_histo=function(input1,input2){
  histo=fread(input1)
  histo_add=fread(input2)
  histo=merge(histo,histo_add,by="N_number",all.x=TRUE)
  histo[histo==""]=NA
  histo=histo[!is.na(N_number)]
  #extract dercription columns
  description_cols=c("WHO2016_classification","WHO2016_classification_comment","Pseudopalisading necrosis","Radiation necrosis","Vascular fibrosis","Vascular proliferation","Calcification","Blood remnants")
  #check which labels are actually there and adapt translations accordingly
  message("Labels found: ")
  labels=apply(histo[,c(description_cols),with=FALSE],2,function(x){sort(unique(x))})
  print(labels)
  #Prepare translations (number to word) the order is what makes the correct match!!
  #1 GBM classic,2 GBM giant cell,3 Gliosarcoma,4 Anaplastic Oligodendroglioma,5 GBM epithelioid,6 Radiation necrosis,7 Diffuse astro II,8 Oligo II
  #0 None,1 Spindle-celled,2 Small-celled,3 Adenoid,4 Clear-celled,5 Giant cell,6 Reactive gliosis,7 Gemistocytic
  translations=list(`WHO2016_classification`=c("GBM classic","GBM giant cell","Gliosarcoma","Anaplastic Oligodendroglioma","Radiation necrosis","Diffuse astro II","Oligo II"),`WHO2016_classification_comment`=c("None","Spindle-celled","Small-celled","Adenoid","Clear-celled","Giant cell","Reactive gliosis","Gemistocytic"),`Pseudopalisading necrosis`=c("Absent","Present","Abundant"),`Radiation necrosis`=c("Absent","Present"),`Vascular fibrosis`=c("Absent","Present","Abundant"),`Vascular proliferation`=c("Absent","Present","Abundant"),`Calcification`=c("Absent","Present","Abundant"),`Blood remnants`=c("Absent","Present","Abundant"))
  
  message("Translations used:" )
  print(translations)
  
  histo_description=histo[,c("N_number",description_cols),with=FALSE][,lapply(.SD, as.character), by=N_number]
  #keep first entry in the comment column
  histo_description[,WHO2016_classification_comment:=gsub(",.*","",WHO2016_classification_comment),]
  
  #actually perform the translation
  for (col in names(translations)){
    print(col)
    histo_description[,eval(col):=translations[[col]][factor(get(col))]]
  }
  
  histo_description=histo_description[apply(histo_description,1,function(x){sum(is.na(x))})<8]
  
  #now create actual immunohistochemistry table
  histo_long=melt(histo[,-description_cols,with=FALSE],id.vars = c("N_number","area of HE"),variable.name="measurement",value.name = "counts")
  histo_long[,measurement_group:=unlist(lapply(measurement,function(x){unlist(strsplit(as.character(x),"\\+| "))[1]})),]
  
  #split actual measurement and reference measurement
  histo_long_count=histo_long[!grepl("reference",measurement)]
  histo_long_ref=histo_long[grepl("reference",measurement)]
  
  histo_long_combined=merge(histo_long_count,histo_long_ref,by=c("N_number","area of HE","measurement_group"),suffixes=c("",".ref"),all=TRUE)
  
  histo_long_combined[,counts.norm:=counts/(counts.ref+0.1),]
  histo_long_combined[measurement_group%in%c("HLA-DR","CD34","cell"),counts.norm:=counts/`area of HE`,]
  
  histo_wide=dcast(histo_long_combined,N_number~measurement_group,value.var="counts.norm",fill=NA)  
  return(list(histo_description=histo_description,histo_annotation=as.data.table(histo_wide)))  
}

#prepare clinical annotations (transform patient centered table into sample centered table)
prep_clin_annot=function(input){
  clinical_annot=fread(input)
  
  clinical_annot[clinical_annot==""]=NA
  #introduce patient ID (should be provided by clinicians/ discussed)
  clinical_annot[,patID:=paste0("pat_",formatC(c(1:nrow(clinical_annot)), width = 3, format = "d", flag = "0")),]
  #add 7th surgery column sothat the very last progression without surgery is still included
  clinical_annot[,c("7th surgery N-No","Date of 6th progression","7th surgery","Age at 7th surg"):=NA,]
  
  #replace ambiguous center
  clinical_annot[Center=="1+3",Center:="3",]
  
  #rample centered
  clinical_annot_long1=melt(clinical_annot,measure.vars = c("1st surgery N-No","2nd surgery N-No","3rd surgery N-No","4th surgery N-No","5th surgery N-No","6th surgery N-No","7th surgery N-No"),variable.name="SampleEvent",value.name = "N_number")
  
  mapping=list("1st surgery N-No"=c(NA,"Date of 1st progression","1st surgery","Age at 1st surg",NA,NA,NA),
               "2nd surgery N-No"=c("Date of 1st progression",NA,"2nd surgery","Age at 2nd surg","First-line treatment","RTX","Dose (Gy)"),
               "3rd surgery N-No"=c("Date of 2nd progression",NA,"3rd surgery","Age at 3rd surg","2nd line CTX","2nd line RTX","Dose (Gy).1"),
               "4th surgery N-No"=c("Date of 3rd progression",NA,"4th surgery","Age at 4th surg","3rd line CTX","3rd line RTX","Dose (Gy).2"),
               "5th surgery N-No"=c("Date of 4th progression",NA,"5th surgery","Age at 5th surg","4th line CTX","4th line RTX","Dose (Gy).3"),
               "6th surgery N-No"=c("Date of 5th progression",NA,"6th surgery","Age at 6th surg","5th line CTX","5th line RTX","Dose (Gy).4"),
               "7th surgery N-No"=c("Date of 6th progression",NA,"7th surgery","Age at 7th surg","6th line CTX","6th line RTX","Dose (Gy).5"))
  
  addColumns=c("ProgressionDate","FirstProgressionDate","SurgeryDate","Age","CTXTreatment","RTXTreatment","RTXDose")
  clinical_annot_long1[,eval(addColumns):="NA",]
  
  
  for (i in 1:length(mapping)){
    print(names(mapping[i]))
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][1]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[1]):=as.character(get(mapping[[i]][1])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),ProgressionDate:="NA",]}
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][2]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[2]):=as.character(get(mapping[[i]][2])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),FirstProgressionDate:="NA",]}
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][3]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[3]):=as.character(get(mapping[[i]][3])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),SurgeryDate:="NA",]}
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][4]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[4]):=as.character(get(mapping[[i]][4])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),Age:="NA",]}
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][5]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[5]):=as.character(get(mapping[[i]][5])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),CTXTreatment:="NA",]}
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][6]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[6]):=as.character(get(mapping[[i]][6])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),RTXTreatment:="NA",]}
    if (with(clinical_annot_long1,is(try(get(mapping[[i]][7]),TRUE)))[1]!="try-error"){clinical_annot_long1[SampleEvent==names(mapping[i]),eval(addColumns[7]):=as.character(get(mapping[[i]][7])),]}else{clinical_annot_long1[SampleEvent==names(mapping[i]),RTXDose:="NA",]}
  }
  
  #clean up
  clinical_annot_long1[clinical_annot_long1=="NA"]=NA
  clinical_annot_long1[,unique(na.omit(unlist(mapping))):=NULL,]
  clinical_annot_long1[,grep("Re-resection",names(clinical_annot_long1),value = TRUE):=NULL,]
  clinical_annot_long1=clinical_annot_long1[!(is.na(N_number)&is.na(ProgressionDate)&is.na(CTXTreatment))]
  
  clinical_annot_long1[,surgery:=as.numeric(substr(SampleEvent,1,1)),]
  clinical_annot_long1[,treatment:=as.numeric(substr(SampleEvent,1,1))-1,]
  clinical_annot_long1[,Age:=as.numeric(Age),]
  clinical_annot_long1[is.na(Age),Age:=round((as.Date(ProgressionDate,format="%d.%m.%Y")-as.Date(DateOfBirth,format="%d.%m.%Y"))/365),]
  
  clinical_annot_long1[,ProgressionDate:=as.Date(ProgressionDate,'%d.%m.%Y'),]
  clinical_annot_long1[,FirstProgressionDate:=as.Date(FirstProgressionDate,'%d.%m.%Y'),]
  clinical_annot_long1[,SurgeryDate:=as.Date(SurgeryDate,'%d.%m.%Y'),]
  clinical_annot_long1[,`DateOfDeath_LastFollow-up`:=as.Date(`DateOfDeath_LastFollow-up`,'%d.%m.%Y'),]
  clinical_annot_long1[,DateOfBirth:=as.Date(DateOfBirth,'%d.%m.%Y'),]
  
  clinical_annot_long1[,NoOfProgressions_pat:=sum(!is.na(ProgressionDate)),by=patID]
  clinical_annot_long1[,NoOfProgressions_tum:=unlist(lapply(SurgeryDate,function(x){sum(na.omit(ProgressionDate<=x))})),by=patID]
  clinical_annot_long1[,TreatmentDate:=ProgressionDate,]
  clinical_annot_long1[,TreatmentDate:=ifelse(is.na(TreatmentDate)&SampleEvent!="1st surgery N-No",max(c(TreatmentDate,SurgeryDate),na.rm = TRUE)+90,TreatmentDate),by=patID]
  
  #translations (number to word)
  translations=list(
    Center=c("1"="MedUni Vienna","2"="Linz","3"="Salzburg","4"="Innsbruck","5"="Klagenfurt" ,"6"="Wr.Neustadt", "7"="St.P?lten", "8"="Rudolfstiftung","9"="Graz"),
    Sex=c("1"="m","2"="f"),
    IDH=c("0"="wt","1"="mut"),
    #TumorPhenotype=c("0"="stable","1"="shift sarcoma-like","2"="shift monstrocellular","3"="shift PNET-like"),
    TumorPhenotype=c("0"="other","1"="classic to sarcoma"), 
    Shape_shift_type=c("0"="stable", "1"="classic to sarcoma", "2"="classic to monstrocellular","3"="classic to small celled",  "4"="sarcoma to classic"),
    Shape_shift=c("0"="stable","1"="shape shift"),
    VitalStatus=c("0"="lost",  "1"="dead", "2"="alive"),
    StuppComplete=c("0"="no","1"="yes","2"="not applicable/other treatment","3"="unknown"),
    StuppDiscontinued=c("0"="no","1"="yes","2"="not applicable","3"="unknown"),
    ReasonForDiscontd=c("0"="none or not applicable","1"="tumor progression","2"="hematotoxicity","3"="nausea","4"="unknown","5"="hepatotoxicity","6"="infection"),
    NoOfAdjuvantTMZCyclesCompleted=c("0"="6","1"="concomitant RT/TMZ phase","2"="1","3"="not applicable","4"="3","5"="5","6"="unknown","7"="2","8"="4","9"="during" ),
    AlternativeTreatment=c("0"="none","1"="bevacizumab","2"="re-resection","3"="other CTX","4"="no treatment","5"="imatinib","6"="unknown","7"="TMZ rechallenge"),
    Bevacizumab=c("0"="no","1"="between 1st and 2nd surgery","2"="after 2nd surgery"),
    RTXTreatment=c("0"="no RTX","1"="RTX","2"="Gamma Knife")
  )
  #perform the actual translation
  for (i in 1:length(translations)){
    name=names(translations[i])
    translation=translations[[i]]
    print(name)
    if (!(name%in%names(clinical_annot_long1))){message("missing");next}
    for (j in 1:length(translation)){
      target=names(translation[j])
      value=  translation[j]
      clinical_annot_long1[,eval(name):=as.character(get(name)),]
      clinical_annot_long1[get(name)==target,eval(name):=value,]
    } 
  }
  clinical_annot_long1
  
  #CTX treatment (extra because values have multiple meanings)
  firstCTX=c("1"="Stupp regimen","2"="Stupp+Vaccination", "3"="No therapy or denied","4"="Fotemustine-Dacarbacine","5"="Nordic Glioma/RTX only","6"="PCV", "7"="Stupp+Cilengitide","8"="Unknown","9"="Stupp+INN-Doxorubicine","10"="RTX only","11"="Watchful waiting")
  
  followingCTX=c("0"="none","1"="TMZ","2"="Beva","3"="Beva+Irinotecan", "4"="Imatinib","5"="Lomustine","6"="Sunitinib","7"="TMZ+Imatinib","8"="Beva+Irinotecan+Imatinib","9"="TMZ+Beva", "10"="Thalidomide", "11"="Beva+Lomustine",  "12"="Cediranib","13"="PCV","14"="Novo-TTF","15"="Lomustine+INN-pazopanib","16"="Beva+Imatinib","17"="Fotemustine","18"="Beva+Fotemustine","19"="Beva+Inn-Doxorubicine","20"="Cilengitide","21"="TMZ+Fotemustine","22"="INN-pazopani","23"="Unknown/Best supportive care","24"="Cisplatin+Etoposide","25"="Pazopanib")
  
  clinical_annot_long1[,CTXTreatment:=as.character(CTXTreatment),]
  for (i in 1:length(firstCTX)){
    target=names(firstCTX[i])
    value=  firstCTX[i]
    print(target)
    clinical_annot_long1[treatment==1&CTXTreatment==target,CTXTreatment:=value,]
  }
  
  for (i in 1:length(followingCTX)){
    target=names(followingCTX[i])
    value=  followingCTX[i]
    print(target)
    clinical_annot_long1[treatment!=1&CTXTreatment==target,CTXTreatment:=value,]
  }
  return(clinical_annot_long1)
}
