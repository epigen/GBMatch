library(project.init)
library(pheatmap)
project.init2("GBMatch")

# Parameters
preProcessing <- "RAW"
minCoverage <- 5
minCntPerAnnot <- 4

if(!is.na(minCoverage)){
  minCoverage <- as.numeric(minCoverage)
}

# OUT FOLDER --------------------------------------------------------------
out <- NA
dir.create(dirout('11.3-DiffMeth_groups/'))
if(!is.na(minCoverage)){
  out = paste0('11.3-DiffMeth_groups/', preProcessing , "_", minCoverage, "/")
} else {
  out = paste0('11.3-DiffMeth_groups/', preProcessing, "/")
}
dir.create(dirout(out))

# # START  ----------------------------------------------------------------
cacheName <- "rrbsCg"

if(!file.exists(dirout(out, "Matrix.RData"))){
  
  # LOAD DATA ---------------------------------------------------------------
  (load(paste0(getOption("PROCESSED.PROJECT"), "/RCache/", cacheName, ".RData")))
  ret$region <- paste0(ret$chr,"_", ret$start,"_",ret$end)

  if(!is.na(minCoverage)){
    ret <- ret[readCount >= minCoverage]
  }
  
  # CLEAN UP
  value.name <- ifelse("methyl" %in% colnames(ret), "methyl", "PDRa")
  ret <- unique(ret[,c("region", "id", value.name), with=F])
  samples <- unique(ret$id)
  probes <- ret$region
  str(probes)
  x <- table(probes)
  probes.keep <- names(which(x > length(samples) * 0.3))
  ret <- ret[region %in% probes.keep]
  probes <- unique(ret$region)
  
  # TO MATRIX
  ret <- dcast.data.table(ret, region ~ id, value.var=value.name)
  ret.row.names <- ret$region
  ret <- as.matrix(ret[,-"region",with=F])
  row.names(ret) <- ret.row.names
  meth_data_mat <- ret
  rm(list="ret")
  meth_data_mat <- meth_data_mat[sort(rownames(meth_data_mat)),]
  
  save(meth_data_mat, file=dirout(out, "Matrix.RData"))
} else {
  load(dirout(out, "Matrix.RData"))
}

# Prepare annotation for diff analysis
annotation_orig=fread(paste0(dirout(),"01.1-combined_annotation/","annotation_combined_final.tsv"))
annotation_orig <- annotation_orig[category %in% c("GBMatch","GBmatch_val") & IDH == "wt"]
colnames(annotation_orig)
(load(paste0(dirout(),"01.1-combined_annotation/","column_annotation_combined.RData")))
immuno <- column_annotation_combined_clean$histo_immuno
#fcts <- c("category","classic T1",
#  "cT1 flare up",
#  "T2 diffus", 
#  "T2 circumscribed",
#  "Siteofsurgery", 
#  "progression_location", 
#  "WHO2016_classification", 
#  "surgery",
#  "sub_group",
#  "Sex")
# fcts <- c(immuno, fcts)

fcts <- c("category","sub_group","CD163","CD68","EZH2","MIB","Enhancing (mm3)")




# Do differential analysis
res <- data.table(cpg=gsub("_$", "",rownames(meth_data_mat)))

cohorts=list("GBMatch","GBmatch_val",c("GBMatch","GBmatch_val"))

for (sel_cohort in cohorts){ 

  for(factorOfInterest in fcts){
    message(factorOfInterest) 
    tryCatch({
      # ANNOTATION --------------------------------------------------------------
      annotation=annotation_orig[category%in%sel_cohort]
      
      if (sel_cohort==c("GBMatch","GBmatch_val")){
        annotation=annotation_orig[category%in%sel_cohort&surgery==1]
      }
      
      aDat <- annotation[!is.na(get(factorOfInterest)),c(factorOfInterest, "N_number_seq"), with=F]
      if(nrow(aDat) == 0){
        print(paste("Not enough samples to analyze ", factorOfInterest, " moving on"))
        next
      }
      if(nrow(unique(aDat[,factorOfInterest,with=F]))<2){
        print("Only one group. Skipping.")
        next
      }
      
      
      if(factorOfInterest == "sub_group"){
        aDat <- annotation[!is.na(get(factorOfInterest))][sub_group_prob >= 0.8][,c(factorOfInterest, "N_number_seq"), with=F]
      }
      if(factorOfInterest %in% immuno){
        aDat <- annotation[!is.na(get(factorOfInterest))]
        qu <- quantile(aDat[[factorOfInterest]], probs=c(0.2,0.8))
        aDat <- aDat[get(factorOfInterest) <= qu[1] | get(factorOfInterest) >= qu[2]][,c(factorOfInterest, "N_number_seq"), with=F]
        aDat[[factorOfInterest]] <- sapply(aDat[[factorOfInterest]], function(x) ifelse(x <= qu[1], "LOW", "HIGH"))
      }
      colnames(aDat)[1] <- gsub("(\\.|\\_|\\-|\\%|\\s)", "", make.names(gsub("_", "",colnames(aDat)[1])))
      factorOfInterest <- gsub("(\\.|\\_|\\-|\\%|\\s)", "", make.names(gsub("_", "",factorOfInterest)))
      dim(aDat)
      dim(aDat <- aDat[N_number_seq %in% colnames(meth_data_mat)])
      
      # those groups with N > 4
      groups <- aDat[,.N, by=factorOfInterest][N >= minCntPerAnnot][[factorOfInterest]]
      if(length(groups) > 1){
        i1 <- 1
        i2 <- 2
        for(i1 in 1:(length(groups)-1)){
          for(i2 in (i1+1):length(groups)){
            message(groups[i1], " vs ", groups[i2])
            m1 <- meth_data_mat[,colnames(meth_data_mat) %in% aDat[get(factorOfInterest) == groups[i1]]$N_number_seq]
            m2 <- meth_data_mat[,colnames(meth_data_mat) %in% aDat[get(factorOfInterest) == groups[i2]]$N_number_seq]
            medianDiff <- apply(m1, 1, median, na.rm=TRUE) - apply(m2, 1, median, na.rm=TRUE)
            # medianDiff[medianDiff == 0] <- 10
            print(Sys.time())
            pval <- c()
            i <- 1
            for(i in 1:nrow(meth_data_mat)){
              if((sum(!is.na(m1[i,])) >= minCntPerAnnot & sum(!is.na(m2[i,])) >= minCntPerAnnot) & abs(medianDiff[i]) >= 0.1){
                  pval[i] <- tryCatch({
                    wilcox.test(m1[i,], m2[i,])$p.value * sign(medianDiff[i]) 
                    }, error = function(e){
                      print(e)
                      return(1)
                    })
              } else {
                pval[i] <- NA # these have not enough values and aren't tested
              }
            }
            print(Sys.time())
            pval[is.na(pval)] <- NA # get rid of NaN
            res[[paste0(factorOfInterest, "_", groups[i1], "_", groups[i2])]] <- pval
          }
        }
      }
    }, error = function(e) print(e))
  }
  write.table(res, file=dirout(out,paste0("Pvalues_",paste0(sel_cohort,collapse="_"),".tsv")),sep="\t", quote=F, row.names=F)
}
print("DONE")