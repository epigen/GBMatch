#Original author: N. Fortelny, adapted by J. Klughammer
library(project.init)
project.init2("GBMatch")
library(pheatmap)

out="02-meth_overview/"
dir.create(dirout(out))

annot=fread(file.path(getOption("PROCESSED.PROJECT"),"results_analysis/01.1-combined_annotation/annotation_combined.tsv"))


#### GET PROBES FROM TURCAN ET AL
table1 <- data.table(read.csv2("metadata/Turcan_etAl_2012/Table1.csv", skip = 1))
table6 <- data.table(read.csv2("metadata/Turcan_etAl_2012/Table6.csv", skip = 1))
colnames(table1)[3:4] <- c("Foldchange", "Delta")
colnames(table6)[3:4] <- c("Foldchange", "Pvalue")
quantile(table6$Foldchange)

# Probes with highest change mut / wt 
turcan <- list(
	table6_maxFC = table6[order(-abs(Foldchange))][1:3000],
	table6_pnFC = rbind(
		table6[Foldchange < 6][order(-Foldchange)][1:2000], # positive but matching about the negative foldchange
		table6[order(Foldchange)][1:2000] # negative
	)
)
sapply(turcan, function(x) quantile(x$Foldchange))

# map probes to chromosomes
source("src/99-illumina450k_liftover.R") 
probe.map <- Illumina450k_to_hg38_forTurcanEtAl

turcan <- lapply(turcan, function(table)
	cbind(
		table,
		probe.map[match(table$Probeset.ID, probe.map$ID),c("chr", "pos","strand", "ID"), with = F]
		))

# ID mapping worked
stopifnot(all(sapply(turcan, function(table) all(table$Probeset.ID == table$ID))))
# type turcan to see that at least chromosome position was transferable

#### GET OUR RRBS VALUES FOR THE TURCAN PROBES

# make a column cPos as "chromosome position"
turcan <- lapply(turcan, function(x){ x$cPos <- paste0(x$chr,"_", x$pos);return(x) })
# Just to control, how many do I get if i match +/- 1
turcan <- lapply(turcan, function(x){ x$cPos_m1 <- paste0(x$chr,"_", x$pos-1);return(x) })
turcan <- lapply(turcan, function(x){ x$cPos_p1 <- paste0(x$chr,"_", x$pos+1);return(x) })


# Extract data with the "original positions" --> FULL DATA 
loadCaches("rrbsCg")
rrbsData=rrbsCg
rm(rrbsCg)
if(!"regionID"%in%colnames(rrbsData)){rrbsData[,regionID:=paste0(chr,"_",start),]}	
sapply(turcan, function(x) rrbsData[regionID %in% x$cPos, .N ])
# table1 table6
#  36058  69654
sapply(turcan, function(x) rrbsData[regionID %in% x$cPos_m1, .N])
sapply(turcan, function(x) rrbsData[regionID %in% x$cPos_p1, .N])
rrbsData.turcan <- lapply(turcan, function(x) rrbsData[regionID %in% x$cPos])

# LIMIT TO THOSE CpGs with not too many NAs
sapply(rrbsData.turcan, function(x) quantile(x[,.N,by=regionID]$N))
sapply(rrbsData.turcan, function(x) length(unique(x$regionID)))
# keep only those which have > 100 values per CpG
rrbsData.turcan2 <- lapply(rrbsData.turcan, function(x) x[,regionCnt:=.N, by=regionID][regionCnt >= 60])
sapply(rrbsData.turcan2, function(x) quantile(x[,.N,by=id]$N))
rrbsData.turcan2 <- lapply(rrbsData.turcan2, function(x) x[,idCnt:=.N, by=id][idCnt >= 20])
sapply(rrbsData.turcan, function(x) quantile(x[,.N,by=regionID]$N))
sapply(rrbsData.turcan2, function(x) quantile(x[,.N,by=regionID]$N))


# make wide matrices
rrbsData.turcan2.wide <- lapply(rrbsData.turcan2, function(x)
	reshape(x[, c("id", "regionID", "methyl"), with=F],timevar="id", idvar="regionID", direction="wide"))
wideDat.turcan <- lapply(rrbsData.turcan2.wide, function(x1)	{
	# rownames(x1) <- x1$regionID
	# x1 <- x1[,-"regionID", with=F]
	colnames(x1) <- gsub("methyl\\.", "", colnames(x1))
	return(x1)})

# only two probes in the maxFC table are FC < 0, thus not visible in the heatmap
table(turcan$table6_maxFC[cPos %in% wideDat.turcan$table6_maxFC$regionID]$Foldchange > 0)
# So below I will only do it for the ncFC to see up and down probes

# For each matrix get the annotations and plot the heatmap
for(x.nam in c("table6_pnFC")){ 
	# x1 contains the region ID column (first column), xPlot does not
	x1 <- wideDat.turcan[[x.nam]]
	xPlot <- x1[,-"regionID", with=F]
	xPlot <- as.data.frame(xPlot)
	rownames(xPlot) <- x1$regionID
	annot_red=annot[N_number_seq %in% colnames(x1)]
	annot_red=annot_red[match(colnames(x1),N_number_seq),]

	# Column annotation
	annot_col=with(annot_red[!is.na(N_number_seq)&!duplicated(N_number_seq)&!category%in%c("GLASS","control","GBmatch_rcl")], data.frame(Cycles=round(enrichmentCycles),IDH=IDH,sex=sex,row.names=N_number_seq))
  
	# Row annotation
	annot_row=data.frame(foldchange= turcan[[x.nam]][match(x1$regionID, turcan[[x.nam]]$cPos)]$Foldchange)
	rownames(annot_row) = rownames(xPlot)
	
	# COLORS
	colors=list(
		IDH=c(mut="#ff9289",wt="#00d8e0"),
		Cycles=colorRampPalette(c("green", "red"))( 11 ),
		sex=c(m="#00d8e0",f="#ff9289"),
		foldchange=c("blue", "red")
	)

  #only use samples selected in annot_col
  xPlot=xPlot[,row.names(annot_col)]
  
  
	# CLUSTERING with NAs
	colDist <- dist(t(xPlot))
	colDist[is.na(colDist)] <- max(colDist, na.rm=T)
	stopifnot(sum(is.na(colDist)) == 0)
	rowDist <- dist(xPlot[,rownames(subset(annot_col, IDH == "mut"))])
	rowDist[is.na(rowDist)] <- max(rowDist, na.rm=T)
	stopifnot(sum(is.na(rowDist)) == 0)
	
	# Display rows based on FC by Turcan
	xPlot_ordered <- xPlot[turcan[[x.nam]][cPos %in% rownames(xPlot)][order(-Foldchange)]$cPos,row.names(annot_col[order(annot_col$IDH,annot_col$Cycles),])]
  xPlot <- xPlot[turcan[[x.nam]][cPos %in% rownames(xPlot)][order(-Foldchange)]$cPos,]	

	# HEATMAP
	pdf(dirout(out,"IDH_heatmap_", x.nam,".pdf"),width=30,height=22)
	pheatmap(xPlot,
		show_rownames=F,
		annotation_col=annot_col,
		annotation_row=annot_row,
		annotation_colors=colors,
		color=colorRampPalette(c("blue" ,"red"))(20),
		kmeans_k=NA,
		clustering_distance_cols=colDist,
		cluster_rows=F)
  pheatmap(xPlot_ordered,
   show_rownames=F,
   annotation_col=annot_col,
   annotation_row=annot_row,
   annotation_colors=colors,
   color=colorRampPalette(c("blue" ,"red"))(20),
   kmeans_k=NA,
   clustering_distance_cols=colDist,
   cluster_rows=F,cluster_cols=F)

	dev.off()
}

