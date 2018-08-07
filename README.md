# GBMatch
This repository contains code and metadata for the glioblastoma progression study (https://www.biorxiv.org/content/early/2017/08/09/173864). 
The data analysis code in [src](src) relies on the processed data (methylation calls, PDR calls, statistics) being available in the directory structure produced by a data processing pipeline based on Pypiper (https://pypiper.readthedocs.io) and Looper (https://looper.readthedocs.io) as specified in the [project configuration file](metadata/GBMatch.yaml).
It further relies on the [projectInit R package](https://github.com/databio/projectInit) to read the pipeline results 
into R for further analysis. The scripts found in this repository are meant to be executed in order and produce the plots and tables 
presented in the corresponding publication.

Addititional data can be found on the supplementary website: http://glioblastoma-progression.computational-epigenetics.org

### Linking figures, scrips, and results:

Additionally to the links between figures and scripts listed below, the positions in the scripts where the respective plots are produced are marked with *####Figure X*. 

The combined annotation tables used in nearly every analysis are produced using the scripts *01.1-combine_annotation.R* (combining the different types of collected annotations) and *10.1_add_annotation.R* (adding annotations produced through analysis e.g. transcriptional subtypes) and stored in the folder *01.1-combined_annotation*.

Figure|Scripts|Results
------|------|--------------
1a| Overview schematic --> No script|No analysis||
1b| 01.2-sample_stats.R|01.2-sample_stats
1c| 02.4-meth_overview_tracks.R|02-meth_overview
1d| 02.2-IDH1_TurcanEtAl.R|02-meth_overview
2a| Overview schematic --> No script|No analysis
2b| 08.1-GBM_classifier.R|08.1-GBM_classifier
2c| 08.1-GBM_classifier.R|08.1-GBM_classifier
2d| 08.1-GBM_classifier.R|08.1-GBM_classifier
2e| 08.1-GBM_classifier.R + 13.1-survival_analysis.R|13.1-survival
2f| 08.1-GBM_classifier.R + 11.3-DiffMeth_groups.R + 11.4-LOLA.R|11.3-DiffMeth_groups
2g| 08.1-GBM_classifier.R + 11.3-DiffMeth_groups.R + 11.4-LOLA.R|11.3-DiffMeth_groups
2h| Overview schematic --> No script|No analysis
2i| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|09-dipScore + 12.1-Annot_Associations_follow-up
3a| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
3b| Histology --> No script|No analysis
3c| 13.1-survival_analysis.R|13.1-survival
3d| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
3e| Histology --> No script|No analysis
3e| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
3g| 08.2-meth_pred.R|08.2-meth_pred
4a| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
4b| 13.1-survival_analysis.R|13.1-survival
4c| 08.2-meth_pred.R|08.2-meth_pred
4d| 08.2-meth_pred.R|08.2-meth_pred
4e| 08.2-meth_pred.R|08.2-meth_pred
4f| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
4g| Histology --> No script|No analysis
4h| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
4i| 08.2-meth_pred.R|08.2-meth_pred
4j| 13.1-survival_analysis.R|13.1-survival
5a| Overview schematic --> No script|No analysis
5b| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R|07-meth_heterogeneity
5c| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R|07-meth_heterogeneity
5d| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R+13.1-survival_analysis.R|13.1-survival
5e| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
6a| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
6b| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
6c| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
6d| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R + 13.1-survival_analysis.R|13.1-survival
6e| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
6f| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R + 13.1-survival_analysis.R|13.1-survival
6g| 08.2-meth_pred.R|08.2-meth_pred
6h| 08.2-meth_pred.R|08.2-meth_pred
S1a| 01.2-sample_stats.R|01.2-sample_stats
S1b| 01.2-sample_stats.R|01.2-sample_stats
S1c| 01.2-sample_stats.R|01.2-sample_stats
S1d| 01.2-sample_stats.R|01.2-sample_stats
S1e| 02.4-meth_overview_regions.R|02-meth_overview
S1f| 02.4-meth_overview_regions.R|02-meth_overview
S1g| 02.4-meth_overview_regions.R|02-meth_overview
S1h| 02.4-meth_overview_regions.R|02-meth_overview
S2a| 02.4-meth_overview_tracks.R|02-meth_overview
S2b| 02.4-meth_overview_tracks.R|02-meth_overview
S2c| 02.3-MGMT_meth.R|02-meth_overview
S2d| 02.3-MGMT_meth.R + 13.1-survival_analysis.R|13.1-survival
S3a| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R|03-CopywriteR
S3b| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R + 13.1-survival_analysis.R|13.1-survival
S3c| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R|03-CopywriteR
S4a| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R + 03.4-validate_CopywrtiteR.R|03-CopywriteR
S4b| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R + 03.4-validate_CopywrtiteR.R|03-CopywriteR
S5a| Overview schematic --> No script|No analysis
S5b| 08.1-GBM_classifier.R + 08.1-GBM_classifier_validation.R|08.1-GBM_classifier
S5c| 08.1-GBM_classifier.R + 08.1-GBM_classifier_validation.R|08.1-GBM_classifier
S5d| 08.1-GBM_classifier.R + 08.1-GBM_classifier_validation.R|08.1-GBM_classifier
S6a| Overview schematic --> No script|No analysis
S6b| 08.1-GBM_classifier.R|08.1-GBM_classifier
S6c| 08.1-GBM_classifier.R|08.1-GBM_classifier
S6d| 08.1-GBM_classifier.R|08.1-GBM_classifier
S6e| 08.1-GBM_classifier.R|08.1-GBM_classifier
S6f| 08.1-GBM_classifier.R|08.1-GBM_classifier
S6g| 08.1-GBM_classifier.R|08.1-GBM_classifier
S6h| 08.1-GBM_classifier.R + 13.1-survival_analysis.R|13.1-survival
S7a| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S7b| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S7c| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S7d| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S7e| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|09-dipScore
S7f| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S8a| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S8b| Histology --> No script|No analysis
S8c| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S8d| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S8e| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S8f| 13.1-survival_analysis.R|13.1-survival
S8g| 08.2-meth_pred.R|08.2-meth_pred
S9a| MR images --> No script|No analysis
S9b| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S9c| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S9d| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S9e| 13.1-survival_analysis.R|13.1-survival
S9f| MR images --> No script|No analysis
S9g| 13.1-survival_analysis.R|13.1-survival
S10a| Overview schematic --> No script|No analysis
S10b| 08.2-meth_pred.R|08.2-meth_pred
S10c| 08.2-meth_pred.R|08.2-meth_pred
S10d| 08.2-meth_pred.R|08.2-meth_pred
S10e| 08.2-meth_pred.R|08.2-meth_pred
S10f| 08.2-meth_pred.R|08.2-meth_pred
S10g| 08.2-meth_pred.R|08.2-meth_pred
S11a| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S11b| 13.1-survival_analysis.R|13.1-survival
S11c| 08.2-meth_pred.R|08.2-meth_pred
S11d| 08.2-meth_pred.R|08.2-meth_pred
S11e| 08.2-meth_pred.R|08.2-meth_pred
S11f| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S12a| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R|07-meth_heterogeneity
S12b| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
S12c| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R + 13.1-survival_analysis.R|13.1-survival
S12d| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
S12e| 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 13.1-survival_analysis.R|13.1-survival
S12f| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
S12g| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
S12h| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
13| Screenshots from the data explorer (http://glioblastoma-progression.computational-epigenetics.org) -->  No script|No analysis
