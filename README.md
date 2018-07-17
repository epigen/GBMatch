# GBMatch
This repository contains code and metadata for the glioblastoma progression study (https://www.biorxiv.org/content/early/2017/08/09/173864). 
The data analysis code in src relies on the processed data (methylation calls, PDR calls, statistics) being available in the directory structure 
produced by a data processing  pipeline based on Pypiper (https://pypiper.readthedocs.io) and Looper (https://looper.readthedocs.io) as specified in the 
[project configuration file](metadata/GBMatch.yaml).
It further relies on the [projectInit R package](https://github.com/databio/projectInit) to read the pipeline results 
into R for further analysis. The scripst found in this repository are meant to be executed in order and produce the plots and tables 
presented in the corresponding publication.

Addititional data can be found on the supplementary website: http://glioblastoma-progression.computational-epigenetics.org

### Linking figures, scrips, and results:

The combined annotation tables used in nearly every analysis are produced using the scripts 01.1-combine_annotation.R (combining the different types of collected annotations) and 10.1_add_annotation.R (adding annotations produced through analysis e.g. transcriptional subtypes) and stored in the folder 01.1-combined_annotation.

|Figure|Scripts|Results|
|------|------|--------------|
Fig. 1a| Overview schematic --> No script|No analysis
Fig. 1b| 01.2-sample_stats.R|01.2-sample_stats
Fig. 1c| 02.4-meth_overview_tracks.R|02-meth_overview
Fig. 1d| 02.2-IDH1_TurcanEtAl.R|02-meth_overview
Fig. 2a| Overview schematic --> No script|No analysis
Fig. 2b| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. 2c| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. 2d| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. 2e| 08.1-GBM_classifier.R + 13.1-survival_analysis.R|13.1-survival
Fig. 2f| 08.1-GBM_classifier.R + 11.3-DiffMeth_groups.R + 11.4-LOLA.R|11.3-DiffMeth_groups
Fig. 2g| 08.1-GBM_classifier.R + 11.3-DiffMeth_groups.R + 11.4-LOLA.R|11.3-DiffMeth_groups
Fig. 2h| Overview schematic --> No script|No analysis
Fig. 2i| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|09-dipScore + 12.1-Annot_Associations_follow-up
Fig. 3a| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. 3b| Histology --> No script|No analysis
Fig. 3c| 13.1-survival_analysis.R|13.1-survival
Fig. 3d| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. 3e| Histology --> No script|No analysis
Fig. 3e| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. 3g| 08.2-meth_pred.R|08.2-meth_pred
Fig. 4a| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. 4b| 13.1-survival_analysis.R|13.1-survival
Fig. 4c| 08.2-meth_pred.R|08.2-meth_pred
Fig. 4d| 08.2-meth_pred.R|08.2-meth_pred
Fig. 4e| 08.2-meth_pred.R|08.2-meth_pred
Fig. 4f| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. 4g| Histology --> No script|No analysis
Fig. 4h| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. 4i| 08.2-meth_pred.R|08.2-meth_pred
Fig. 4j| 13.1-survival_analysis.R|13.1-survival
Fig. 5a| Overview schematic --> No script|No analysis
Fig. 5b| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R|07-meth_heterogeneity
Fig. 5c| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R|07-meth_heterogeneity
Fig. 5d| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R+13.1-survival_analysis.R|13.1-survival
Fig. 5e| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
Fig. 6a| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
Fig. 6b| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
Fig. 6c| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
Fig. 6d| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R + 13.1-survival_analysis.R|13.1-survival
Fig. 6e| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R|11.2-diffMeth_single
Fig. 6f| 11.1-diffMeth.R + 11.2-prom_diffMeth_single.R + 13.1-survival_analysis.R|13.1-survival
Fig. 6g| 08.2-meth_pred.R|08.2-meth_pred
Fig. 6h| 08.2-meth_pred.R|08.2-meth_pred
Fig. S1a| 01.2-sample_stats.R|01.2-sample_stats
Fig. S1b| 01.2-sample_stats.R|01.2-sample_stats
Fig. S1c| 01.2-sample_stats.R|01.2-sample_stats
Fig. S1d| 01.2-sample_stats.R|01.2-sample_stats
Fig. S1e| 02.4-meth_overview_regions.R|02-meth_overview
Fig. S1f| 02.4-meth_overview_regions.R|02-meth_overview
Fig. S1g| 02.4-meth_overview_regions.R|02-meth_overview
Fig. S1h| 02.4-meth_overview_regions.R|02-meth_overview
Fig. S2a| 02.4-meth_overview_tracks.R|02-meth_overview
Fig. S2b| 02.4-meth_overview_tracks.R|02-meth_overview
Fig. S2c| 02.3-MGMT_meth.R|02-meth_overview
Fig. S2d| 02.3-MGMT_meth.R + 13.1-survival_analysis.R|13.1-survival
Fig. S3a| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R|03-CopywriteR
Fig. S3b| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R + 13.1-survival_analysis.R|13.1-survival
Fig. S3c| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R|03-CopywriteR
Fig. S4a| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R + 03.4-validate_CopywrtiteR.R|03-CopywriteR
Fig. S4b| 03.2-run_CopywrtiteR.R + 03.3-plot_CopywrtiteR.R + 03.4-validate_CopywrtiteR.R|03-CopywriteR
Fig. S5a| Overview schematic --> No script|No analysis
Fig. S5b| 08.1-GBM_classifier.R + 08.1-GBM_classifier_validation.R|08.1-GBM_classifier
Fig. S5c| 08.1-GBM_classifier.R + 08.1-GBM_classifier_validation.R|08.1-GBM_classifier
Fig. S5d| 08.1-GBM_classifier.R + 08.1-GBM_classifier_validation.R|08.1-GBM_classifier
Fig. S6a| Overview schematic --> No script|No analysis
Fig. S6b| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. S6c| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. S6d| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. S6e| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. S6f| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. S6g| 08.1-GBM_classifier.R|08.1-GBM_classifier
Fig. S6h| 08.1-GBM_classifier.R + 13.1-survival_analysis.R|13.1-survival
Fig. S7a| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S7b| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S7c| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S7d| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S7e| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|09-dipScore
Fig. S7f| 08.1-GBM_classifier.R + 09.1-dipScore.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S8a| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S8b| Histology --> No script|No analysis
Fig. S8c| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S8d| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S8e| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S8f| 13.1-survival_analysis.R|13.1-survival
Fig. S8g| 08.2-meth_pred.R|08.2-meth_pred
Fig. S9a| MR images --> No script|No analysis
Fig. S9b| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S9c| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S9d| 08.1-GBM_classifier.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S9e| 13.1-survival_analysis.R|13.1-survival
Fig. S9f| MR images --> No script|No analysis
Fig. S9g| 13.1-survival_analysis.R|13.1-survival
Fig. S10a| Overview schematic --> No script|No analysis
Fig. S10b| 08.2-meth_pred.R|08.2-meth_pred
Fig. S10c| 08.2-meth_pred.R|08.2-meth_pred
Fig. S10d| 08.2-meth_pred.R|08.2-meth_pred
Fig. S10e| 08.2-meth_pred.R|08.2-meth_pred
Fig. S10f| 08.2-meth_pred.R|08.2-meth_pred
Fig. S10g| 08.2-meth_pred.R|08.2-meth_pred
Fig. S11a| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S11b| 13.1-survival_analysis.R|13.1-survival
Fig. S11c| 08.2-meth_pred.R|08.2-meth_pred
Fig. S11d| 08.2-meth_pred.R|08.2-meth_pred
Fig. S11e| 08.2-meth_pred.R|08.2-meth_pred
Fig. S11f| 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S12a| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R|07-meth_heterogeneity
Fig. S12b| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R + 12.1-Annotation_Association_folow-up.R|12.1-Annot_Associations_follow-up
Fig. S12c| 05.1-pdr_caches.R + 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 07.1-meth_heterogeneity.R + 13.1-survival_analysis.R|13.1-survival
Fig. S12d| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
Fig. S12e| 06.1-submit_methclone.R + 06.2-analyse_methclone.R + 13.1-survival_analysis.R|13.1-survival
Fig. S12f| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
Fig. S12g| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
Fig. S12h| 06.1-submit_methclone.R + 06.2-analyse_methclone.R|06-methclone
Fig. 13| Screenshots from the data explorer on the supplementary website -->  No script|No analysis
