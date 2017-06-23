# GBMatch
This repository contains code and metadata for the glioblastoma progression study. 
The data analysis code in src relies on the processed data (methylation calls, PDR calls, statistics) being available in the directory structure 
produced by the analysis pipeline based on Pypiper (http://databio.org/pypiper) and Looper (http://databio.org/looper) as specified in the 
[project configuration file](metadata/GBMatch.yaml).
It further relies on the [project.init R package](https://github.com/nsheff/project.init) to read the pipeline results 
into R for further analysis. The scripst found in this repository are meant to be executed in order and produce the plots and tables 
presented in the corresponding publication.
