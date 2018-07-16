# GBMatch
This repository contains code and metadata for the glioblastoma progression study. 
The data analysis code in src relies on the processed data (methylation calls, PDR calls, statistics) being available in the directory structure 
produced by a data processing  pipeline based on Pypiper (https://pypiper.readthedocs.io) and Looper (https://looper.readthedocs.io) as specified in the 
[project configuration file](metadata/GBMatch.yaml).
It further relies on the [projectInit R package](https://github.com/databio/projectInit) to read the pipeline results 
into R for further analysis. The scripst found in this repository are meant to be executed in order and produce the plots and tables 
presented in the corresponding publication.
