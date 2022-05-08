path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
source("2.0_wrangling_raw_data.r") # cleaning data
source("3.0_get_tail_analysis_res.r") # tail analysis
source("stability_and_TL.R")