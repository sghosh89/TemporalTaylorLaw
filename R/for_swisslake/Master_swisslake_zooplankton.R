path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
source("lakedata_cleaning_for_ZP.R")# clean the raw data
source("tail_analysis_zooplankton.R")# tail analysis results
