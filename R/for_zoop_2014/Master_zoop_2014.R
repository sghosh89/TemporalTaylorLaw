path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
source("cleaning_zoop_2014.R")# clean the raw data
# we could not find a single lake with 20 years sampling where 15 species 
# were present at least, so we will exclude this zoop2014 data from the analysis
