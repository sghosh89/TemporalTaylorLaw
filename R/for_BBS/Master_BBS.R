path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
source("data_wrangling.R") # wrangling data and find good routes
