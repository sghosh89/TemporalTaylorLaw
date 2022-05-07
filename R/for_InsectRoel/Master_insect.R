path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
source("wrangling_data.R") # clean and prepare data
source("insect_ta.R") # tail analysis

# estimate stability and TL slope, portfolio effects etc.
source("stability_and_TL.R")
