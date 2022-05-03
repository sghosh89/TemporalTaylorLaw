path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

#------------ plot-level analysis ---------------
source("get_BioTIME_data.R") # reading data
source("freshwater_plotlevel.R") 
source("terrestrial_plotlevel.R") 