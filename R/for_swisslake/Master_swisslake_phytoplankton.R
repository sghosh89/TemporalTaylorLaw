path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
source("lakedata_cleaning_phytoplankton.R") # clean the raw data (>=20 yrs)
source("select_sp_forBlake.R") # genus aggregation
source("get_input_spmat_phytoplankton.R") # get input matrix for tail analysis (>=15 sp.)
source("tail_analysis_phytoplankton.R") # tail analysis results