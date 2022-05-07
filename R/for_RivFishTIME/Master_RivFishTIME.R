path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# wrangling data, find good routes, get community matrix minimum for 20 years, 
# minimum 15 species included sampled for atleast 70% of sampling period
source("wrangling_data.R") 

# computing tail-dep. between rows for a community matrix (year by species) and summarized
source("RivFishTIME_ta.R") 

# estimate stability and TL slope, portfolio effects etc.
source("stability_and_TL.R")

