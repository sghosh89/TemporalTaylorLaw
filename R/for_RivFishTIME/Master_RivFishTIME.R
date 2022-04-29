path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# wrangling data, find good routes, get community matrix minimum for 20 years, 
# species included sampled for atleast 70% of sampling period
source("data_wrangling.R") 

# computing tail-dep. between rows for a community matrix (year by species) and summarized
source("RivFishTIME_ta.R") 
