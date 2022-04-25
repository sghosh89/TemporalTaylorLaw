path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# wrangling data, find good routes, get community matrix minimum for 20 years, 
# species included sampled for atleast 70% of sampling period
source("data_wrangling.R") 