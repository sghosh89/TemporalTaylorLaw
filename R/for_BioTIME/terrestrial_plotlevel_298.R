rm(list=ls())
source("tail_analysis.R")
source("monthly_rarefy.R")
library(tidyverse)
`%notin%` <- Negate(`%in%`)
xxm<-readRDS("../../DATA/for_BioTIME/BioTIME_public_private_metadata.RDS")
grid_terres<-readRDS("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/bt_terres_min20yr_rawdata.RDS")
df<-readRDS("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/table_for_map.RDS")
df<-df%>%filter(site==298)
# single lat-lon is reported, but multiple plots

#----------- create result folder for wrangle data -------------------------
resloc<-"../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/298/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}
#--------------------------------------------------------------------------------

site<-df$site
x<-grid_terres%>%filter(STUDY_ID==site)
x<-x%>%mutate(newsite=paste("STUDY_ID_",site,"_PLOT_",PLOT,sep=""))
newsite<-sort(unique(x$newsite))

# check if each newsite visited for >20 years?
tt<-x%>%group_by(newsite)%>%summarize(n=n_distinct(YEAR))%>%ungroup()

# include sites which are sampled > 20 years
tt<-tt%>%filter(n>=20)

