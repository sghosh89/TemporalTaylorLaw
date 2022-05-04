rm(list=ls())
source("tail_analysis.R")
source("monthly_rarefy.R")
library(tidyverse)
`%notin%` <- Negate(`%in%`)
xxm<-readRDS("../../DATA/for_BioTIME/BioTIME_public_private_metadata.RDS")
grid_terres<-readRDS("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/bt_terres_min20yr_rawdata.RDS")
df<-readRDS("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/table_for_map.RDS")
df<-df%>%filter(site==221)

# multiple lat-lon is reported, sampled for different monthly freq in year

#----------- create result folder for wrangle data -------------------------
resloc<-"../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/221/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}
#--------------------------------------------------------------------------------

site<-df$site
x<-grid_terres%>%filter(STUDY_ID==site)
x<-x%>%mutate(newsite=paste("STUDY_ID_",site,"_LAT_",LATITUDE,"_LON_",LONGITUDE,sep=""))
newsite<-sort(unique(x$newsite))

# check if each newsite visited for >20 years?
tt<-x%>%group_by(newsite)%>%summarize(n=n_distinct(YEAR))%>%ungroup()

# include sites which are sampled > 20 years
tt<-tt%>%filter(n>=20)

#update
x_allsite<- x %>% filter(newsite %in% tt$newsite)
newsite<-tt$newsite

t2<-x_allsite%>%group_by(newsite,YEAR)%>%summarize(n=n_distinct(MONTH))%>%ungroup()

