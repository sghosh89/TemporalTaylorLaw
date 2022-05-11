rm(list=ls())
library(tidyverse)
`%notin%` <- Negate(`%in%`)
bt<-readRDS("../../DATA/for_BioTIME/BioTIME_public_private_query_data_metadata.RDS")
data_pt_thrs<-20 # 20 years minimum

grid <- bt %>% 
  dplyr::select(STUDY_ID, PLOT, DAY, MONTH, YEAR, 
                GENUS_SPECIES, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS,
                CLIMATE, REALM, TAXA, ABUNDANCE_TYPE, BIOMASS_TYPE, 
                LATITUDE, LONGITUDE, CENT_LAT, CENT_LONG, NUMBER_LAT_LONG, SUMMARY_METHODS)
colnames(grid)[6:8] <- c('Species', 'Abundance', 'Biomass')

grid<-grid%>%group_by(STUDY_ID)%>%filter(n_distinct(YEAR)>=data_pt_thrs) %>% ungroup()
grid_terres<-grid%>%filter(REALM=="Terrestrial")

#===================== generate results folder for terrestrial ===============
resloc<-"../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}
saveRDS(grid_terres,paste(resloc,"bt_terres_min20yr_rawdata.RDS",sep=""))
grid_terres<-readRDS("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/bt_terres_min20yr_rawdata.RDS")
#============================================================================
# now watch each STUDY_ID
site<-sort(unique(grid_terres$STUDY_ID))
df<-data.frame(site=site)
df$nyr<-NA
df$nPLOT<-NA
df$NUMBER_LAT_LONG<-NA
df$LATmin<-NA
df$LATmax<-NA
df$LONmin<-NA
df$LONmax<-NA
df$monthlyfreqsamp<-NA
df$n_methods<-NA

library(htmltools) 
library(htmlwidgets)
library(leaflet) 

for(i in 1:nrow(df)){
  dat<-grid_terres%>%filter(STUDY_ID==site[i])
  df$nyr[i]<-length(unique(dat$YEAR))
  df$nPLOT[i]<-list(unique(dat$PLOT))
  df$NUMBER_LAT_LONG[i]<-unique(dat$NUMBER_LAT_LONG)
  df$LATmin[i]<-min(unique(dat$LATITUDE))
  df$LATmax[i]<-max(unique(dat$LATITUDE))
  df$LONmin[i]<-min(unique(dat$LONGITUDE))
  df$LONmax[i]<-max(unique(dat$LONGITUDE))
  t1<-dat%>%group_by(YEAR)%>%summarise(n_distinct(MONTH))%>%ungroup()
  df$monthlyfreqsamp[i]<-list(range(t1$`n_distinct(MONTH)`))
  df$n_methods[i]<-list(unique(dat$SUMMARY_METHODS))
  
  #---------- save sampling sites on map ----------
  dat<-dat%>%dplyr::select(STUDY_ID,LATITUDE,LONGITUDE)%>%distinct()
  
  sitemap<-leaflet(dat) %>% addTiles() %>%
    addMarkers(~LONGITUDE, ~LATITUDE, label = ~htmlEscape(STUDY_ID))
  f<-paste("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/samplingsite_",
           dat$STUDY_ID[1],".html",sep="")
  htmlwidgets::saveWidget(sitemap, 
                          file.path(normalizePath(dirname(f)),basename(f)))
}
saveRDS(df,"../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/table_for_map.RDS")

#-------- create res folder ------------
resloc2<-"../../Results/for_BioTIME/Terrestrial_plotlevel/"
if(!dir.exists(resloc2)){
  dir.create(resloc2)
}

#-----------------------------
source("terrestrial_plotlevel_18.R")
source("terrestrial_plotlevel_39.R")
source("terrestrial_plotlevel_42.R")
source("terrestrial_plotlevel_46.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_47.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_54.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_56.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_59.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_63.R") # no plots have>=15 spdragonfly: lake ecosystem, still terrestrial?
source("terrestrial_plotlevel_67.R") # no plots have >=15 sp.
#source("terrestrial_plotlevel_195.R") # STUDY_ID=195 is BBS data - so we excluded here
source("terrestrial_plotlevel_214.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_215.R")
# source("terrestrial_plotlevel_221.R") # not a single lat-lon sampled atleast for 20 yrs, if want to include this study aggregate all
source("terrestrial_plotlevel_243.R")# no plots have >=15 sp.
#source("terrestrial_plotlevel_298.R") # not a single lat-lon sampled atleast for 20 yrs, if want to include this study aggregate all
# source("terrestrial_plotlevel_300.R") # STUDY_ID=300 is BioTIMEx data for landis_2018 - so we excluded here
source("terrestrial_plotlevel_301.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_308.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_311.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_333.R")
source("terrestrial_plotlevel_339.R")
source("terrestrial_plotlevel_355.R")
#source("terrestrial_plotlevel_356.R")# not a single lat-lon sampled atleast for 20 yrs, if want to include this study aggregate all
#source("terrestrial_plotlevel_360.R")# all raresp, warnings!
source("terrestrial_plotlevel_361.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_363.R")
source("terrestrial_plotlevel_366.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_413.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_414.R") #
source("terrestrial_plotlevel_416.R") # no plots have >=15 sp.
source("terrestrial_plotlevel_420.R")# no plots have >=15 sp.
#source("terrestrial_plotlevel_483.R") # all raresp, warnings! # no plots have >=15 sp.
#source("terrestrial_plotlevel_497.R") # all raresp, warnings! # no plots have >=15 sp.
source("terrestrial_plotlevel_528.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_529.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_530.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_531.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_532.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_533.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_534.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_536.R")# no plots have >=15 sp.
source("terrestrial_plotlevel_538.R")
#source("terrestrial_plotlevel_540.R") # # no site left with min 20 year sampling


df<-readRDS("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/table_for_map.RDS")
df_included<-df%>%filter(site%notin%c(46,47,54,56,59,63,67,
                                      195,214,243,221,298,
                                      300,301,308,311,
                                      356,360,361,366,
                                      413,416,420,483,497,
                                      528,529,530,531,532,533,534,536,540))
saveRDS(df_included,"../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/table_for_map_selected.RDS")

#--------------- Do a summary stats for terrestrial sites ------------------
source("compute_avg_cor.R")
summary_table<-c()
for (i in c(1:length(df_included$site))){
  resloc<-paste("../../DATA/for_BioTIME/wrangled_data/Terrestrial_plotlevel/",df_included$site[i],"/",sep="")
  newsitelist<-readRDS(paste(resloc,"newsite.RDS",sep=""))
  
  if(length(newsitelist)==1){
    tempo<-newsitelist==df_included$site[i]
    if(tempo==T){
      readpath<-resloc
      resloc2<-paste("../../Results/for_BioTIME/Terrestrial_plotlevel/",df_included$site[i],"/",df_included$site[i],"/",sep="")
    }else{
      readpath<-paste(resloc,newsitelist,"/",sep="")
      resloc2<-paste("../../Results/for_BioTIME/Terrestrial_plotlevel/",df_included$site[i],"/",newsitelist,"/",sep="")
    }
    st<-readRDS(paste(resloc2,"summary_df.RDS",sep=""))
    st$STUDY_ID<-df_included$site[i]
    st$newsite<-newsitelist
    bigM<-readRDS(paste(readpath,"spmat.RDS",sep=""))
    st$initR<-ncol(bigM$spmat)
    mat<-bigM$spmat
    id<-which(colnames(mat)=="raresp")
    if(length(id)>0){
      mat<-mat[,-id]
    }
    res<-compute_avg_cor(mat=mat)
    saveRDS(res,paste(resloc2,"avg_cor_metric.RDS",sep=""))
    st<-cbind(st,res)
    
    summary_table<-rbind(summary_table,st)
  }else{
    for(j in 1:length(newsitelist)){
      readpath<-paste(resloc,newsitelist[j],"/",sep="")
      resloc2<-paste("../../Results/for_BioTIME/Terrestrial_plotlevel/",df_included$site[i],"/",newsitelist[j],"/",sep="")
      st<-readRDS(paste(resloc2,"summary_df.RDS",sep=""))
      st$STUDY_ID<-df_included$site[i]
      st$newsite<-newsitelist[j]
      bigM<-readRDS(paste(readpath,"spmat.RDS",sep=""))
      st$initR<-ncol(bigM$spmat)
      
      mat<-bigM$spmat
      id<-which(colnames(mat)=="raresp")
      if(length(id)>0){
        mat<-mat[,-id]
      }
      res<-compute_avg_cor(mat=mat)
      saveRDS(res,paste(resloc2,"avg_cor_metric.RDS",sep=""))
      st<-cbind(st,res)
      
      summary_table<-rbind(summary_table,st)
    }
  }
}
# reorganize
summary_table<-summary_table%>%dplyr::select(STUDY_ID, newsite,initR,nsp,nyr,nint,nind,npos,nL,nU,nneg,L,U,
                                             avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp)

saveRDS(summary_table,"../../Results/for_BioTIME/Terrestrial_plotlevel/summary_table.RDS")



