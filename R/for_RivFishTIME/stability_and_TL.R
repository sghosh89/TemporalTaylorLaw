rm(list=ls())
library(tidyverse)
#===================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_RF<-readRDS("../../Results/for_RivFishTIME/summary_table_detail_version.RDS")

sm_RF$cv_com<-NA
sm_RF$stability<-NA
sm_RF$vr_LdM<-NA
sm_RF$TLintercept<-NA
sm_RF$TLslope.z<-NA
sm_RF$TLslope.z.lowCI<-NA
sm_RF$TLslope.z.upCI<-NA
sm_RF$pe_avg_cv<-NA
sm_RF$pe_mv<-NA

for(i in 1:nrow(sm_RF)){
  # read input data
  m<-readRDS(paste("../../DATA/for_RivFishTIME/wrangled_data/",sm_RF$siteid[i],"/commonspecies_timeseries.RDS",sep=""))
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_RF$cv_com[i]<-cv_com
  sm_RF$stability[i]<- 1/cv_com
  sm_RF$vr_LdM[i]<-ecofolio::synchrony(m)
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_RF$TLintercept[i]<-res$c
  sm_RF$TLslope.z[i]<-res$z
  sm_RF$TLslope.z.lowCI[i]<-res$z.l
  sm_RF$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_RF$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_RF$pe_mv[i]<-res2
}

sm_RF$REALM<-"Freshwater"
sm_RF$TAXA <-"Fish"
sm_RF$ORGANISMS <-"Fish"
sm_RF<-rename(sm_RF, newsite = siteid) # each siteid is renamed as newsite to be nested within the hydrobasin 
sm_RF<-rename(sm_RF, STUDY_ID = HydroBasin) # Hydrobasin renamed as STUDY_ID
sm_RF$source<-"RivFishTIME"

#reorganize
sm_RF<-sm_RF%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                 initR,nsp, nind,npos,nL,nU,nneg,L,U,vr_LdM,
                                 avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                 stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                 TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_RF,"../../Results/for_RivFishTIME/stability_and_TL_est_RF.RDS")







