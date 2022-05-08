rm(list=ls())
library(tidyverse)
#===================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_swisslake_zoo<-readRDS("../../Results/for_swisslake/summary_table_zooplankton.RDS")
sm_swisslake_zoo$TAXA<-"Freshwater invertebrates"
sm_swisslake_zoo$ORGANISMS<-"Zooplankton"
sm_swisslake_zoo$newsite<-sm_swisslake_zoo$siteid
sm_swisslake_zoo$STUDY_ID<-c("lake hallwil","lake greifensee","lake baldegg")

sm_swisslake_zoo$cv_com<-NA
sm_swisslake_zoo$stability<-NA
sm_swisslake_zoo$vr_LdM<-NA
sm_swisslake_zoo$TLintercept<-NA
sm_swisslake_zoo$TLslope.z<-NA
sm_swisslake_zoo$TLslope.z.lowCI<-NA
sm_swisslake_zoo$TLslope.z.upCI<-NA
sm_swisslake_zoo$pe_avg_cv<-NA
sm_swisslake_zoo$pe_mv<-NA

for(i in 1:nrow(sm_swisslake_zoo)){
  # read input data
  m<-readRDS(paste("../../DATA/for_swisslake/wrangled_data/zooplankton/input_mat_for_tail_analysis_zoo_",sm_swisslake_zoo$siteid[i],".RDS",sep=""))
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_swisslake_zoo$cv_com[i]<-cv_com
  sm_swisslake_zoo$stability[i]<- 1/cv_com
  sm_swisslake_zoo$vr_LdM[i]<-ecofolio::synchrony(m)
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_swisslake_zoo$TLintercept[i]<-res$c
  sm_swisslake_zoo$TLslope.z[i]<-res$z
  sm_swisslake_zoo$TLslope.z.lowCI[i]<-res$z.l
  sm_swisslake_zoo$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_swisslake_zoo$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_swisslake_zoo$pe_mv[i]<-res2
}
sm_swisslake_zoo$source<-"SwissLakeZoo"
sm_swisslake_zoo$REALM<-"Freshwater"
#reorganize
sm_swisslake_zoo<-sm_swisslake_zoo%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                                         initR,nsp,nyr,nind,npos,nL,nU,nneg,L,U,vr_LdM,
                                                         avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                                         stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                                         TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_swisslake_zoo,"../../Results/for_swisslake/stability_and_TL_est_swisslakezoo.RDS")







