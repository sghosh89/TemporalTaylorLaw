rm(list=ls())
library(tidyverse)
#===================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_swisslake_phyto<-readRDS("../../Results/for_swisslake/summary_table_phytoplankton.RDS")
sm_swisslake_phyto$TAXA<-"Freshwater plants" # phytoplanktons are tagged as invertebrates in BioTIME?
sm_swisslake_phyto$ORGANISMS<-"Phytoplankton"
sm_swisslake_phyto$newsite<-sm_swisslake_phyto$siteid
sm_swisslake_phyto$STUDY_ID<-c("lake walensee","lake zurich","lake luzern","lake zurich",
                               "lake sempach","lake hallwil","lake baldegg","lake greifensee")

sm_swisslake_phyto$cv_com<-NA
sm_swisslake_phyto$stability<-NA
sm_swisslake_phyto$TLintercept<-NA
sm_swisslake_phyto$TLslope.z<-NA
sm_swisslake_phyto$TLslope.z.lowCI<-NA
sm_swisslake_phyto$TLslope.z.upCI<-NA
sm_swisslake_phyto$pe_avg_cv<-NA
sm_swisslake_phyto$pe_mv<-NA

for(i in 1:nrow(sm_swisslake_phyto)){
  # read input data
  siteid<-sm_swisslake_phyto$siteid[i]
  m<-readRDS(paste("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_",siteid,".RDS",sep=""))
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_swisslake_phyto$cv_com[i]<-cv_com
  sm_swisslake_phyto$stability[i]<- 1/cv_com
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_swisslake_phyto$TLintercept[i]<-res$c
  sm_swisslake_phyto$TLslope.z[i]<-res$z
  sm_swisslake_phyto$TLslope.z.lowCI[i]<-res$z.l
  sm_swisslake_phyto$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_swisslake_phyto$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_swisslake_phyto$pe_mv[i]<-res2
}
sm_swisslake_phyto$source<-"SwissLakePhyto"
sm_swisslake_phyto$REALM<-"Freshwater"
#reorganize
sm_swisslake_phyto<-sm_swisslake_phyto%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                 initR,nsp, nind,npos,nL,nU,nneg,L,U,
                                 avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                 stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                 TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_swisslake_phyto,"../../Results/for_swisslake/stability_and_TL_est_swisslakephyto.RDS")







