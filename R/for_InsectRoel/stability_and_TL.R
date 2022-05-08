rm(list=ls())
library(tidyverse)
#===================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_insect<-readRDS("../../Results/for_insectRoel/summary_table_detail_version.RDS")

sm_insect$cv_com<-NA
sm_insect$stability<-NA
sm_insect$vr_LdM<-NA
sm_insect$TLintercept<-NA
sm_insect$TLslope.z<-NA
sm_insect$TLslope.z.lowCI<-NA
sm_insect$TLslope.z.upCI<-NA
sm_insect$pe_avg_cv<-NA
sm_insect$pe_mv<-NA

for(i in 1:nrow(sm_insect)){
  # read input data
  mypath<-paste("../../DATA/for_insectRoel/wrangled_data/",sm_insect$STUDY_ID[i],"/",sm_insect$newsite[i],"/",sep="")
  m<-readRDS(paste(mypath,"inputmat_for_tailanal.RDS",sep=""))
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_insect$cv_com[i]<-cv_com
  sm_insect$stability[i]<- 1/cv_com
  sm_insect$vr_LdM[i]<-ecofolio::synchrony(m)
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_insect$TLintercept[i]<-res$c
  sm_insect$TLslope.z[i]<-res$z
  sm_insect$TLslope.z.lowCI[i]<-res$z.l
  sm_insect$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_insect$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_insect$pe_mv[i]<-res2
}
sm_insect$source<-"InsectRoel"

#reorganize
sm_insect<-sm_insect%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                 initR,nsp,nyr,nind,npos,nL,nU,nneg,L,U,vr_LdM,
                                 avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                 stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                 TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_insect,"../../Results/for_insectRoel/stability_and_TL_est_insect.RDS")







