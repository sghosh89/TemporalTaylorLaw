rm(list=ls())
library(tidyverse)
#===================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_BBS<-readRDS("../../Results/for_BBS/summary_table_detail_version.RDS")

sm_BBS$cv_com<-NA
sm_BBS$stability<-NA
sm_BBS$vr_LdM<-NA # variance ratio by LdM, synchrony metric in [0,1]
sm_BBS$TLintercept<-NA
sm_BBS$TLslope.z<-NA
sm_BBS$TLslope.z.lowCI<-NA
sm_BBS$TLslope.z.upCI<-NA
sm_BBS$pe_avg_cv<-NA
sm_BBS$pe_mv<-NA

for(i in 1:nrow(sm_BBS)){
  # read input data
  mypath<-paste("../../DATA/for_BBS/wrangled_data/",sm_BBS$siteid[i],"/",sep="")
  m<-readRDS(paste(mypath,"input_mat_for_tailanal.RDS",sep=""))
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_BBS$cv_com[i]<-cv_com
  sm_BBS$stability[i]<- 1/cv_com
  sm_BBS$vr_LdM[i]<-ecofolio::synchrony(m)
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_BBS$TLintercept[i]<-res$c
  sm_BBS$TLslope.z[i]<-res$z
  sm_BBS$TLslope.z.lowCI[i]<-res$z.l
  sm_BBS$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_BBS$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_BBS$pe_mv[i]<-res2
}

sm_BBS$REALM<-"Terrestrial"
sm_BBS$TAXA <-"Birds"
sm_BBS$ORGANISMS <-"Birds"
sm_BBS<-rename(sm_BBS, newsite = siteid) # each siteid is renamed as newsite to be nested within the stratum 
sm_BBS<-rename(sm_BBS, STUDY_ID = Stratum_name) # stratum name renamed as STUDY_ID
sm_BBS$source<-"BBS"

#reorganize
sm_BBS<-sm_BBS%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                 initR,nsp, nind,npos,nL,nU,nneg,L,U,vr_LdM,
                                 avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                 stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                 TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_BBS,"../../Results/for_BBS/stability_and_TL_est_BBS.RDS")







