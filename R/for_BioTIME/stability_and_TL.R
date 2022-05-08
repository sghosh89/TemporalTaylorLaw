rm(list=ls())
library(tidyverse)
#===================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_frs<-readRDS("../../Results/for_BioTIME/Freshwater_plotlevel/summary_table.RDS")
sm_ter<-readRDS("../../Results/for_BioTIME/Terrestrial_plotlevel/summary_table.RDS")
sm_BT<-rbind(sm_frs,sm_ter)
xxm<-readRDS("../../DATA/for_BioTIME/BioTIME_public_private_metadata.RDS")
xxm<-xxm%>%dplyr::select(STUDY_ID,TAXA,REALM,ORGANISMS)
sm_BT<-inner_join(sm_BT,xxm,"STUDY_ID")

sm_BT$cv_com<-NA
sm_BT$stability<-NA
sm_BT$vr_LdM<-NA # variance ratio by LdM, synchrony metric in [0,1]
sm_BT$TLintercept<-NA
sm_BT$TLslope.z<-NA
sm_BT$TLslope.z.lowCI<-NA
sm_BT$TLslope.z.upCI<-NA
sm_BT$pe_avg_cv<-NA
sm_BT$pe_mv<-NA

for(i in 1:nrow(sm_BT)){
  # read input data
  if(sm_BT$STUDY_ID[i]==sm_BT$newsite[i]){
    ns<-""
  }else{
    ns<-paste(sm_BT$newsite[i],"/",sep="")
  }
  mypath<-paste("../../DATA/for_BioTIME/wrangled_data/",
                sm_BT$REALM[i],"_plotlevel/",
                sm_BT$STUDY_ID[i],"/",ns,sep="")
  m<-readRDS(paste(mypath,"input_tailanal.RDS",sep=""))
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_BT$cv_com[i]<-cv_com
  sm_BT$stability[i]<- 1/cv_com
  sm_BT$vr_LdM[i]<-ecofolio::synchrony(m)
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_BT$TLintercept[i]<-res$c
  sm_BT$TLslope.z[i]<-res$z
  sm_BT$TLslope.z.lowCI[i]<-res$z.l
  sm_BT$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_BT$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_BT$pe_mv[i]<-res2
}
sm_BT$source<-"BioTIME"

#reorganize
sm_BT<-sm_BT%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                 initR,nsp,nyr,nind,npos,nL,nU,nneg,L,U,vr_LdM,
                                 avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                 stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                 TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_BT,"../../Results/for_BioTIME/stability_and_TL_est_insect.RDS")







