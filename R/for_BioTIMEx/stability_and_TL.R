rm(list=ls())
library(tidyverse)
#===================
# now call the function
get_inputloc_table<-function(dataset_idset){
  STUDY_ID<-c()
  newsite<-c()
  inputloc<-c()
  resloc<-c()
  for(i in 1:length(dataset_idset)){
    dataset_id<-dataset_idset[i]
    inputmatpath<-paste("../../DATA/for_BioTIMEx/wrangled_data/",dataset_id,"/",sep="")
    if(!file.exists(paste("../../DATA/for_BioTIMEx/wrangled_data/",dataset_id,"/sitelist.RDS",sep=""))){
      readF<- list.files(inputmatpath, pattern = "_inputmatrix_tailanal.RDS", full.names = TRUE)
      outpath<-paste("../../Results/for_BioTIMEx/",dataset_id,"/",sep="")
      STUDY_ID<-c(STUDY_ID,dataset_id)
      newsite<-c(newsite,dataset_id)
      inputloc<-c(inputloc,readF)
      resloc<-c(resloc,outpath)
    }else{
      sitelist<-readRDS(paste("../../DATA/for_BioTIMEx/wrangled_data/",dataset_id,"/sitelist.RDS",sep=""))
      for(k in 1:length(sitelist)){
        readF<- list.files(inputmatpath, pattern = paste(sitelist[k],"_inputmatrix_tailanal.RDS",sep=""), full.names = TRUE)
        outpath<-paste("../../Results/for_BioTIMEx/",dataset_id,"/",sitelist[k],"/",sep="")
        STUDY_ID<-c(STUDY_ID,dataset_id)
        newsite<-c(newsite,sitelist[k])
        inputloc<-c(inputloc,readF)
        resloc<-c(resloc,outpath)
      }
    }
    print(i)
  }
  
  inputloc_table<-as.data.frame(cbind(STUDY_ID,newsite,inputloc,resloc))
  
  return(inputloc_table)
}

dataset_idset<-c("cumbrian_phyto",
                 "gross_2016","oneida_phytopl_1975")
inputloc_table<-get_inputloc_table(dataset_idset=dataset_idset)

# get initial richness
inputloc_table$tempoloc<-NA
for(i in 1:nrow(inputloc_table)){
  inputloc_table$tempoloc[i]<-paste(strsplit(inputloc_table$inputloc[i],"/")[[1]][1:6],collapse="/")
  inputloc_table$tempoloc[i]<-paste(inputloc_table$tempoloc[i],"/",sep="")
}

inputloc_table$initR<-NA
i<-1
#M<-readRDS(inputloc_table$inputloc[i])
bigM<-read.csv(paste(inputloc_table$tempoloc[i],"allrawdata.csv",sep=""))
bigM_BLEL<-bigM%>%filter(SITE=="BLEL")
inputloc_table$initR[i]<-length(unique(bigM_BLEL$ID_SPECIES))
i<-2
#M<-readRDS(inputloc_table$inputloc[i])
bigM_ESTH<-bigM%>%filter(SITE=="ESTH")
inputloc_table$initR[i]<-length(unique(bigM_ESTH$ID_SPECIES))
i<-3
#M<-readRDS(inputloc_table$inputloc[i])
bigM_SBAS<-bigM%>%filter(SITE=="SBAS")
inputloc_table$initR[i]<-length(unique(bigM_SBAS$ID_SPECIES))
i<-4
inputloc_table$initR[i]<-129 # I know it from file data wrangling/gross_2016.r
i<-5
bigM<-read.csv(paste(inputloc_table$tempoloc[i],"oneida_phytopl_1975_grouped_phytoplankton_list_1975to2013_BM.csv",sep=""))
bigM<-bigM%>%filter(include==1)
inputloc_table$initR[i]<-length(unique(bigM$species))
saveRDS(inputloc_table,"../../Results/for_BioTIMEx/inputloc_table.RDS")
#====================
summary_table<-c()
pathlist <- inputloc_table$resloc
for(i in 1:length(pathlist)){
  tempo<-readRDS(paste(pathlist[i],"summary_df.RDS",sep=""))
  y<-readRDS(paste(pathlist[i],"avg_cor_metric.RDS",sep=""))
  z<-cbind(tempo,y)
  summary_table<-rbind(summary_table,z)
}
summary_table<-cbind(STUDY_ID=inputloc_table$STUDY_ID,
                     newsite=inputloc_table$newsite,
                     initR=inputloc_table$initR,
                     summary_table)
summary_table$REALM<-NA
summary_table$TAXA<-NA
summary_table$ORGANISMS<-NA

id<-which(summary_table$STUDY_ID%in%c("cumbrian_phyto","oneida_phytopl_1975"))
summary_table$REALM[id]<-"Freshwater"
summary_table$TAXA[id]<-"Freshwater plants"
summary_table$ORGANISMS[id]<-"Phytoplankton"

id<-which(summary_table$STUDY_ID=="gross_2016")
summary_table$REALM[id]<-"Terrestrial"
summary_table$TAXA[id]<-"Terrestrial plants"
summary_table$ORGANISMS[id]<-"Plant"
saveRDS(summary_table,"../../Results/for_BioTIMEx/summary_table.RDS")
#==================================================
#library(devtools)
# you also need to have RTools installed 
#devtools::install_github('seananderson/ecofolio')
library(ecofolio)
#==================
sm_BTx<-summary_table
sm_BTx$cv_com<-NA
sm_BTx$stability<-NA
sm_BTx$vr_LdM<-NA # variance ratio by LdM, synchrony metric in [0,1]
sm_BTx$TLintercept<-NA
sm_BTx$TLslope.z<-NA
sm_BTx$TLslope.z.lowCI<-NA
sm_BTx$TLslope.z.upCI<-NA
sm_BTx$pe_avg_cv<-NA
sm_BTx$pe_mv<-NA

for(i in 1:nrow(sm_BTx)){
  # read input data
  m<-readRDS(inputloc_table$inputloc[i])
  id<-which(colnames(m)%in%c("raresp","covsp"))
  if(length(id)>0){
    m<-m[,-id]
  }
  # compute stability
  tot_biomass<-apply(m, MARGIN=1, FUN=sum)
  cv_com<-ecofolio::cv(tot_biomass)
  sm_BTx$cv_com[i]<-cv_com
  sm_BTx$stability[i]<- 1/cv_com
  sm_BTx$vr_LdM[i]<-ecofolio::synchrony(m)
  # now compute Taylor's slope
  res<-fit_taylor(x=m, ci = T, na.rm = F)
  sm_BTx$TLintercept[i]<-res$c
  sm_BTx$TLslope.z[i]<-res$z
  sm_BTx$TLslope.z.lowCI[i]<-res$z.l
  sm_BTx$TLslope.z.upCI[i]<-res$z.u
  res1<-pe_avg_cv(x=m, ci = F, na.rm = F) # default not detrended
  sm_BTx$pe_avg_cv[i]<-res1
  res2<-pe_mv(x=m, ci = F, na.rm = F) # default not detrended
  sm_BTx$pe_mv[i]<-res2
}
sm_BTx$source<-"BioTIMEx"
#reorganize
sm_BTx<-sm_BTx%>%dplyr::select(c(source,STUDY_ID,newsite,REALM,TAXA,ORGANISMS,
                                 initR,nsp,nyr,nind,npos,nL,nU,nneg,L,U,vr_LdM,
                                 avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp,
                                 stability,cv_com,TLslope.z,TLslope.z.lowCI,
                                 TLslope.z.upCI,TLintercept,pe_avg_cv,pe_mv))

saveRDS(sm_BTx,"../../Results/for_BioTIMEx/stability_and_TL_est_BioTIMEx.RDS")








