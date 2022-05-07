rm(list=ls())
#------------------------------------
source("tail_analysis.R")
source("compute_avg_cor.R")
#----------------------------------
library(tidyverse)
#--------- read data ---------------------------------
xm<-read.csv("../../DATA/for_insectRoel/20yrFreshwater_Metadata.csv")
x<-readRDS("../../DATA/for_insectRoel/20yrFreshwaterData 202106.rds")
#x<-readRDS("../../DATA/for_insectRoel/20yrFreshwaterData.rds")
#x<-x%>%filter(Number>0)
xtbl<-x%>%dplyr::distinct(Rank,.keep_all=T)
xtbl<-xtbl[order(xtbl$Rank),]
xtbl<-xtbl%>%dplyr::select(Rank,Level)
#---------------------------------------
resloc<-"../../Results/for_insectRoel/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}

Datasource_ID<-setdiff(sort(unique(x$Datasource_ID)),63)
bad_did<-c()
for(i in 1:length(Datasource_ID)){
  did<-paste(resloc,Datasource_ID[i],"/",sep="")
  if(!dir.exists(did)){
    dir.create(did)
  }
  pidlist<-readRDS(paste("../../DATA/for_insectRoel/wrangled_data/",Datasource_ID[i],"/Plot_ID_list.RDS",sep=""))
  badpidlist<-readRDS(paste("../../DATA/for_insectRoel/wrangled_data/",Datasource_ID[i],"/bad_pidlist.RDS",sep=""))
  goodpidlist<-setdiff(pidlist,badpidlist)
  if(length(goodpidlist)!=0){
    saveRDS(goodpidlist,paste("../../Results/for_insectRoel/",Datasource_ID[i],"/goodpidlist.RDS",sep=""))
    for(j in 1:length(goodpidlist)){
      pid<-paste(did,goodpidlist[j],"/",sep="")
      if(!dir.exists(pid)){
        dir.create(pid)
      }
    }
  }else{
    bad_did<-c(bad_did,Datasource_ID[i])
  }
}
saveRDS(bad_did,paste("../../Results/for_insectRoel/bad_did_",bad_did,".RDS",sep=""))

#----------- Now compute and plot the tail stats ---------------------
Datasource_ID<-setdiff(Datasource_ID,bad_did)
for(i in 1:length(Datasource_ID)){
    did<-Datasource_ID[i]
    goodpidlist<-readRDS(paste("../../Results/for_insectRoel/",did,"/goodpidlist.RDS",sep=""))
    for(j in 1:length(goodpidlist)){
      pid<-goodpidlist[j]
      resloc_output<-paste(resloc,did,"/",pid,"/",sep="")
      resloc_input<-paste("../../DATA/for_insectRoel/wrangled_data/",did,"/",pid,"/",sep="")
      
      df<-readRDS(paste(resloc_input,"inputmat_for_tailanal.RDS",sep="")) # dataframe with species timeseries along column
      #----------- tail analysis ----------------
      
      res<-tail_analysis(mat = df, resloc = resloc_output, nbin = 2)
      cat("------- i= ",i," j=",j," did =", did, " pid = ",pid," ----------\n")
      resxx<-compute_avg_cor(mat=df)
      saveRDS(resxx,paste(resloc_output,"avg_cor_metric.RDS",sep=""))
      
    }
}

#--------------- Do a summary stats for all ------------------
summary_table<-c()
didlist<-c()
pidlist<-c()
for(i in 1:length(Datasource_ID)){
  did<-Datasource_ID[i]
  goodpidlist<-readRDS(paste("../../Results/for_insectRoel/",did,"/goodpidlist.RDS",sep=""))
  for(j in 1:length(goodpidlist)){
  pid<-goodpidlist[j]
  didlist<-c(didlist,did)
  pidlist<-c(pidlist,pid)
  resloc_input<-paste("../../Results/for_insectRoel/",did,"/",pid,"/",sep="")
  st<-readRDS(paste(resloc_input,"summary_df.RDS",sep=""))
  y<-readRDS(paste(resloc_input,"avg_cor_metric.RDS",sep=""))
  z<-cbind(st,y)
  summary_table<-rbind(summary_table,z)
  }
}
summary_table<-cbind(STUDY_ID=didlist,newsite=pidlist,summary_table)
saveRDS(summary_table,"../../Results/for_insectRoel/summary_table.RDS")

source("get_operational_richness.R")
summary_table$initR<-NA
for(i in 1:nrow(summary_table)){
  initR<-readRDS(paste("../../DATA/for_insectRoel/wrangled_data/",
                       summary_table$STUDY_ID[i],"/",summary_table$newsite[i],
                       "/initial_richness.RDS",sep=""))
  summary_table$initR[i]<-initR
}

# reorganize
summary_table<-summary_table%>%dplyr::select(STUDY_ID,newsite,initR,nsp,nyr,nint,nind,npos,nL,nU,nneg,L,U,
                                             avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp)
saveRDS(summary_table,"../../Results/for_insectRoel/summary_table.RDS")


metadata<-xm%>%dplyr::select(Plot_ID,REALM=Realm,TAXA=Taxonomic_scope,ORGANISMS=Taxonomic_scope,Latitude,Longitude)
length(unique(xm$Plot_ID))==nrow(xm)
summary_table<-inner_join(summary_table,metadata,by=c("newsite"="Plot_ID"))
summary_table$TAXA<-"Freshwater invertebrates"
#summary_table<-summary_table%>%filter(f_nind!=1)
saveRDS(summary_table,"../../Results/for_insectRoel/summary_table_detail_version.RDS")

#################################################################################################
