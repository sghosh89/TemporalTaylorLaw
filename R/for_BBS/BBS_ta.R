rm(list=ls())
#------------------------------------
source("tail_analysis.R")
#----------------------------------
library(tidyverse)
#---------------------------------------
resloc<-"../../Results/for_BBS/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}
#-----------------------------------------
# prepare the metadata
fshort_list<-readRDS("../../DATA/for_BBS/wrangled_data/sourcefile_list.RDS")
uroutes<-readRDS("../../DATA/for_BBS/wrangled_data/unique_routes_all.RDS")
uroutes<-data.frame(Country_State_Route=uroutes)
x_meta<-read.csv("../../DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/routes.csv")
x_meta<-x_meta%>%unite("Country_State_Route",CountryNum,StateNum,Route,sep="_")
metadata<-inner_join(uroutes,x_meta,by="Country_State_Route")%>%
                            rename(Stratum_code=Stratum)

bbs_strata1<-read.csv("../../DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/BBS_physiographic_strataname.csv")
bbs_strata1<-bbs_strata1%>%dplyr::select(Stratum_code=Stratum,Stratum_name=Name,Stratum_area=Area.Km2)
#bbs_strata2<-read.csv("../../DATA/for_BBS/raw_data/BBSdata_accessed_03dec2020/BBS_physiographic_strataname_statemap.csv")
metadata<-inner_join(metadata,bbs_strata1,by="Stratum_code")
saveRDS(metadata,"../../DATA/for_BBS/wrangled_data/unique_routes_all_metadata.RDS")

#=================== create results folder for each study sites/routes ==================

for(i in 1:nrow(uroutes)){
  k<-paste(resloc,uroutes$Country_State_Route[i],sep="")
  if(!dir.exists(k)){
    dir.create(k)
  }
}

#------------ Now compute and plot the tail stats ---------------------

for(i in 1:nrow(uroutes)){
  siteid<-uroutes$Country_State_Route[i]
  resloc_output<-paste(resloc,siteid,"/",sep="")
  
  resloc_input<-paste("../../DATA/for_BBS/wrangled_data/",siteid,"/",sep="")
  df<-readRDS(paste(resloc_input,"input_mat_for_tailanal.RDS",sep="")) # dataframe with species timeseries along column
  
  #----------- analysis without raresp ----------------
  id<-which(colnames(df)=="raresp")
  if(length(id)>0){
    df<-df[,-id]
  }
  res<-tail_analysis(mat = df, resloc = resloc_output, nbin = 2)
  cat("---------- i= ",i," routeid=",siteid," ----------\n")
}
#--------------------------------------------------------------------
# compute between years and between species spearman correlations (only significant)
# for each input year by species matrix

source("compute_avg_cor.R") # minimum 5 species, minimum 5 years 
# are required to estimate this metric, we are considering always >=20 years, and >=15 species 
# so it's okay for years
for(i in 1:nrow(uroutes)){
  siteid<-uroutes$Country_State_Route[i]
  resloc_output<-paste(resloc,siteid,"/",sep="")
  
  resloc_input<-paste("../../DATA/for_BBS/wrangled_data/",siteid,"/",sep="")
  df<-readRDS(paste(resloc_input,"input_mat_for_tailanal.RDS",sep="")) # dataframe with species timeseries along column

  id<-which(colnames(df)=="raresp")
  if(length(id)>0){
    df<-df[,-id]
  }
  res<-compute_avg_cor(mat=df)
  saveRDS(res,paste(resloc_output,"avg_cor_metric.RDS",sep=""))
  cat("---------- i= ",i," routeid=",siteid," ----------\n")
}

#--------------- Do a summary stats for all routes ------------------
summary_table<-c()
for(i in 1:nrow(uroutes)){
  siteid<-uroutes$Country_State_Route[i]
  resloc_input<-paste(resloc,siteid,"/",sep="")
  x<-readRDS(paste(resloc_input,"summary_df.RDS",sep=""))
  y<-readRDS(paste(resloc_input,"avg_cor_metric.RDS",sep=""))
  z<-cbind(x,y)
  summary_table<-rbind(summary_table,z)
}
summary_table<-cbind(siteid=uroutes$Country_State_Route,summary_table)

# to get initial richness
for(i in 1:nrow(summary_table)){
  siteid<-summary_table$siteid[i]
  resloc_input<-paste("../../DATA/for_BBS/wrangled_data/",siteid,"/",sep="")
  bigM<-readRDS(paste(resloc_input,"sourcefile.RDS",sep=""))
  summary_table$initR[i]<-length(unique(bigM$AOU))
}
# reorganize
summary_table<-summary_table%>%dplyr::select(siteid,initR,nsp,nyr,nint,nind,npos,nL,nU,nneg,L,U,
                                             avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp)

saveRDS(summary_table,"../../Results/for_BBS/summary_table.RDS")

summary_table<-inner_join(summary_table,metadata,by=c("siteid"="Country_State_Route"))
saveRDS(summary_table,"../../Results/for_BBS/summary_table_detail_version.RDS")


