rm(list=ls())
#------------------------------------
source("tail_analysis.R")
#----------------------------------
library(tidyverse)
#---------------------------------------
resloc<-"../../Results/for_RivFishTIME/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}
#-----------------------------------------
# plot the sampling sites
good_TimeSeriesID_q3q4<-readRDS("../../DATA/for_RivFishTIME/wrangled_data/good_TimeSeriesID_q3q4.RDS")
x<-read.csv("../../DATA/for_RivFishTIME/raw_data/RivFishTIME_accessed_08dec2020/1873_2_RivFishTIME_SurveyTable.csv") # a dataframe
x_meta<-read.csv("../../DATA/for_RivFishTIME/raw_data/RivFishTIME_accessed_08dec2020/1873_2_RivFishTIME_TimeseriesTable.csv")
z<-x %>% distinct(TimeSeriesID, .keep_all = TRUE)%>%dplyr::select(TimeSeriesID,UnitAbundance)
x_meta<-inner_join(z,x_meta,by="TimeSeriesID")

x_meta<-x_meta%>%filter(TimeSeriesID%in%good_TimeSeriesID_q3q4)

library(maps)
wd<-map_data("world")
wd<-wd%>%filter(long<50 & lat>-50)
g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
g1<-g1+geom_polygon(data=wd, aes(x=long, y=lat, group=group), colour="gray90", fill="gray90")
g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
             legend.position="none",axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
g1<-g1+geom_point(data=x_meta,aes(y=Latitude,x=Longitude),color="blue",alpha=0.1)+
  ggtitle(paste("RivFishTIME: ",nrow(x_meta)," sites: min 20 years",sep=""))
g1
ggsave(paste(resloc,"sites_on_map.pdf",sep =""),
       width = 20, height = 10, units = "cm")
#=================== create results folder for each study sites ==================

for(i in 1:length(good_TimeSeriesID_q3q4)){
  k<-paste(resloc,good_TimeSeriesID_q3q4[i],sep="")
  if(!dir.exists(k)){
    dir.create(k)
  }
}

#------------ Now compute and plot the tail stats ---------------------

for(i in 1:length(good_TimeSeriesID_q3q4)){
  siteid<-good_TimeSeriesID_q3q4[i]
  resloc_output<-paste(resloc,siteid,"/",sep="")
  
  resloc_input<-paste("../../DATA/for_RivFishTIME/wrangled_data/",siteid,"/",sep="")
  df<-readRDS(paste(resloc_input,"commonspecies_timeseries.RDS",sep="")) # dataframe with species timeseries along column
  
  #----------- analysis without raresp ----------------
  id<-which(colnames(df)=="raresp")
  if(length(id)>0){
    df<-df[,-id]
  }
  res<-tail_analysis(mat = df, resloc = resloc_output, nbin = 2)
  cat("---------- i= ",i," siteid=",siteid," ----------\n")
}

#--------------------------------------------------------------------
# compute between years and between species spearman correlations (only significant)
# for each input year by species matrix

source("compute_avg_cor.R") # minimum 5 species, minimum 5 years 
# are required to estimate this metric, we are considering always >=20 years, and >=15 species 
# so it's okay for years
for(i in c(1:length(good_TimeSeriesID_q3q4))){
  
  siteid<-good_TimeSeriesID_q3q4[i]
  resloc_output<-paste(resloc,siteid,"/",sep="")
  
  resloc_input<-paste("../../DATA/for_RivFISHTIME/wrangled_data/",siteid,"/",sep="")
  df<-readRDS(paste(resloc_input,"commonspecies_timeseries.RDS",sep="")) # dataframe with species timeseries along column
  
  id<-which(colnames(df)=="raresp")
  if(length(id)>0){
    df<-df[,-id]
  }
  
  res<-compute_avg_cor(mat=df)
  saveRDS(res,paste(resloc_output,"avg_cor_metric.RDS",sep=""))
  
}

#--------------- Do a summary stats  ------------------
summary_table<-c()
for (i in c(1:length(good_TimeSeriesID_q3q4))){
  siteid<-good_TimeSeriesID_q3q4[i]
  resloc_input<-paste(resloc,siteid,"/",sep="")
  x<-readRDS(paste(resloc_input,"summary_df.RDS",sep=""))
  y<-readRDS(paste(resloc_input,"avg_cor_metric.RDS",sep=""))
  z<-cbind(x,y)
  summary_table<-rbind(summary_table,z)
}
summary_table<-cbind(siteid=good_TimeSeriesID_q3q4,summary_table)

# to get initial richness
for(i in 1:nrow(summary_table)){
  siteid<-summary_table$siteid[i]
  resloc_input<-paste("../../DATA/for_RivFishTIME/wrangled_data/",siteid,"/",sep="")
  bigM<-readRDS(paste(resloc_input,"allspecies_timeseries.RDS",sep=""))
  summary_table$initR[i]<-ncol(bigM)
}
# reorganize
summary_table<-summary_table%>%dplyr::select(siteid,initR,nsp,nyr,nint,nind,npos,nL,nU,nneg,L,U,
                                             avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp)

saveRDS(summary_table,"../../Results/for_RivFishTIME/summary_table.RDS")

summary_table<-inner_join(summary_table,x_meta,by=c("siteid"="TimeSeriesID"))
saveRDS(summary_table,"../../Results/for_RivFishTIME/summary_table_detail_version.RDS")




