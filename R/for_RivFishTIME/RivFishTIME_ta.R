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







