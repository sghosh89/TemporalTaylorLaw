rm(list=ls())
#------------------------------------
source("tail_analysis.R")
source("compute_avg_cor.R")
#----------------------------------
library(tidyverse)
#---------------------------------------
resloc<-"../../Results/for_swisslake/"
if(!dir.exists(resloc)){
  dir.create(resloc)
}
#--------------------------------------------

# ============= for lake 1 =======================
resloc1<-"../../Results/for_swisslake/L1WA/"
if(!dir.exists(resloc1)){
  dir.create(resloc1)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L1WA.RDS")
tail_analysis(mat=mat, resloc=resloc1, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc1,"avg_cor_metric.RDS",sep=""))

# =========== for lake 2 ===========================
resloc2<-"../../Results/for_swisslake/L2UZ/"
if(!dir.exists(resloc2)){
  dir.create(resloc2)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L2UZ.RDS")
tail_analysis(mat=mat, resloc=resloc2, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc2,"avg_cor_metric.RDS",sep=""))

# =========== for lake 3 ===========================
resloc3<-"../../Results/for_swisslake/L3LU/"
if(!dir.exists(resloc3)){
  dir.create(resloc3)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L3LU.RDS")
tail_analysis(mat=mat, resloc=resloc3, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc3,"avg_cor_metric.RDS",sep=""))

# =========== for lake 4 ===========================
resloc4<-"../../Results/for_swisslake/L4LZ/"
if(!dir.exists(resloc4)){
  dir.create(resloc4)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L4LZ.RDS")
tail_analysis(mat=mat, resloc=resloc4, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc4,"avg_cor_metric.RDS",sep=""))

# =========== for lake 5 ===========================
resloc5<-"../../Results/for_swisslake/L5SE/"
if(!dir.exists(resloc5)){
  dir.create(resloc5)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L5SE.RDS")
tail_analysis(mat=mat, resloc=resloc5, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc5,"avg_cor_metric.RDS",sep=""))

# =========== for lake 6 ===========================
resloc6<-"../../Results/for_swisslake/L6HA/"
if(!dir.exists(resloc6)){
  dir.create(resloc6)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L6HA.RDS")
tail_analysis(mat=mat, resloc=resloc6, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc6,"avg_cor_metric.RDS",sep=""))

# =========== for lake 7 ===========================
resloc7<-"../../Results/for_swisslake/L7BA/"
if(!dir.exists(resloc7)){
  dir.create(resloc7)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L7BA.RDS")
tail_analysis(mat=mat, resloc=resloc7, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc7,"avg_cor_metric.RDS",sep=""))

# =========== for lake 8 ===========================
resloc8<-"../../Results/for_swisslake/L8GR/"
if(!dir.exists(resloc8)){
  dir.create(resloc8)
}

mat<-readRDS("../../DATA/for_swisslake/wrangled_data/input_mat_for_tail_analysis_L8GR.RDS")
tail_analysis(mat=mat, resloc=resloc8, nbin=2)
res<-compute_avg_cor(mat=mat)
saveRDS(res,paste(resloc8,"avg_cor_metric.RDS",sep=""))

#####################################################################################
# Now, do the summary results

resloc_list<-c(resloc1,resloc2,resloc3,resloc4,resloc5,resloc6,resloc7,resloc8)
summary_df<-c()
for(i in 1:8){
  resl<-resloc_list[i]
  df<-readRDS(paste(resl,"summary_df.RDS",sep=""))
  y<-readRDS(paste(resl,"avg_cor_metric.RDS",sep=""))
  z<-cbind(df,y)
  summary_df<-rbind(summary_df,z)
}
summary_df$siteid<-c("L1WA","L2UZ","L3LU","L4LZ","L5SE","L6HA","L7BA","L8GR")

summary_df$initR<-NA
for(i in 1:nrow(summary_df)){
  bigM<-read.csv(paste("../../DATA/for_swisslake/wrangled_data/species_list_",summary_df$siteid[i],"_c_blake_sorted.csv",sep=""))
  bigM<-bigM%>%filter(include==1)
  summary_df$initR[i]<-length(unique(bigM$species))
}
# reorganize
summary_df<-summary_df%>%dplyr::select(siteid,initR,nsp,nyr,nint,nind,npos,nL,nU,nneg,L,U,
                                             avg_cor_btw_yr,avg_cor_pos_btw_sp,avg_cor_neg_btw_sp)

saveRDS(summary_df,"../../Results/for_swisslake/summary_table_phytoplankton.RDS")












