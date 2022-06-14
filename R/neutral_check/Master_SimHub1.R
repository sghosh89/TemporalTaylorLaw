rm(list=ls())
###################
resloc<-"../../Results/neutral_check/"
if(!dir.exists(resloc)){dir.create(resloc)}
#####################################################
set.seed(seed=123)
#####################################################
source("./mySimHub1.R")

nrep=100
S=30
j=100
cycles=1e7
len=1e2
resloc<-paste("../../Results/neutral_check/mySimHub1_S_",S,"_j_",j,"_cycles_",cycles,"/",sep="")
if(!dir.exists(resloc)){dir.create(resloc)}

call_mySimHub1(resloc=resloc,nrep=nrep,S=S,j=j,cycles=cycles,len=len)

####################################

# do tail dep analysis, avg cor for each replicate of neutral model
source("taildep_neutral.R")
resloc_ta<-paste(resloc,"tail_analysis_res/",sep="")
if(!dir.exists(resloc_ta)){dir.create(resloc_ta)}
resloc_input<-paste("../../Results/neutral_check/mySimHub1_S_",S,"_j_",j,"_cycles_",cycles,"/",sep="")

taildep_neutral(resloc_ta = resloc_ta, resloc_input = resloc_input, nrep=nrep)

####################################

# now summarize res for all replicates
summary_table<-c()
for(i in 1:nrep){
  resloc_input<-paste(resloc_ta,i,"/",sep="")
  x<-readRDS(paste(resloc_input,"summary_df.RDS",sep=""))
  y<-readRDS(paste(resloc_input,"avg_cor_metric.RDS",sep=""))
  z<-cbind(x,y)
  summary_table<-rbind(summary_table,z)
}
summary_table<-cbind(nrep=1:nrep,summary_table)
saveRDS(summary_table,paste(resloc,"summary_table.RDS",sep=""))

####################################
# compare stats with empirical data
library(tidyverse)
######################
sm_all<-summary_table
range(sm_all$avg_cor_btw_yr) # 0.461105 0.762551
sm_all$net_taildep<-sm_all$L+sm_all$U
sm_all$xx<-sm_all$net_taildep/sm_all$nint
hist(sm_all$xx, 100, xlim=c(-0.2,0.1))

# from empirical data (TemporalTaylorLaw folder)
emp_df<-readRDS("../../Results/gather_res/TaylorEstimate_alldata.RDS")
emp_df$net_taildep<-emp_df$L+emp_df$U
emp_df$nint<-emp_df$nyr*(emp_df$nyr - 1)*0.5
emp_df$yy<-emp_df$net_taildep/emp_df$nint

df<-emp_df%>%filter(nsp==30) #21 obs
range(df$avg_cor_btw_yr) #0.5656070 0.8494451
range(df$nyr) # 20 26
hist(df$yy, 100, col="black", add=T)

#make comparison between two distributions:
kst<-ks.test(sm_all$xx,df$yy) # they are different

gp<-ggplot() +
  geom_histogram(aes(xx, fill = "neutral"),colour="black", alpha = .2, data = sm_all)+
  geom_density(aes(xx, fill = "neutral"), alpha = .2, data = sm_all) +
  geom_histogram(aes(yy, fill = "empirical"),colour="black", alpha = .2, data = df)+
  geom_density(aes(yy, fill = "empirical"), alpha = .2, data = df) +
  scale_fill_manual(name = paste("k=",round(kst$statistic,2),", p=",round(kst$p.value,2),sep="") , 
                    values = c(neutral = "gray", empirical = "purple"))+
  xlab("net_tail_dep/nint")+xlim(-0.3,0.1)+
  theme_classic()
gp
ggsave(paste(resloc,"compare_distribn.pdf",sep=""), height=6,width=6)

##############################
# ok, now check the rank abundance curve for all sp 
# you sampled in the last 26 timesteps
source("check_rank_abund.R")
check_rank_abund(resloc=resloc)
##############################
# we also repeated the analysis with S=30, j=100, cycles=1e5
# to check the robustness of the results

# AND we also repeated the analysis with S=30, j=100, cycles=1e7
# to check the robustness of the results









