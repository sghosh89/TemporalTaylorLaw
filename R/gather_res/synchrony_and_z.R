# Ok, we have found earlier that stability-diversity relationship is affected 
# by heterogeniety in z
# But what about stability-synchrony relationship?

# so, we will start with some data visualization
# plot for stability vs. variance ratio (0:asynchrony to 1:synchrony) 
# for different z region (0<=z<2, 2<=z<=4)

rm(list=ls()) 
library(tidyverse)
#==========================
sm_all<-readRDS("../../Results/gather_res/TaylorEstimate_alldata.RDS")
sm_all$pe_diff<-sm_all$pe_avg_cv - sm_all$pe_mv
sm_all$TA<-sm_all$L+abs(sm_all$U)#total tail asymmetry

# --------------------- for overall syn between species ------------------------
df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,vr_LdM)
df_p0<-df%>%ggplot(aes(x=vr_LdM,y=stability, col=REALM))+
  geom_jitter(position = position_jitter(width = 0.02, height = 0.02), alpha=0.3)+ 
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(values=c("skyblue", "green"))+
  ylab("Stability, iCV")+xlab("Synchrony between sp, VR")+
  theme_classic()
df_p0

# first for 0<=z<2
#df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,vr_LdM)%>%filter(TLslope.z>=0 & TLslope.z<2)
#df_p1<-df%>%ggplot(aes(x=vr_LdM,y=stability, col=REALM))+
#  geom_jitter(position = position_jitter(width = 0.02, height = 0.02), alpha=0.3)+
#  scale_x_continuous(trans = 'log2') +
#  scale_y_continuous(trans = 'log2')+
#  geom_smooth(method="lm",se=T)+
#  scale_color_manual(values=c("skyblue", "green"))+
#  ylab("Stability, iCV")+xlab("Synchrony between sp, VR")+
#  theme_classic()
#df_p1

# now, z in [2,4]
#df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,vr_LdM)%>%filter(TLslope.z>=2 & TLslope.z<4)
#df_p2<-df%>%ggplot(aes(x=vr_LdM,y=stability, col=REALM))+
#  geom_jitter(position = position_jitter(width = 0.02, height = 0.02), alpha=0.3)+ 
#  scale_x_continuous(trans = 'log2') +
#  scale_y_continuous(trans = 'log2')+
#  geom_smooth(method="lm",se=T)+
#  scale_color_manual(values=c("skyblue", "green"))+
#  ylab("Stability, iCV")+xlab("Synchrony between sp, VR")+
 # theme_classic()
#df_p2
sm_all$asyn_by_syn<-abs(sm_all$avg_cor_neg_btw_sp)/sm_all$avg_cor_pos_btw_sp
df<-sm_all%>%select(REALM,TAXA,TLslope.z,avg_cor_btw_yr,vr_LdM,asyn_by_syn)
df_p0<-df%>%ggplot(aes(x=asyn_by_syn,y=TLslope.z, col=REALM))+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+ 
  #scale_x_continuous(trans = 'log2') +
  #scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(values=c("skyblue", "green"))+ylab("z")+
  theme_classic()
df_p0 # same pattern as VR_LdM vs z plot

# ---------------- for tail dep. syn between years-------------------------

df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,TA)
df_p0<-df%>%ggplot(aes(x=TA,y=stability, col=REALM))+
  geom_jitter(position = position_jitter(width = 0.02, height = 0.02), alpha=0.3)+ 
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(values=c("skyblue", "green"))+
  ylab("Stability, iCV")+xlab("Synchrony between years, Tail-dep.")+
  theme_classic()
df_p0

# first for 0<=z<2
#df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,TA)%>%filter(TLslope.z>=0 & TLslope.z<2)
#df_p1<-df%>%ggplot(aes(x=TA,y=stability, col=REALM))+
#  geom_jitter(position = position_jitter(width = 0.02, height = 0.02), alpha=0.3)+
#  scale_x_continuous(trans = 'log2') +
#  scale_y_continuous(trans = 'log2')+
#  geom_smooth(method="lm",se=T)+
#  scale_color_manual(values=c("skyblue", "green"))+
#  ylab("Stability, iCV")+xlab("Synchrony, Tail-dep.")+
#  theme_classic()
#df_p1

# now, z in [2,4]
#df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,TA)%>%filter(TLslope.z>=2 & TLslope.z<4)
#df_p2<-df%>%ggplot(aes(x=TA,y=stability, col=REALM))+
#  geom_jitter(position = position_jitter(width = 0.02, height = 0.02), alpha=0.3)+ 
#  scale_x_continuous(trans = 'log2') +
#  scale_y_continuous(trans = 'log2')+
#  geom_smooth(method="lm",se=T)+
#  scale_color_manual(values=c("skyblue", "green"))+
#  ylab("Stability, iCV")+xlab("Synchrony, tail-dep.")+
#  theme_classic()
#df_p2
#=====================================================
# difference btw patterns detected!!!
#=====================================================
# my hypothesis: z can vary depending on tail-dep btw years 
# irrespective of synchrony/asynchrony ratio between species

# ----------- avg cor btw years vs z plot --------------------------
p1<-sm_all%>%ggplot(aes(y=avg_cor_btw_yr,x=TLslope.z, col=REALM))+
  #ylim(0,4)+xlim(0,1)+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+
  geom_smooth(method="lm")+geom_vline(xintercept=2,linetype='dotted', col = 'gray4')+
  scale_color_manual(values=c("skyblue", "green"))+
  xlab("Taylor's slope, z")+ylab("rho, avg. cor btw years")+
  theme_classic()
p1

# ----------- synchrony vs z plot --------------------------

# now plot histogram of synchrony vs. z
df1<-sm_all%>%select(REALM,TAXA,TLslope.z,stability,vr_LdM)%>%
  filter(TLslope.z>=0 & TLslope.z<2)
df2<-sm_all%>%select(REALM,TAXA,TLslope.z,stability,vr_LdM)%>%
  filter(TLslope.z>=2 & TLslope.z<=4)

p1<-sm_all%>%ggplot(aes(y=vr_LdM,x=TLslope.z, col=REALM))+
  #ylim(0,4)+xlim(0,1)+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+
  geom_smooth(method="lm")+geom_vline(xintercept=2,linetype='dotted', col = 'gray4')+
  scale_color_manual(values=c("skyblue", "green"))+
  xlab("Taylor's slope, z")+ylab("Synchrony btw sp, VR")+
  theme_classic()
p1 # it showed consistent result as predicted by Ives 2003 Nature paper
# with incresing asynchrony 1<z<2 for terrestrial, freshwater shows 
# lots of variation and no consistent pattern.

#sm_all$asyn_by_syn<-abs(sm_all$avg_cor_neg_btw_sp)/sm_all$avg_cor_pos_btw_sp
#df<-sm_all%>%select(REALM,TAXA,TLslope.z,avg_cor_btw_yr,L_by_U,vr_LdM,asyn_by_syn)
#df_p0<-df%>%ggplot(aes(x=asyn_by_syn,y=TLslope.z, col=REALM))+
#  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+ 
#    geom_smooth(method="lm",se=T)+
#  scale_color_manual(values=c("skyblue", "green"))+ylab("z")+
#  theme_classic()
#df_p0 # same pattern as VR_LdM vs z plot

#--------- let's see what it would look for tail-dependent synchrony ----------
p2<-sm_all%>%ggplot(aes(y=TA,x=TLslope.z, col=REALM))+
  #ylim(0,4)+xlim(0,15)+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+
  geom_smooth(method="lm")+geom_vline(xintercept=2,linetype='dotted', col = 'gray4')+
  scale_color_manual(values=c("skyblue", "green"))+
  xlab("Taylor's slope, z")+ylab("Synchrony between years, Tail-dependent")+
  #scale_y_continuous(trans = 'log2') +
  theme_classic()
p2
# so for freshwater with increasing tail-dep synchrony z sharply decreases.
# striking difference between patterns of z vs synchrony for both realm:
# z increases as overall syn increases both realm
# z decreases as tail dep syn increases for freshwater, for terrestrial it's same

# curious? check distribution of L/U for both realm
sm_all$L_by_U<-sm_all$L/abs(sm_all$U)
range(sm_all$L_by_U) # we checked, most of the time it's upper tail dep. for both realms

sm_all_T<-sm_all%>%filter(REALM=="Terrestrial")
sm_all_F<-sm_all%>%filter(REALM=="Freshwater")
hist(sm_all_T$nyr,col=rgb(0,1,0,0.3),breaks=100)
hist(sm_all_F$nyr,col=rgb(0,0,1,0.3),add=T,breaks=100)
#=======================================
# I want to test something
# say, if 0.6<=rho<=0.7, rho= between-year avg cor.,
# then visualize how Terrestrial vs. Freshwater tail dep. syn 
# spread

sm_all$net_taildep<-sm_all$L+sm_all$U
mydat<- sm_all%>%filter(avg_cor_btw_yr>=0.5 & avg_cor_btw_yr<=0.7)
# also keep a common range of species
mydat<-mydat%>%filter(nsp<=50)
mydat_T<-mydat%>%filter(REALM=="Terrestrial") #756 obs
mydat_F<-mydat%>%filter(REALM=="Freshwater") #85 obs

hist(abs(mydat_F$U), breaks=30,col=rgb(0,0,1,0.3), 
     ylim=c(0,20), main="Freshwater", xlab="L, |U|")
hist(mydat_F$L, breaks=30,col=rgb(1,0,0,0.3), add=T)

hist(abs(mydat_T$U), breaks=30,col=rgb(0,0,1,0.3), 
     ylim=c(0,20), main="Terrestrial", xlab="L, |U|")
hist(mydat_T$L, breaks=30,col=rgb(1,0,0,0.3), add=T)

# this shows terrestrial are highly upper tail dep, freshwater is still balanced

hist(mydat_T$avg_cor_btw_yr, breaks=30,col=rgb(0,1,0,0.3), 
     ylim=c(0,20), main="", xlab="rho")
hist(mydat_F$avg_cor_btw_yr, breaks=30,col=rgb(0,0,1,0.3), add=T)

# this shows terrestrial are highly upper tail dep
# but we need some rarefaction/randomization to compare equal number of sample from the above data
source("test_empirical_z_taildep.R")






