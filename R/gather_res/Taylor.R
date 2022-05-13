rm(list=ls())
library(tidyverse)
#======================= gather all results ======================
#----------------------------- for BioTIME -----------------------------------------------------
sm_BioTIME<-readRDS("../../Results/for_BioTIME/stability_and_TL_est_insect.RDS")
#----------------------------------- for_BioTIMEx -------------------------------------------
sm_BioTIMEx<-readRDS("../../Results/for_BioTIMEx/stability_and_TL_est_BioTIMEx.RDS")
all(colnames(sm_BioTIME)==colnames(sm_BioTIMEx))==T
#------------------------------------- for BBS -----------------------------------------------
sm_BBS<-readRDS("../../Results/for_BBS/stability_and_TL_est_BBS.RDS")
all(colnames(sm_BioTIME)==colnames(sm_BBS))==T
#-------------------------------- for RivFishTIME ---------------------------------------------
sm_RF<-readRDS("../../Results/for_RivFishTIME/stability_and_TL_est_RF.RDS")
all(colnames(sm_BioTIME)==colnames(sm_RF))==T
#--------------------------------------- for swisslake phyto -----------------------------------------------------
sm_swisslake_phyto<-readRDS("../../Results/for_swisslake/stability_and_TL_est_swisslakephyto.RDS")
all(colnames(sm_BioTIME)==colnames(sm_swisslake_phyto))==T
#--------------------------------------- for swisslake zoo -----------------------------------------------------
sm_swisslake_zoo<-readRDS("../../Results/for_swisslake/stability_and_TL_est_swisslakezoo.RDS")
all(colnames(sm_BioTIME)==colnames(sm_swisslake_zoo))==T
#------------------------------- for insectRoel ----------------------------------------------
sm_insect<-readRDS("../../Results/for_insectRoel/stability_and_TL_est_insect.RDS")
all(colnames(sm_BioTIME)==colnames(sm_insect))==T
#-------- bind together -------------
sm_all<-rbind(sm_BioTIME,sm_BioTIMEx,sm_BBS,sm_RF,sm_swisslake_phyto,sm_swisslake_zoo,sm_insect)
sm_all$TAXA<-tolower(sm_all$TAXA)
sm_all$ORGANISMS<-tolower(sm_all$ORGANISMS)
sm_all$REALM<-as.factor(sm_all$REALM)
sm_all$TAXA<-as.factor(sm_all$TAXA)

# only considers significant estimates
sm_all$sig<-ifelse(
  (sm_all$TLslope.z.lowCI*sm_all$TLslope.z.upCI)>0,1,0)
sm_all<-sm_all%>%filter(sig==1) #all 1763 obs were significant

sm_all$avg_cor_neg_btw_sp[is.na(sm_all$avg_cor_neg_btw_sp)]<-0 # no neg cor found

saveRDS(sm_all,"../../Results/gather_res/TaylorEstimate_alldata.RDS")

# total 1763 data
#-------------------------------------------------------------------
table(sm_all$REALM) # 111 freshwater, 1652 terrestrial
#------------------------ NOW DO PLOTTING ----------------------------------------
# first histogram of taylor's slope by REALM
range(sm_all$TLslope.z) # 0.9057848 to 3.5075759
# 1763 community data

sm_all%>%group_by(REALM)%>%summarise(n=sum(TLslope.z>0 & TLslope.z<=2))
sm_all%>%group_by(REALM)%>%summarise(n=sum(TLslope.z>2 & TLslope.z<=4))

#p1 <- sm_all %>%
#  ggplot( aes(x=TLslope.z, fill=REALM)) +
#  geom_histogram( color="NA", alpha=0.6, position = 'identity') +
#  geom_vline(aes(xintercept = 1))+
#  geom_vline(aes(xintercept = 2))+
#  scale_fill_manual(values=c("skyblue", "green")) + facet_wrap(~REALM)+
#  labs(fill="")+theme_bw()
#p1 # 60% for freshwater and 90% for terrestrial distributions within slope [1,2]

br <- c(0,1,2,3,4)
# for terrestrial
p1T<-sm_all%>%filter(REALM=="Terrestrial")%>%
  ggplot( aes(TLslope.z, stat(density))) +
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 breaks=br, fill="green") +
  geom_text(
    aes(label = round(stat(count) / sum((count)), 3)), 
    stat = 'bin', vjust = -0.5, breaks = br
  )+ ylim(c(0, 1))+
  #scale_y_continuous(labels = scales::percent)+
  theme_classic()
p1T

p1F<-sm_all%>%filter(REALM=="Freshwater")%>%
  ggplot( aes(TLslope.z, stat(density))) +
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 breaks=br, fill="skyblue") +
  geom_text(
    aes(label = round(stat(count) / sum((count)), 3)), 
    stat = 'bin', vjust = -0.5, breaks = br
  )+ ylim(c(0, 1))+
  #scale_y_continuous(labels = scales::percent)+
  theme_classic()
p1F

# now histogram of taylor's slope by TAXA
# terrestrial taxa
p2T <- sm_all %>% filter(REALM=="Terrestrial")%>%
  ggplot( aes(TLslope.z, stat(density))) +
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 breaks=br, fill="green") +
  geom_text(
    aes(label = round(stat(count) / sum((count)), 3)), 
    stat = 'bin', vjust = -0.1, breaks = br
  )+
  facet_grid(~TAXA)+
  labs(fill="")+theme_classic()
p2T
#freshwater taxa
p2F <- sm_all %>% filter(REALM=="Freshwater")%>%
  ggplot( aes(TLslope.z, stat(density))) +
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 breaks=br, fill="skyblue") +
  geom_text(
    aes(label = round(stat(count) / sum((count)), 3)), 
    stat = 'bin', vjust = -0.1, breaks = br
  )+
  facet_grid(~TAXA)+
  labs(fill="")+theme_classic()
p2F

#------------ Ok, now plot the iCV vs. species number for different z ---------
#sm_all1 #2439 obs.
df<-sm_all%>%select(REALM,TAXA,TLslope.z,stability,nsp)
df_p0<-df%>%ggplot(aes(x=nsp,y=stability, col=REALM))+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+ 
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(values=c("skyblue", "green"))+ylab("Stability, iCV")+
  theme_classic()
df_p0

# first for 0<=z<=2
df<-sm_all%>%select(REALM,TAXA,TLslope.z,stability,nsp)%>%
                filter(TLslope.z>=0 & TLslope.z<=2)
#1676 obs
df_p1<-df%>%ggplot(aes(x=nsp,y=stability, col=REALM))+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
    scale_color_manual(values=c("skyblue", "green"))+ylab("Stability, iCV")+ 
    theme_classic()
df_p1

# now, z in [2,4]
df<-sm_all%>%select(REALM,TAXA,TLslope.z,stability,nsp)%>%
  filter(TLslope.z>2 & TLslope.z<=4)
df_p2<-df%>%ggplot(aes(x=nsp,y=stability, col=REALM))+
  geom_jitter(position = position_jitter(width = 0, height = 0), alpha=0.3)+ 
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(values=c("skyblue", "green"))+ylab("Stability, iCV")+ 
  theme_classic()
df_p2

#---- plot aquatic data only, for 0<z<=2 and 2<z<=4 region in same plot ------

df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,nsp)%>%
    filter(TLslope.z>0 & TLslope.z<=4 & REALM=="Freshwater")
# 111 observations

# creating two categories for two regions
df$group<-ifelse(df$TLslope.z<=2,0,1) 
p1<-ggplot(df,aes(x=nsp,y=stability, col=as.factor(group), group=as.factor(group)))+
  geom_jitter(position = position_jitter(width = 0.15, height = 0.15), 
              alpha=0.3)+
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(labels = c("z<2", "z>2"), values=c("skyblue", "blue"))+ 
  theme_classic()+labs(title = "Freshwater", x = "nsp", y = "stability, iCV", color = "Taylor's Slope")
p1

df<-sm_all%>%select(REALM,TAXA,TLslope.z,cv_com,stability,nsp)%>%
  filter(TLslope.z>0 & TLslope.z<=4 & REALM=="Terrestrial")
# 1652 observations

# creating two categories for two regions
df$group<-ifelse(df$TLslope.z<=2,0,1) 
p1<-ggplot(df,aes(x=nsp,y=stability, col=as.factor(group), group=as.factor(group)))+
  geom_jitter(position = position_jitter(width = 0.01, height = 0.05), 
              alpha=0.3)+
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T)+
  scale_color_manual(labels = c("z<2", "z>2"),values=c("green", "green4"))+
  labs(title = "Terrestrial", x = "nsp", y = "stability, iCV", color = "Taylor's Slope")+
  theme_classic()
p1

#==============================================================
# we will now plot the Portfolio Effect (PE): average CV and
# with mean variance scaling CV
#==============================================================
 # for average CV and mean-var CV distribution
#plot(sm_all1$pe_avg_cv,sm_all1$pe_mv)
library(reshape2)
id<-which(colnames(sm_all)%in%c("pe_avg_cv","pe_mv","STUDY_ID","newsite"))
df<-sm_all%>%select(id)
df$idvar<-paste(df$STUDY_ID,df$newsite,sep=".")
df<-df%>%select(idvar,pe_avg_cv,pe_mv)
df = melt(df,id.vars = "idvar")
pg1<-ggplot(data = df) +
  geom_histogram(aes(x = value, y=(..count..)/sum(..count..), fill=variable), 
                 alpha=0.3, binwidth=2, position="identity")

br <- c(0:14)
dfT<-sm_all%>%filter(REALM=="Terrestrial")
ggplot(data = dfT) +
  geom_histogram(aes(x = pe_avg_cv, y=(..count..)/sum(..count..)), 
                 alpha=0.4, fill ="green",breaks=br) +
  geom_histogram(aes(x = pe_mv, y=(..count..)/sum(..count..)), 
                 alpha=0.4, fill ="green4",breaks=br)+
  xlab("Portfolio Effct")+ylab("Frequency")+ylim(c(0, 1))+theme_classic()
  
dfF<-sm_all%>%filter(REALM=="Freshwater")
ggplot(data = dfF) +
  geom_histogram(aes(x = pe_avg_cv, y=(..count..)/sum(..count..)), 
                 alpha=0.4, fill ="blue",breaks=br) +
  geom_histogram(aes(x = pe_mv, y=(..count..)/sum(..count..)), 
                 alpha=0.4, fill ="skyblue",breaks=br)+
  xlab("Portfolio Effct")+ylab("Frequency")+ylim(c(0, 1))+theme_classic()








