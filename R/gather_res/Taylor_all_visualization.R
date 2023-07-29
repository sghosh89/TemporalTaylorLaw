rm(list=ls())
library(tidyverse)
library(gridExtra)
`%notin%` <- Negate(`%in%`)
#======================= gather all results ======================
#----------------------------- for BioTIME -----------------------------------------------------
sm_BioTIME<-readRDS("../../Results/for_BioTIME/stability_and_TL_est_BioTIME.RDS")
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

checktab1<-sm_all%>%group_by(source,STUDY_ID)%>%summarise(nsites=n_distinct(newsite),
                                                          nuq_stab=n_distinct(stability))%>%ungroup()
# In the above table, whenever you find nsites!=nuq_stab, there are some duplicated sites
duplicated_tab<-checktab1%>%filter(nsites!=nuq_stab)

# So, we need to choose only any one site for a given STUDY_ID
sm_all_good<-sm_all%>%filter(STUDY_ID%notin%duplicated_tab$STUDY_ID)
sm_all_bad<-sm_all%>%filter(STUDY_ID%in%duplicated_tab$STUDY_ID)
sm_all_good2<-sm_all_bad%>%distinct(STUDY_ID,.keep_all = T)# distinct function always choose the first row

sm_all_good<-rbind(sm_all_good,sm_all_good2)
sm_all_good<-sm_all_good%>%arrange(source)#1722 observations in total

saveRDS(sm_all_good,"../../Results/gather_res/TaylorEstimate_alldata.RDS")

# total 1754 data
#-------------------------------------------------------------------
table(sm_all_good$REALM) # 110 freshwater, 1644 terrestrial
#------------------------ NOW DO PLOTTING ----------------------------------------
# first histogram of taylor's slope by REALM
range(sm_all_good$TLslope.z) # 0.529 to 3.5075759
# 1754 community data

br <- c(0,1,2,3,4)
# for terrestrial
gp1<-sm_all_good%>%
  ggplot( aes(TLslope.z, stat(density))) +
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 breaks=br) +
  geom_text(
    aes(label = round(stat(count) / sum((count)), 3)), 
    stat = 'bin', vjust = -0.5, breaks = br
  )+ ylim(c(0, 1))+xlab("Taylor's slope, z")+ylab("Density")+
  #scale_y_continuous(labels = scales::percent)+
  theme_bw()+theme(text = element_text(size = 14),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))
  
gp1

# now histogram of taylor's slope by TAXA
gp2 <- sm_all_good %>%
  ggplot( aes(TLslope.z, stat(density))) +
  geom_histogram(aes(y = stat(count) / sum(count)), 
                 breaks=br) +
  geom_text(
    aes(label = round(stat(count) / sum((count)), 3)), 
    stat = 'bin', vjust = -0.1, breaks = br
  )+
  facet_grid(~TAXA)+
  labs(fill="")+theme_classic()
gp2

#------------ Ok, now plot the iCV vs. species number for different z ---------
df<-sm_all_good%>%select(REALM,TAXA,TLslope.z,stability,nsp,pe_avg_cv,pe_mv)
df$zgt2<-ifelse(df$TLslope.z<2,0,1)
df$zgt2<-as.factor(df$zgt2)
  
gp3<-df%>%ggplot(aes(x=nsp,y=stability,fill=zgt2))+
  geom_point(alpha=0.3,shape=21)+ 
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')+
  geom_smooth(method="lm",se=T,col="black")+
  scale_fill_manual("",values=c("green","magenta"),labels = c("z<2", "z>2"))+
  ylab("Stability (=1/CV)")+xlab("Richness")+
  theme_bw()+theme(text = element_text(size = 14),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+ 
  theme(legend.position = c(0.84,0.13))
gp3

#=============================================================
# we will now plot the Portfolio Effect (PE): average CV and
# with mean variance scaling CV
#==============================================================
# for average CV and mean-var CV distribution
gp4<-ggplot(df,aes(y=pe_avg_cv,x=pe_mv,fill=zgt2))+
  geom_point(alpha=0.3,shape=21)+
  geom_abline(intercept = 0, slope = 1, color="gray55", 
              linetype="dotted", size=0.5)+
  #geom_smooth(method="lm",se=T,col="black")+
  scale_fill_manual("Taylor's slope",values=c("green","magenta"),labels = c("z<2", "z>2"))+
  xlim(c(0,13))+ylim(c(0,13))+
  xlab("Portfolio effect with mean-variance scaling")+ylab("Portfolio effect without mean-variance scaling")+
  theme_bw()+theme(text = element_text(size = 14),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  theme(legend.position = c(0.8,0.5))

df<-df[c("TLslope.z","pe_avg_cv","pe_mv")]
df<-gather(df, PE, value,2:3)
gp5<-df %>%
  ggplot( aes(x=value, fill=PE)) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity') +
  scale_fill_manual(values=c("red", "blue"),labels=c("Without mean-variance scaling",
                                                       "With mean variance scaling")) +
  labs(fill="")+xlab("Portfolio effect, PE")+ylab("Count")+
  theme_bw()+theme(text = element_text(size = 14),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  theme(legend.position = c(0.5,0.85))
gp5

gp1<-gp1+annotate(geom="text",x=Inf, y=Inf, label="A", size=7, vjust=2, hjust=2)
gp3<-gp3+annotate(geom="text",x=Inf, y=Inf, label="B", size=7, vjust=2, hjust=2)
gp5<-gp5+annotate(geom="text",x=Inf, y=Inf, label="C", size=7, vjust=2, hjust=2)
gp4<-gp4+annotate(geom="text",x=Inf, y=Inf, label="D", size=7, vjust=2, hjust=2)


pdf("../../Results/Prelim_res_plot/Taylor_all_visualization.pdf", width = 9, height = 9)
grid.arrange(gp1,gp3,gp5,gp4,nrow = 2)
dev.off()



