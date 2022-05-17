rm(list=ls()) 
library(tidyverse)
#==========================
sm_all<-readRDS("../../Results/gather_res/TaylorEstimate_alldata.RDS")
sm_all$net_taildep<-sm_all$L+sm_all$U
sm_all$TA<-sm_all$L+abs(sm_all$U)
sm_all$nint<-sm_all$nyr*(sm_all$nyr - 1)*0.5
sm_all$rank_td<-rank(sm_all$net_taildep/sm_all$nint) # lowest value= rank 1

# see how nsp across realm distributed
ggplot() + 
  geom_density(data=sm_all, aes(x=nsp, group=REALM, fill=REALM),alpha=0.5, adjust=2) + 
  xlab("Richness") +geom_vline(xintercept=45)+ 
  ylab("Density")+scale_fill_manual(values=c("skyblue","green"))+
  theme_classic()

# So, I choose range 15 to 45 (until the vertical line), avg nsp=30

library(RColorBrewer)
nsplist<-c(15:45)
df<-sm_all%>%filter(nsp%in%nsplist)
table(df$REALM)
range(df$avg_cor_btw_yr) #0.4354706 0.8989881
ggplot() + 
  geom_density(data=df, aes(x=avg_cor_btw_yr, group=REALM, fill=REALM),alpha=0.5, adjust=2) + 
  xlab("rho") +geom_vline(xintercept=0.5)+
  ylab("Density")+scale_fill_manual(values=c("skyblue","green"))+
  theme_classic()


#df<-df%>%filter(avg_cor_btw_yr>=0.5)
#df_T<-df%>%filter(REALM=="Freshwater")

#range(df$avg_cor_btw_yr) #0.5241136 0.8989881

mean(df$nsp) #~30 sp
mean(df$nyr)# ~26 yr

# considering by species
ggplot(df, aes(x = net_taildep/nint, y = TLslope.z, 
               color=as.factor(nsp), group=as.factor(nsp))) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=F)+facet_wrap(~nsp)+
  #scale_color_brewer(palette = "Reds")+
  #labs(color="Richness")+
  theme_classic()+theme(legend.position = "none")

# considering by groups: 2 REALM
ggplot(df, aes(x = net_taildep/nint, y = TLslope.z,
               color=REALM, group=REALM)) + 
  geom_point(alpha=0.3,col=as.factor(df$nsp)) + 
  geom_smooth(method="lm",se=T, aes(fill=REALM))+#facet_grid(vars(REALM))+
  #scale_color_gradient(low = "white", high = "red")+
  scale_color_manual(values=c("skyblue","green"))+
  scale_fill_manual(values=c("skyblue","darkolivegreen2"))+
  theme_classic()#+theme(legend.position = "none")

# considering all point
ggplot(df, aes(x = net_taildep/nint, y = TLslope.z
               )) + 
  geom_point(alpha=0.3,col=as.factor(df$nsp)) + 
  geom_smooth(method="lm",se=T,col="red",aes(fill="red"))+#facet_grid(vars(REALM))+
  theme_classic()+theme(legend.position = "none")

# make plot for z vs. net_taildep/nint for different range of species
df$bins<-cut(df$nsp, breaks=c(14,20,25,30,35,45)) # make range 
df0<-df%>%select(nsp,avg_cor_btw_yr,nyr,TLslope.z,REALM,bins)
tbl<-df0%>%group_by(bins)%>%summarise(avg_z=mean(TLslope.z))%>%ungroup()
tbl
plot(tbl$bins,tbl$avg_z)
ggplot(df0, aes(x=bins, y=TLslope.z, col=factor(bins),alpha=0.3)) + 
  geom_point() +
  geom_point(data=tbl, aes(x=bins, y=avg_z), color='black', cex=3, pch=15)+
  xlab("Richness")+
  theme_classic()+theme(legend.position = "none")



############### EXTRA PLOT ###########################
# now plot z vs. rho for both realm with error bar
df$color<-ifelse(df$REALM=="Terrestrial","green","skyblue")

op<-par(mar=c(3,3,1,1),mgp=c(2,1,0))
plot(df$avg_cor_btw_yr, df$TLslope.z, pch=19, xlab="rho", ylab="z",
     col=df$color, ylim=c(0.1,3.5))
arrows(x0=df$avg_cor_btw_yr, y0=df$TLslope.z.lowCI, 
       x1=df$avg_cor_btw_yr, y1=df$TLslope.z.upCI, 
       code=3, angle=90, length=0.05, col=df$color)
abline(h=2, col="gray", lty=3)
par(op)
