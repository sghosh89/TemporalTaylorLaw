rm(list=ls()) 
library(tidyverse)
library(ggpubr)
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

# richness distribution


dens <- density(sm_all$nsp)
df <- data.frame(x=dens$x, y=dens$y)
q<-quantile(sm_all$nsp,c(0.25,0.75))
q
range(sm_all$nsp)
df$quant <- factor(findInterval(df$x,q))
gR<-ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + 
  scale_x_continuous(breaks=q) + scale_fill_manual(values=c("cyan","lightsteelblue1","dodgerblue"))+ 
  xlab("Richness, R")+ 
  ylab("Density")+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   legend.position="none",
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  annotate(geom="label",x=0, y=0.04, label="R = [15,34]", size=3, vjust=1, hjust=0, fill="cyan")+
  annotate(geom="label",x=30, y=0.04, label="R = [35,53]", size=3, vjust=1, hjust=0,fill="lightsteelblue1")+
annotate(geom="label",x=60, y=0.04, label="R = [54,89]", size=3, vjust=1, hjust=0,fill="dodgerblue")+
  annotate(geom="text",x=40, y=0.03, label="n = 1754 communities", size=4, vjust=1, hjust=0)

gRho<-ggplot() + 
    geom_density(data=sm_all, aes(x=avg_cor_btw_yr),fill="orange", alpha=0.4, adjust=1) + 
    xlab("Average correlation between years")+
  ylab("Density")+
    theme_bw()+theme(text = element_text(size = 12),axis.text = element_text(size = 12),
                     plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                     panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))
  

mean(sm_all$nsp) #~44 sp
mean(sm_all$nyr)# ~23 yr

# considering by species
ggplot(sm_all, aes(x = net_taildep/nint, y = TLslope.z, 
               color=as.factor(nsp), group=as.factor(nsp))) + 
  geom_point(alpha=0.3) + 
  geom_smooth(method="lm",se=F)+facet_wrap(~nsp)+
  theme_classic()+theme(legend.position = "none")

# considering by groups: 2 REALM
ggplot(sm_all, aes(x = net_taildep/nint, y = TLslope.z,
               color=REALM, group=REALM)) + 
  geom_point(alpha=0.3,col=as.factor(sm_all$nsp)) + 
  geom_smooth(method="lm",se=T, aes(fill=REALM))+
  scale_color_manual(values=c("skyblue","green"))+
  scale_fill_manual(values=c("skyblue","darkolivegreen2"))+
  theme_classic()

# considering all point
gTA<-ggplot(sm_all, aes(x = net_taildep/nint, y = TLslope.z)) + ylab("Taylor's slope, z")+
  geom_point(alpha=0.3,col="black",pch=21,fill="black") +
  geom_smooth(method="loess",se=T,col="red",aes(fill="red"))+
  theme_bw()+theme(text = element_text(size = 12),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),legend.position="none",
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))


# make plot for z vs. net_taildep/nint for different range of species
q<-quantile(sm_all$nsp,c(0.25,0.75))
sm_all<-sm_all %>%
  mutate(richness_cat = case_when(nsp < q[1] ~ "low R, <50%CI",
                                           nsp > q[2] ~ "high R, >50%CI",
                                           nsp <= q[2] & nsp >= q[1] ~ "mid R, in 50%CI"
                                            )) %>%
  mutate(richness_col = case_when(nsp < q[1] ~ "cyan",
                                  nsp > q[2] ~ "lightsteelblue1",
                                  nsp <= q[2] & nsp >= q[1] ~ "dodgerblue"))

gz_TA_R<-ggplot(sm_all, aes(x = net_taildep/nint, y = TLslope.z)) + 
  geom_hline(yintercept=2,linetype="dotted",col="black")+
  ylab("Taylor's slope, z")+xlab("Dependence in ranks (upper to lower)")+
  geom_point(alpha=0.8,fill=sm_all$richness_col,pch=21) +
  geom_smooth(method="lm",se=T,col="red",fill="red")+
  theme_bw()+theme(text = element_text(size = 12),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  stat_cor(method = "pearson") 


# plot for variance ratio vs z
gz_VR<-ggplot(sm_all, aes(x = vr_LdM, y = TLslope.z)) + 
  geom_hline(yintercept=2,linetype="dotted",col="black")+
  geom_vline(xintercept=0.5,linetype="dotted",col="black")+
  ylab("Taylor's slope, z")+xlab("Synchrony, VR")+xlim(c(0,1))+
  geom_point(alpha=0.8,fill=sm_all$richness_col,pch=21) +
  theme_bw()+theme(text = element_text(size = 12),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))

# plot for avg. cor btw years vs z
gz_rho<-ggplot(sm_all, aes(x = avg_cor_btw_yr, y = TLslope.z)) + 
  geom_hline(yintercept=2,linetype="dotted",col="black")+
  #geom_vline(xintercept=0.5,linetype="dotted",col="black")+
  ylab("Taylor's slope, z")+xlab("Average correlation between years")+
  geom_point(alpha=0.8,fill=sm_all$richness_col,pch=21) +
  theme_bw()+theme(text = element_text(size = 12),axis.text = element_text(size = 12),
                   plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
                   panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))



gR<-gR+annotate(geom="text",x=90, y=0.04, label="A", size=7, vjust=1, hjust=0,fill="dodgerblue")
gz_rho<-gz_rho+annotate(geom="text",x=0.85, y=4, label="B", size=7, vjust=1, hjust=0,fill="dodgerblue")
gz_VR<-gz_VR+annotate(geom="text",x=0.95, y=3.5, label="C", size=7, vjust=1, hjust=0,fill="dodgerblue")
gz_TA_R<-gz_TA_R+annotate(geom="text",x=0.01, y=3.5, label="D", size=7, vjust=1, hjust=0,fill="dodgerblue")

pdf("../../Results/Prelim_res_plot/test_empirical_z_taildep.pdf", width = 7, height = 7)
grid.arrange(gR,gz_rho,
             gz_VR,gz_TA_R,nrow = 2)
dev.off()


#hist(sm_all$TLslope.z, breaks=100) # nearly normal distribn
#hist(sm_all$net_taildep/sm_all$nint, breaks=100) # normal distribn









