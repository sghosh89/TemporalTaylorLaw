rm(list=ls()) 
#------------
#library(remotes)
#install_github("seananderson/ecofolio")
#-------------
set.seed(seed=123)
library(ecofolio)
library(adiv)
library(ggplot2)
library(gridExtra)
library(dplyr)
library("RColorBrewer") 

#z in the range of 0.6 to 2.8
calcsd2<- function (mu,z) {sigmasq<- 0.01*mu^z} # variance ~ constant * mean^z

gen_com_stab<-function(nyr, nsp, mu, z){
  
  eachmat<-matrix(NA, nrow=nyr, ncol=nsp)
  rownames(eachmat)<-paste("year",1:nyr,sep="")
  colnames(eachmat)<-paste("sp",1:nsp,sep="")
  
  for(i in 1:nsp){
    mui<-mu[i]
    eachmat[,i]<- rnorm(n=nyr, mean=mui, sd=sqrt(calcsd2(mui,z)))
  }
  pe.avg.cv<-pe_avg_cv(x=eachmat, ci = F, na.rm = F) # default not detrended
  pe.mv<-pe_mv(x=eachmat, ci = F, na.rm = F) # default not detrended
  
  tot_ts<-apply(eachmat,MARGIN = 1, FUN=sum)
  res<-ecofolio::cv(tot_ts)
  icv<-1/res
  
  zfit<-fit_taylor(x=eachmat)
  
  m1<-rbind(eachmat,apply(eachmat,MARGIN=2,FUN=sum))
  m1<-tail(m1,1)
  m1<-as.data.frame(m1)
  evenness<-adiv::specieseve(eachmat, method = "SmithWilson", tol = 1e-8)
  SmithWilson<-evenness[1,1]
  
  spmat<-eachmat
  return(list(spmat=spmat,
              icv=icv,
              fittaylor_c=zfit$c,
              fittaylor_z=zfit$z,
               pe.avg.cv=pe.avg.cv,
               pe.mv=pe.mv,
              SmithWilson=SmithWilson))
}

# now generate data for CV vs. nsp for different z
nyr<-100
z_all<-seq(from=1,to=3,by=0.1)
nsplist<-c(30, 50, 70) # vary nsp 
res<-data.frame(z=NA*numeric(length(z_all)*length(nsplist)),
                nsp=NA*numeric(length(z_all)*length(nsplist)),
                icv=NA*numeric(length(z_all)*length(nsplist)),
                zfit=NA*numeric(length(z_all)*length(nsplist)),
                pe.avg.cv=NA*numeric(length(z_all)*length(nsplist)),
                pe.mv=NA*numeric(length(z_all)*length(nsplist)),
                SmithWilson=NA*numeric(length(z_all)*length(nsplist)))
res$z<-rep(z_all,length(nsplist))
res$nsp<-rep(nsplist,each=length(z_all))

for(i in 1:nrow(res)){
  
  z<-res$z[i]
  nsp<-res$nsp[i]
  mu<-c(1:nsp)/10
  
  tempo<-gen_com_stab(nyr=nyr,nsp=nsp,mu=mu,z=z)
  
  res$icv[i]<-tempo$icv
  res$zfit[i]<-tempo$fittaylor_z
  res$pe.avg.cv[i]<-tempo$pe.avg.cv
  res$pe.mv[i]<-tempo$pe.mv
  res$SmithWilson[i]<-tempo$SmithWilson
}
#res$nsp<-as.factor(res$nsp)

# z and zfit same?
res2<-res%>%filter(nsp==70)
ggplot(data=res2,aes(x=z,y=zfit))+geom_point()+geom_abline(slope=1,intercept=0)


# plot stability vs z for different richness
g1<-ggplot(data=res, aes(x=z,y=icv,col=as.factor(nsp)))+geom_smooth(se=F)+
  xlab("Taylor's law slope, z")+ylab("Stability, iCV")+
  theme_bw()+
  scale_colour_brewer("Richness",palette = "Set2")+
  #scale_color_manual("Richness",values=c("lightpink1","magenta","purple"))+
  #scale_fill_manual("Richness",values=c("lightpink1","magenta","purple"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+ 
  theme(legend.position = c(0.8,0.5))
g1

# plot stability vs richness for different z
rese<-res%>%filter(z%in%c(1,2,3))%>%arrange(z)
g2<-ggplot(data=rese, aes(x=nsp,y=icv,col=as.factor(z)))+geom_point()+geom_line()+
  xlab("Richness")+ylab("Stability, iCV")+
  theme_bw()+
  scale_color_manual("z", values=c("darkblue","goldenrod4","magenta"))+
  theme(text = element_text(size = 16),axis.text = element_text(size = 16),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  theme(legend.position = c(0.2,0.8))

# evenness: no change with z
g2e<-ggplot(data=res, aes(x=SmithWilson,y=icv,col=z))+geom_point()+
  xlab("Evenness")+ylab("Stability, iCV")+
  theme_bw()+
  scale_fill_brewer("z",palette = "Blues")+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  theme(legend.position = "top")


# plot PE vs z for different richness
gPE<-ggplot(data=res, aes(x=z,y=pe.avg.cv,col=as.factor(nsp)))+geom_smooth(se=F)+
  xlab("Taylor's law slope, z")+ylab("Portfolio effect, PE")+
 geom_vline(xintercept = 2, linetype="dotted", col="gray")+
  theme_bw()+
  geom_smooth(aes(y=pe.mv),linetype = "dashed",se=F)+
  scale_colour_brewer("Richness",palette = "Set2")+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        legend.position="none",
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
annotate(geom="text",x=1.8, y=40, 
         label="Solid: PE without mean-variance scaling,\n Dashed: PE with mean-variance scaling", size=4.5)+
annotate(geom="text",x=1.5, y=15, label="Overestimate PE", size=5)+
  annotate(geom="text",x=2.5, y=0, label="Underestimate PE", size=5)
gPE

# Figure 2
pdf("../../Results/Prelim_res_plot/conceptual_fig2.pdf", width = 10, height = 5)
g1<-g1+annotate(geom="text",x=Inf, y=Inf,hjust=1.5,vjust=1.2, label="A", size=8)
gPE<-gPE+annotate(geom="text",x=-Inf, y=Inf,hjust=-0.5,vjust=1.2, label="B", size=8)
grid.arrange(g1, gPE, nrow = 1)
dev.off()

pdf("../../Results/Prelim_res_plot/conceptual_fig2inset.pdf", width = 3.8, height = 3.8)
grid.arrange(g2, nrow = 1)
dev.off()


#==========================================

# now generate data for CV vs. nsp for different z
nyr<-100
z_all<-seq(from=1,to=3,by=0.1)
nsp<-70

df_R70_z123<-data.frame(z=z_all,
                        zfit=NA*numeric(length(z_all)),
                        iCV=NA*numeric(length(z_all)),
                        delta=NA*numeric(length(z_all))) 

for(i in 1:nrow(df_R70_z123)){
  
  z<-df_R70_z123$z[i]
  mu<-c(1:nsp)/10
  
  tempo<-gen_com_stab(nyr=nyr,nsp=nsp,mu=mu,z=z)
  
  df_R70_z123$iCV[i]<-tempo$icv
  df_R70_z123$zfit[i]<-tempo$fittaylor_z
  spmat<-tempo$spmat
  
  log.v.actual=log(var(apply(spmat, 1, sum)))
  totcommunity<-data.frame(log.m=log(mean(apply(spmat, 1, sum))))# mean of total community biomass
                      
  m <- apply(spmat, 2, mean)
  v <- apply(spmat, 2, var)
  log.m <- log(m)
  log.v <- log(v)
  fit <- lm(log.v ~ log.m)
  logv.predicted<-predict(fit,totcommunity) # predicted as per Taylor's fit
  
  df_R70_z123$delta[i]<-logv.predicted - log.v.actual
  
  if(z==1){
    spmat_z1<-spmat
  }
  
  if(z==2){
    spmat_z2<-spmat
  }
  
  if(z==3){
    spmat_z3<-spmat
  }
  
}


# plot community abundance timeseries for z=1
spind<-c(20,30,40)
mus<-spind/10

a<-spmat_z1[,spind]
sa <- stack(as.data.frame(a))
sa$x <- rep(seq_len(nrow(a)), ncol(a))

dfa<-apply(spmat_z1,MARGIN=1, FUN=sum)
dfa<-data.frame(x=c(1:nyr),tot=dfa, ind="sptot")
totmz1<-log(mean(dfa$tot))
totvz1<-log(var(dfa$tot))

scale=30
gz1<-ggplot(data=sa,aes(x=x,y=values,col=as.factor(ind)))+
  geom_hline(yintercept = c(mus,mean(dfa$tot)/scale), linetype="dotted",col="black")+
  geom_line()+
  theme_bw()+xlab("Time in years")+ylab("Abundance, scaled")+
  scale_color_manual("Richness", values=c("lightskyblue1","deepskyblue","blue"))+
  annotate(geom="text",x=45, y=9, label="Abundance timeseries, z=1", size=5)+
geom_line(data=dfa,aes(x=x,y=tot/scale),size=1.2,col="darkblue")+
  scale_y_continuous(
    # Features of the first axis
    name = "Three individual sp.",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*scale, name="Total community (R=70 sp)"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        legend.position="none",axis.title.y.right = element_text(color = "darkblue", face="bold"),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))

gz1

b<-spmat_z2[,spind]
sb <- stack(as.data.frame(b))
sb$x <- rep(seq_len(nrow(b)), ncol(b))

dfb<-apply(spmat_z2,MARGIN=1, FUN=sum)
dfb<-data.frame(x=c(1:nyr),tot=dfb, ind="sptot")
totmz2<-log(mean(dfb$tot))
totvz2<-log(var(dfb$tot))

gz2<-ggplot(data=sb,aes(x=x,y=values,col=as.factor(ind)))+
  geom_hline(yintercept = c(mus,mean(dfb$tot)/scale), linetype="dotted",col="black")+
    geom_line()+
  theme_bw()+xlab("Time in years")+ylab("Abundance, scaled")+
  scale_color_manual("Richness", values=c("goldenrod1","goldenrod2","goldenrod3"))+
  annotate(geom="text",x=45, y=9, label="Abundance timeseries, z=2", size=5)+
  geom_line(data=dfb,aes(x=x,y=tot/scale),size=1.2,col="goldenrod4")+
  #annotate(geom="text",x=25, y=5.5, label="z=2", size=5)+ 
  scale_y_continuous(
    # Features of the first axis
    name = "Three individual sp.",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*scale, name="Total community (R=70 sp)"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        legend.position="none",axis.title.y.right = element_text(color = "goldenrod4", face="bold"),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))
gz2

b<-spmat_z3[,spind]
sb <- stack(as.data.frame(b))
sb$x <- rep(seq_len(nrow(b)), ncol(b))

dfb<-apply(spmat_z3,MARGIN=1, FUN=sum)
dfb<-data.frame(x=c(1:nyr),tot=dfb, ind="sptot")
totmz3<-log(mean(dfb$tot))
totvz3<-log(var(dfb$tot))

gz3<-ggplot(data=sb,aes(x=x,y=values,col=as.factor(ind)))+
  geom_hline(yintercept = c(mus,mean(dfb$tot)/scale), linetype="dotted",col="black")+
  geom_line()+
  theme_bw()+xlab("Time in years")+ylab("Abundance, scaled")+
  scale_color_manual("Richness", values=c("orchid1","orchid2","orchid3"))+
  annotate(geom="text",x=45, y=9, label="Abundance timeseries, z=3", size=5)+
  geom_line(data=dfb,aes(x=x,y=tot/scale),size=1.2,col="magenta")+
  #annotate(geom="text",x=10, y=5.5, label="z=3", size=5)+ 
  scale_y_continuous(
    # Features of the first axis
    name = "Three individual sp.",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*scale, name="Total community (R=70 sp)"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        legend.position="none",axis.title.y.right = element_text(color = "magenta", face="bold"),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))
gz3





x<-spmat_z1
m <- apply(x, 2, mean)
v <- apply(x, 2, var)
log.m <- log(m)
log.v <- log(v)
dfz1<-data.frame(m=log.m,v=log.v)
highlight_dfz1 <- dfz1[spind,]

x<-spmat_z2
m <- apply(x, 2, mean)
v <- apply(x, 2, var)
log.m <- log(m)
log.v <- log(v)
dfz2<-data.frame(m=log.m,v=log.v)
highlight_dfz2 <- dfz2[spind,]

x<-spmat_z3
m <- apply(x, 2, mean)
v <- apply(x, 2, var)
log.m <- log(m)
log.v <- log(v)
dfz3<-data.frame(m=log.m,v=log.v)
highlight_dfz3 <- dfz3[spind,]

s1<-lm(v~m,data=dfz1)
s2<-lm(v~m,data=dfz2)
s3<-lm(v~m,data=dfz3)


gTL<-ggplot(data=dfz1,aes(x=m,y=v))+geom_point(col="lightskyblue1",pch=1,size=2.8)+
  geom_abline(slope = unname(s1$coefficients[2]), 
              intercept = unname(s1$coefficients[1]),col="darkblue",linetype="dotted")+
  geom_abline(slope = unname(s2$coefficients[2]), 
              intercept = unname(s2$coefficients[1]),col="goldenrod4",linetype="dotted")+
  geom_abline(slope = unname(s3$coefficients[2]), 
              intercept = unname(s3$coefficients[1]),col="magenta",linetype="dotted")+
  geom_smooth(method="lm",se=F,col="darkblue",linetype=1)+ylim(-12,12)+
  theme_bw()+xlab("log(mean)")+ylab("log(variance)")+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        legend.position="none",
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  geom_point(data=highlight_dfz1, 
             aes(x=m,y=v), 
             fill=c("lightskyblue1","deepskyblue","blue"),
             size=3, pch=21)+
  geom_point(data=dfz2, 
             aes(x=m,y=v), 
             color="gold",pch=1,
             size=2.8)+geom_smooth(data=dfz2,method="lm",se=F,col="wheat4",linetype=1)+
  geom_point(data=highlight_dfz2, 
             aes(x=m,y=v), 
             fill=c("goldenrod1","goldenrod2","goldenrod3"),
             size=3,pch=21)+
  geom_point(data=dfz3, 
             aes(x=m,y=v), 
             color="plum",pch=1,
             size=2.8)+
  geom_smooth(data=dfz3,method="lm",se=F,col="magenta",linetype=1)+
  geom_point(data=highlight_dfz3, 
             aes(x=m,y=v), 
             fill=c("orchid1","orchid2","orchid3"),
             size=3,pch=21)+
  annotate("text", x = totmz1, y = totvz1, 
           parse = T, label = "X", col="darkblue")+
  annotate("text", x = totmz3, y = totvz3, 
           parse = T, label = "X", col="magenta")+
  annotate("text", x = totmz2, y = totvz2, 
           parse = T, label = "X", col="goldenrod4")

gTL

scale2=10
gnew<-ggplot(df_R70_z123,aes(x=zfit,y=iCV))+geom_point(col="black")+geom_smooth(se=F,col="black")+
  geom_point(aes(y=delta*scale2),col="red")+geom_smooth(se=F,aes(y=delta*10),col="red")+ 
  scale_y_continuous(
    # Features of the first axis
    name = "Stability, iCV",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*(1/scale2), name="Delta"))+theme_bw()+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        legend.position="none",axis.title.y.right = element_text(color = "red"),
        axis.title.y= element_text(color = "black"),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  xlab("Taylor's slope, z")


pdf("../../Results/Prelim_res_plot/conceptual_fig1.pdf", width = 8, height =10)
grid.arrange(gz1, gz2, gz3, gTL, gnew, ncol = 2, 
             layout_matrix= cbind(c(1,2,3),c(4,5,NA)))
dev.off()





