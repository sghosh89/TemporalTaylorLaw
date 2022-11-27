rm(list=ls()) 
#------------
#library(remotes)
#install_github("seananderson/ecofolio")
#-------------
set.seed(seed=123)
library(ecofolio)
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
  
  spmat<-eachmat
  return(list(spmat=spmat,
              icv=icv,
               pe.avg.cv=pe.avg.cv,
               pe.mv=pe.mv))
}

# now generate data for CV vs. nsp for different z
nyr<-100
z_all<-seq(from=1,to=3,by=0.1)
nsplist<-c(30, 40, 50, 60, 70) # vary nsp from 30, 40, 50
res<-data.frame(z=NA*numeric(length(z_all)*length(nsplist)),
                nsp=NA*numeric(length(z_all)*length(nsplist)),
                icv=NA*numeric(length(z_all)*length(nsplist)),
                pe.avg.cv=NA*numeric(length(z_all)*length(nsplist)),
                pe.mv=NA*numeric(length(z_all)*length(nsplist)))
res$z<-rep(z_all,length(nsplist))
res$nsp<-rep(nsplist,each=length(z_all))

for(i in 1:nrow(res)){
  
  z<-res$z[i]
  nsp<-res$nsp[i]
  mu<-c(1:nsp)/10
  
  tempo<-gen_com_stab(nyr=nyr,nsp=nsp,mu=mu,z=z)
  if(nsp==70 & z==1){
    spmat_z1<-tempo$spmat 
  }else if(nsp==70 & z==2){
    spmat_z2<-tempo$spmat 
  }else if(nsp==70 & z==3){
    spmat_z3<-tempo$spmat 
  }
  res$icv[i]<-tempo$icv
  res$pe.avg.cv[i]<-tempo$pe.avg.cv
  res$pe.mv[i]<-tempo$pe.mv
}
#res$nsp<-as.factor(res$nsp)

# plot stability vs z for different richness
g1<-ggplot(data=res, aes(x=z,y=icv,col=as.factor(nsp)))+geom_point()+geom_smooth(se=F)+
  xlab("Taylor's law slope, z")+ylab("Stability, iCV")+
  theme_bw()+
  scale_colour_brewer("Richness",palette = "Blues")+
  #scale_color_manual("Richness",values=c("lightpink1","magenta","purple"))+
  #scale_fill_manual("Richness",values=c("lightpink1","magenta","purple"))+
  theme(text = element_text(size = 12),axis.text = element_text(size = 12),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+ 
  theme(legend.position = "right")
#g1

# plot stability vs z for different richness
g2<-ggplot(data=res, aes(x=nsp,y=icv,col=z))+geom_point()+
  xlab("Richness")+ylab("Stability, iCV")+
  theme_bw()+
  scale_fill_brewer("z",palette = "Blues")+
  theme(text = element_text(size = 12),axis.text = element_text(size = 12),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  theme(legend.position = "top")

# plot PE vs z for different richness
gPE<-ggplot(data=res, aes(x=z,y=pe.avg.cv,col=as.factor(nsp)))+geom_smooth(se=F)+
  xlab("Taylor's law slope, z")+ylab("Portfolio effect, PE")+
 geom_vline(xintercept = 2, linetype="dotted", col="gray")+
  theme_bw()+
  geom_smooth(aes(y=pe.mv),linetype = "dashed",se=F)+
  scale_colour_brewer("Richness",palette = "Blues")+
  theme(text = element_text(size = 12),axis.text = element_text(size = 12),
        legend.position="none",
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
annotate(geom="text",x=1.5, y=40, label="Solid: PE for avg Cv based,\n
         Dashed: PE for mean-variance scaling", size=3)
#gPE

# plot community abundance timeseries for z=1
spind<-c(10,20,30)
mus<-spind/10

a<-spmat_z1[,spind]
sa <- stack(as.data.frame(a))
sa$x <- rep(seq_len(nrow(a)), ncol(a))
dfa<-apply(a,MARGIN=1, FUN=sum)
dfa<-data.frame(x=c(1:nyr),tot=dfa, ind="sptot")
totmz1<-log(mean(dfa$tot))
totvz1<-log(var(dfa$tot))

gz1<-ggplot(data=sa,aes(x=x,y=values,col=as.factor(ind)))+geom_hline(yintercept = mus, linetype="dotted",col="gray")+
  geom_line()+
  theme_bw()+ylim(c(0.5,7))+xlab("Time in years")+ylab("Abundance, scaled")+
  scale_color_manual("Richness", values=c("palegreen","limegreen","green4","blue"))+
  theme(text = element_text(size = 12),axis.text = element_text(size = 12),
        legend.position="none",
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  annotate(geom="text",x=25, y=3.8, label="Three individual sp", size=5)+
  annotate(geom="text",x=50, y=4.8, label="z=1", size=5)+
geom_line(data=dfa,aes(x=x,y=tot),size=1.2,col="green")+
  annotate(geom="text",x=30, y=6.7, label="Community with 70 sp", size=5)


b<-spmat_z3[,spind]
sb <- stack(as.data.frame(b))
sb$x <- rep(seq_len(nrow(b)), ncol(b))
dfb<-apply(b,MARGIN=1, FUN=sum)
dfb<-data.frame(x=c(1:nyr),tot=dfb, ind="sptot")
totmz3<-log(mean(dfb$tot))
totvz3<-log(var(dfb$tot))

gz3<-ggplot(data=sb,aes(x=x,y=values,col=ind))+geom_hline(yintercept = mus, linetype="dotted")+
  geom_line()+
  theme_bw()+ylim(c(0.5,8))+xlab("Time in years")+ylab("Abundance, scaled")+
  scale_color_manual("Richness", values=c("gold","orange","hotpink"))+
  theme(text = element_text(size = 12),axis.text = element_text(size = 12),
        legend.position="none",
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  annotate(geom="text",x=25, y=4, label="Three individual sp", size=5)+
  annotate(geom="text",x=50, y=4.8, label="z=3", size=5)+
  geom_line(data=dfb,aes(x=x,y=tot),size=1.2,col="magenta")+
  annotate(geom="text",x=30, y=7.7, label="Community with 70 sp", size=5)

x<-spmat_z1
m <- apply(x, 2, mean)
v <- apply(x, 2, var)
log.m <- log(m)
log.v <- log(v)
dfz1<-data.frame(m=log.m,v=log.v)
highlight_dfz1 <- dfz1[spind,]

x<-spmat_z3
m <- apply(x, 2, mean)
v <- apply(x, 2, var)
log.m <- log(m)
log.v <- log(v)
dfz3<-data.frame(m=log.m,v=log.v)
highlight_dfz3 <- dfz3[spind,]


gTL<-ggplot(data=dfz1,aes(x=m,y=v))+geom_point(col="gray",pch=1,size=2.8)+
  geom_smooth(method="lm",se=F,col="green",linetype="dashed")+
  theme_bw()+xlab("log(mean)")+ylab("log(variance)")+
  theme(text = element_text(size = 12),axis.text = element_text(size = 12),
        legend.position="none",
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+
  geom_point(data=highlight_dfz1, 
             aes(x=m,y=v), 
             color=c("palegreen","limegreen","green4"),
             size=3)+
  geom_point(data=dfz3, 
             aes(x=m,y=v), 
             color="black",pch=1,
             size=2.8)+
  geom_smooth(data=dfz3,method="lm",se=F,col="magenta",linetype="dashed")+
  geom_point(data=highlight_dfz3, 
             aes(x=m,y=v), 
             color=c("gold","orange","hotpink"),
             size=3)+
  annotate(geom="text",x=-1.5, y=-4.5, label="slope, z=1", size=5)+
  annotate(geom="text",x=-1, y=-10, label="slope, z=3", size=5)+ 
  annotate("text", x = totmz1, y = totvz1, 
           parse = T, label = "X", col="green")+
  annotate("text", x = totmz3, y = totvz3, 
           parse = T, label = "X", col="magenta")

#gTL

pdf("../../Results/Prelim_res_plot/conceptual_fig1.pdf", width = 12, height = 6)
grid.arrange(gz1, gTL, gz3, g1, gPE,nrow = 2, 
             layout_matrix= rbind(c(1,2,3),c(4,NA,4.5)))
dev.off()



