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
library(mvtnorm)
`%notin%` <- Negate(`%in%`)
#================== with correlation structure =================
# first, check
#nyr<-100
#a<-0.01
#nsp<-70
#z<-2.5
#rho_mean<- 0.5

#mu<-c(1:nsp)/10
#mu<- sample(x=c(20:100), size=nsp, replace = F)

#range(eachmat)
#=========================================

gen_com_stab_w_interaction<-function(nyr=100, nsp=70, a=0.01, mu, z, rho_mean){
  
  #=======================
  lenmu<-length(mu)
  sigma2 <- matrix(lenmu*NA, lenmu, lenmu)
  diag(sigma2)<-(a*mu^z) 
  idmat<-which(upper.tri(sigma2)==1,arr.ind = T)
  idmat<-as.data.frame(idmat)
  
  rho_ij<-rnorm(n=(lenmu^2-lenmu)*0.5,mean=rho_mean,sd=0.001)
  #hist(rho_ij)
  
  for(i in 1:nrow(idmat)){
    ri<-idmat$row[i]
    rj<-idmat$col[i]
    #sigma[ri,rj]<-a*rho_ij[i]*(mu[ri]*mu[rj])^(z/2)
    sigma_i2<-a*(mu[ri])^z
    sigma_j2<-a*(mu[rj])^z
    sigma2[ri,rj]<-rho_ij[i]*sqrt(sigma_i2*sigma_j2)
  }
  sigma2[lower.tri(sigma2)]<-t(sigma2)[lower.tri(sigma2)]
  
  eachmat<-rmvnorm(n=nyr,mean=mu,sigma=sigma2)# sp. timeseries
  #========================
  
  pe.avg.cv<-pe_avg_cv(x=eachmat, ci = F, na.rm = F) # default not detrended
  pe.mv<-pe_mv(x=eachmat, ci = F, na.rm = F) # default not detrended
  
  tot_ts<-apply(eachmat,MARGIN = 1, FUN=sum)
  res<-ecofolio::cv(tot_ts)
  icv<-1/res
  
  zfit<-fit_taylor(x=eachmat)
  
  m1<-rbind(eachmat,apply(eachmat,MARGIN=2,FUN=sum))
  m1<-tail(m1,1)
  m1<-as.data.frame(m1)
  #evenness<-adiv::specieseve(eachmat, method = "SmithWilson", tol = 1e-8)
  #SmithWilson<-evenness[1,1]
  lms<-ecofolio::synchrony(eachmat)# Loreau-Mazancourt synchrony
  
  spmat<-eachmat
  return(list(spmat=spmat,
              icv=icv,
              fittaylor_c=zfit$c,
              fittaylor_z=zfit$z,
              pe.avg.cv=pe.avg.cv,
              pe.mv=pe.mv,
              #SmithWilson=SmithWilson,
              lms=lms))
}

# now generate data for CV vs. nsp for different z
nrep<-50 # 
nyr<-100
z_all<-seq(from=1.5,to=2.5,by=0.1)
nsplist<-c(30,50,70) # vary nsp 
seedlist<-c(123:(nrep+122))
#nsplist<-c(70)
rholist<-seq(from=-0.8,to=0.8,by=0.4)
lendf<-length(z_all)*length(nsplist)*length(rholist)*length(seedlist)
res<-data.frame(z=NA*numeric(lendf),
                nsp=NA*numeric(lendf),
                icv=NA*numeric(lendf),
                zfit=NA*numeric(lendf),
                pe.avg.cv=NA*numeric(lendf),
                pe.mv=NA*numeric(lendf),
                #SmithWilson=NA*numeric(lendf),
                lms=NA*numeric(lendf),
                nrep=NA*numeric(lendf))
res$z<-rep(z_all,length(nsplist))
res$nsp<-rep(nsplist,each=length(z_all))
res$rho<-rep(rholist,each=length(z_all))
res$nrep<-rep(seedlist,each=length(z_all))

badid<-c()
for(i in 1:nrow(res)){
  
  set.seed(seed=res$nrep[i])
  z<-res$z[i]
  nsp<-res$nsp[i]
  mu<-runif(nsp,min=1,max=nsp)/10
  #mu<-c(1:nsp)/10
  rho_mean<-res$rho[i]
    
  tempo<-gen_com_stab_w_interaction(nyr=nyr, nsp=nsp, a=0.01, 
                                    mu, z=z, rho_mean=rho_mean)
  
  res$icv[i]<-tempo$icv
  res$zfit[i]<-tempo$fittaylor_z
  res$pe.avg.cv[i]<-tempo$pe.avg.cv
  res$pe.mv[i]<-tempo$pe.mv
  #res$SmithWilson[i]<-tempo$SmithWilson
  res$lms[i]<-tempo$lms
  
  if(min(tempo$spmat)<0){
    badid<-c(badid,i)
    print(i)
  }
  
}
res$i<-1:nrow(res)
res<-rename(res, Richness=nsp)
resg<-res%>%filter(i%notin%badid)
#======================

saveRDS(resg, "../../Results/Prelim_res_plot/data_for_conceptual_fig2_with_interaction.RDS")
#dd<-resg%>%group_by(Richness,z,rho)%>%summarise(mnicv=mean(icv))
#ggplot(data=dd, aes(x=z,y=mnicv,col=as.factor(Richness)))+
#  geom_point()+
#  geom_line()+
#  facet_wrap(~rho, scales="free", labeller = label_both)

g1<-ggplot(data=resg, aes(x=z,y=icv,col=as.factor(Richness)))+
  #geom_point()+
  geom_smooth(se=F)+
  xlab("Taylor's law slope, z")+ylab("Stability (=1/CV)")+
  theme_bw()+
  scale_colour_brewer("Richness",palette = "Set2")+
  #scale_color_manual("Richness",values=c("#66C2A5","#8DA0CB"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+ 
  theme(legend.position ="right")+
  facet_wrap(~rho, scales="free", labeller = label_both, ncol=1)
g1

g1inset<-ggplot(data=resg, aes(x=Richness,y=icv,col=as.factor(z)))+
  #geom_point()+
  geom_smooth(se=F)+
  xlab("Richness, R")+ylab("Stability (=1/CV)")+
  theme_bw()+
  scale_colour_brewer("z",palette = "RdYlBu", direction=-1)+
  #scale_color_manual("Richness",values=c("#66C2A5","#8DA0CB"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+ 
  theme(legend.position ="right")+
  facet_wrap(~rho, scales="free", labeller = label_both, ncol=1)
g1inset

g2<-ggplot(data=resg, aes(x=z,y=pe.avg.cv,col=as.factor(Richness)))+
  #geom_point()+geom_line()+
  geom_smooth(se=F)+
  xlab("Taylor's law slope, z")+ylab("Portfolio effect, PE")+
  geom_vline(xintercept = 2, linetype="dotted", col="gray40")+
  theme_bw()+
  #geom_point(aes(y=pe.mv), pch=1)+geom_line(aes(y=pe.mv),linetype = "dashed")+
  geom_smooth(aes(y=pe.mv),linetype = "dashed",se=F)+
  scale_colour_brewer("Richness",palette = "Set2")+
  #scale_color_manual("Richness",values=c("#66C2A5","#8DA0CB"))+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14),
        plot.margin = margin(t = 8, r = 9, b = 4, l = 4, unit = "pt"),
        panel.grid = element_line(color = rgb(235, 235, 235, 100, maxColorValue = 255)))+ 
  theme(legend.position =c(0.8,0.2))+
  #facet_grid(Richness~rho, labeller = label_both)
  facet_wrap(~rho, scales="free", labeller = label_both)
g2


pdf("../../Results/Prelim_res_plot/conceptual_fig2stability_with_interaction.pdf", width = 4, height = 12)
g1
dev.off()

pdf("../../Results/Prelim_res_plot/conceptual_fig2stability_with_interaction_inset.pdf", width = 4, height = 12)
g1inset
dev.off()

pdf("../../Results/Prelim_res_plot/conceptual_fig2PE_with_interaction.pdf", width = 10, height = 6)
g2
dev.off()


