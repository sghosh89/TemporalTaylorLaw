# temporal Taylor's slope
#----------------------------------
rm(list=ls()) 
source("copsurrognd.R")
library(mvtnorm)
library(ecofolio)
library(copula)
library(tidyverse)
set.seed(seed=123)

#----------- secondary functions ------
cop_to_marg_trans<-function(x){
  return(qgamma(x,shape=gshape,scale=gscale)) #this is where we choose the marginals,
  #should be something skewed
}

marg_to_cop_trans<-function(x){ #transforms to copula space
  return(pgamma(x,shape=gshape,scale=gscale))
}
#--------------------------------------------

#***constants
N<-22 #length of time series
n<-40 #number of sp = R
numdat<-1000 #the number of surrogates

rlist<-c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) #normal-copula parameter for dependence between years

m_list <- vector(mode = "list", length = length(rlist)) # holds each array in each space
names(m_list)<-paste("r_",rlist,sep="") # for Normal copula
C_list<-SC_list<-m_list # for Clayton and Survival Clayton
for(k in 1:length(rlist)){
  r<-rlist[k]
  #***generate the starting data, which has normal copula structure
  #because it comes originally from a multivariate normal, before transformation
  sig<-matrix(r,N,N)
  diag(sig)<-1
  m<-rmvnorm(n*numdat,mean=rep(0,N),sigma=sig) # each column of m has normal histogram
  
  #dim(m)
  m<-array(m,c(n,numdat,N))
  #dim(m)
  m<-aperm(m,c(1,3,2)) 
  #dim(m)# 2nd dim is normally distributed 
  #e.g. hist(m[,2,],50)
  #plot(m[,1,1],m[,2,1]) # yr1 vs yr2
  
  gshape<-7.5
  gscale<-1 #gamma parameters
  
  # pnorm(m) makes normal to uniform
  # then qgamma(...) makes uniform to gamma distributed
  mg<-cop_to_marg_trans(pnorm(m))
  #dim(mg)
  #plot(mg[,1,1],mg[,2,1]) #yr1 timeseries vs. yr2 timeseries
  #hist(mg[,1,]) #gamma
  
  #***surrogates using a Clayton copula, spearman preserving
  #get the Clayton parameter to use
  ncop<-normalCopula(r,2)
  
  ccop<-claytonCopula(2,2)
  cparam<-iRho(ccop,rho(ncop))
  ccop<-claytonCopula(cparam,N)
  
  #generate the surrogates
  csurrog<-array(NA,dim(mg))
  #nsurrog<-csurrog
  for (counter in 1:dim(mg)[3]){
    csurrog[,,counter]<-copsurrognd(mg[,,counter],targetcop=ccop,numsurrog=1)
    #nsurrog[,,counter]<-copsurrognd(mg[,,counter],targetcop=ncop,numsurrog=1)
  }
  
  #***surrogates using a survival Clayton, spearman preserving
  #take the clayton surrogates and transform them to the 
  #copula space, then flip, then transform back
  scsurrog<-cop_to_marg_trans((1-marg_to_cop_trans(csurrog)))
  
  # make N by n by numdat = year by species by replica
  m1<-aperm(mg,c(2,1,3)) # normal copula
  #dim(m1)
  c1<-aperm(csurrog,c(2,1,3))
  sc1<-aperm(scsurrog,c(2,1,3))
  m_list[[k]]<-m1
  C_list[[k]]<-c1
  SC_list[[k]]<-sc1
}

res<-list(Normalsurrogs=m_list,
          Csurrogs=C_list,
          SCsurrogs=SC_list)
saveRDS(res,paste("../../Results/gather_res/surrogates_nsp_",n,"_nyr_",N,".RDS",sep=""))

z_N_r<-matrix(NA,nrow=numdat,ncol=length(rlist))
colnames(z_N_r)<-paste("r_",rlist,sep="")
rownames(z_N_r)<-paste("surrog_",1:numdat,sep="")
z_C_r<-z_SC_r<-z_N_r
for(i in 1:length(rlist)){
  for(j in 1:numdat){
    res<-fit_taylor(m_list[[i]][,,j])
    z_N_r[j,i]<-res$z
    res2<-fit_taylor(C_list[[i]][,,j])
    z_C_r[j,i]<-res2$z
    res3<-fit_taylor(SC_list[[i]][,,j])
    z_SC_r[j,i]<-res3$z
  }
}

zs_N<-apply(z_N_r,FUN=mean,MARGIN=2)
zs_C<-apply(z_C_r,FUN=mean,MARGIN=2)
zs_SC<-apply(z_SC_r,FUN=mean,MARGIN=2)

sd_zs_N<-apply(z_N_r,FUN=sd,MARGIN=2)
sd_zs_C<-apply(z_C_r,FUN=sd,MARGIN=2)
sd_zs_SC<-apply(z_SC_r,FUN=sd,MARGIN=2)

# z vs. increasing synchrony between species 
pdf(paste("../../Results/gather_res/z_vs_btwyr_cor_simulation_nsp_",n,"_nyr_",N,".pdf",sep=""),
    height=6,width=6)
plot(rlist,zs_N,xlab="r (Between year correlation)",ylab="z (N,C,SC)",type="b",lty=3,
     ylim=c(-1,4),xlim=c(0,1),pch=16)
arrows(x0=rlist, y0=zs_N-sd_zs_N, x1=rlist, y1=zs_N+sd_zs_N, 
       code=3, angle=90, length=0.05)
lines(rlist,zs_C,pch=16,col="red",type="b",lty=3)
arrows(x0=rlist, y0=zs_C-sd_zs_C, x1=rlist, y1=zs_C+sd_zs_C, 
       code=3, angle=90, length=0.05,col="red")
lines(rlist,zs_SC,pch=16,col="blue",type="b",lty=3)
arrows(x0=rlist, y0=zs_SC-sd_zs_SC, x1=rlist, y1=zs_SC+sd_zs_SC, 
       code=3, angle=90, length=0.05,col="blue")
dev.off()
# just check
#plot(C_list[[5]][1,,1],C_list[[5]][2,,1]) # lower tail dep.
#cor(C_list[[5]][1,,1],C_list[[5]][2,,1],method="spearman")

# function to compute average positive and negative 
# pairwise Spearman correlations between columns (species)
source("compute_avg_cor.R")

avg_cor_posmat_N<-matrix(NA,nrow=numdat,ncol=length(rlist))
colnames(avg_cor_posmat_N)<-paste("r_",rlist,sep="")
rownames(avg_cor_posmat_N)<-paste("surrog_",1:numdat,sep="")

avg_cor_posmat_C<-avg_cor_posmat_SC<-avg_cor_posmat_N
avg_cor_mat_C<-avg_cor_mat_SC<-avg_cor_mat_N<-avg_cor_posmat_N
avg_cor_negmat_N<-avg_cor_negmat_C<-avg_cor_negmat_SC<-avg_cor_posmat_N
vr_LdM_N<-vr_LdM_C<-vr_LdM_SC<-avg_cor_posmat_N
# this is not needed, just to check
avg_cor_btwyr_N<-avg_cor_btwyr_C<-avg_cor_btwyr_SC<-avg_cor_posmat_N

for(i in 1:length(rlist)){
  for(j in 1:numdat){
    x1<-m_list[[i]][,,j]
    x2<-C_list[[i]][,,j]
    x3<-SC_list[[i]][,,j]
    res1<-compute_avg_cor(mat=x1)
    res2<-compute_avg_cor(mat=x2)
    res3<-compute_avg_cor(mat=x3)
    
    vr_LdM_N[j,i]<-ecofolio::synchrony(x1)
    vr_LdM_C[j,i]<-ecofolio::synchrony(x2)
    vr_LdM_SC[j,i]<-ecofolio::synchrony(x3)
    
    avg_cor_mat_N[j,i]<-res1$avg_cor_btw_sp
    avg_cor_mat_C[j,i]<-res2$avg_cor_btw_sp
    avg_cor_mat_SC[j,i]<-res3$avg_cor_btw_sp
    
    avg_cor_posmat_N[j,i]<-res1$avg_cor_pos_btw_sp
    avg_cor_posmat_C[j,i]<-res2$avg_cor_pos_btw_sp
    avg_cor_posmat_SC[j,i]<-res3$avg_cor_pos_btw_sp
    
    avg_cor_negmat_N[j,i]<-res1$avg_cor_neg_btw_sp
    avg_cor_negmat_C[j,i]<-res2$avg_cor_neg_btw_sp
    avg_cor_negmat_SC[j,i]<-res3$avg_cor_neg_btw_sp
    
    avg_cor_btwyr_N[j,i]<-res1$avg_cor_btw_yr
    avg_cor_btwyr_C[j,i]<-res2$avg_cor_btw_yr
    avg_cor_btwyr_SC[j,i]<-res3$avg_cor_btw_yr
  }
}

#avg_cor_negmat_N[is.na(avg_cor_negmat_N)]<-0
#avg_cor_negmat_C[is.na(avg_cor_negmat_C)]<-0
#avg_cor_negmat_SC[is.na(avg_cor_negmat_SC)]<-0

avg_cor_matlist<-list(avg_cor_mat_N=avg_cor_mat_N,
                      avg_cor_mat_C=avg_cor_mat_C,
                      avg_cor_mat_SC=avg_cor_mat_SC,
                      avg_cor_posmat_N=avg_cor_posmat_N,
                      avg_cor_posmat_C=avg_cor_posmat_C,
                      avg_cor_posmat_SC=avg_cor_posmat_SC,
                      avg_cor_negmat_N=avg_cor_negmat_N,
                      avg_cor_negmat_C=avg_cor_negmat_C,
                      avg_cor_negmat_SC=avg_cor_negmat_SC,
                      avg_cor_btwyr_N=avg_cor_btwyr_N,
                      avg_cor_btwyr_C=avg_cor_btwyr_C,
                      avg_cor_btwyr_SC=avg_cor_btwyr_SC,
                      vr_LdM_N=vr_LdM_N,
                      vr_LdM_C=vr_LdM_C,
                      vr_LdM_SC=vr_LdM_SC)
saveRDS(avg_cor_matlist,paste(
  "../../Results/gather_res/avg_cor_matlist_nsp_",n,"_nyr_",N,".RDS",sep=""))

#------------ plot avg correlation between sp.-------------------
r1<-avg_cor_matlist$avg_cor_mat_N
r2<-avg_cor_matlist$avg_cor_mat_C
r3<-avg_cor_matlist$avg_cor_mat_SC

d1<-r1 %>% as.data.frame()%>%gather() 
d2<-r2 %>% as.data.frame()%>%gather()
d3<-r3 %>% as.data.frame()%>%gather()

g1<-ggplot(d1, aes(value)) + 
  geom_density() + 
  facet_wrap(~key)+
  xlab("Average correlation between species")+
  geom_density(data=d3,aes(value),col="blue",linetype="dotted",size=0.8)+
  geom_density(data=d2,aes(value),col="red",linetype="dashed")+
  theme_classic()

ggsave(g1, 
       filename = paste("../../Results/gather_res/avg_cor_btw_sp_at_different_r_nsp_",n,"_nyr_",N,".pdf",sep=""),
       device = cairo_pdf, 
       width = 15, height = 10, units = "cm")


# for LM synchrony

r1<-avg_cor_matlist$vr_LdM_N
r2<-avg_cor_matlist$vr_LdM_C
r3<-avg_cor_matlist$vr_LdM_SC

d1<-r1 %>% as.data.frame()%>%gather() 
d2<-r2 %>% as.data.frame()%>%gather()
d3<-r3 %>% as.data.frame()%>%gather()

g1<-ggplot(d1, aes(value)) + 
  geom_density() + 
  facet_wrap(~key)+
  xlab("vr_LdM between species")+
  geom_density(data=d3,aes(value),col="blue",linetype="dotted",size=0.8)+
  geom_density(data=d2,aes(value),col="red",linetype="dashed")+
  theme_classic()
ggsave(g1, 
       filename = paste("../../Results/gather_res/vr_LdM_at_different_r_nsp_",n,"_nyr_",N,".pdf",sep=""),
       device = cairo_pdf, 
       width = 15, height = 10, units = "cm")




# I run this code with species number n =20, n=25, n=30, n=40, and n=60
#


############################################
# pedagog figure
a1<-c(1,3,2,6,7,5,8,4,9,10)
a2<-c(1,2,3,10,4,9,7,6,5,8)

v1<-VineCopula::pobs(a1)
v2<-VineCopula::pobs(a2)
#v1<-1-v1
v2<-1-v2
d<-data.frame(Year=paste("sp",1:10,sep=""), v1, v2)
plot(v1,v2, pch=16, xlim=c(0,1),ylim=c(0,1))
text(d$v1, d$v2-0.05, labels=d$Year)

#######################################

# all these numbers should be similar
#dim(m1)
#cor.test(m1[1,,1],m1[4,,1],method="spearman")
#cor.test(c1[1,,1],c1[4,,1],method="spearman")

#plot(m1[3,,1],m1[4,,1])
#plot(c1[3,,1],c1[4,,1])

#range(m1)
#range(c1)




