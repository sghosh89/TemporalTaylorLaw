# This script checks theoretical prediction that with increasing richness 
# community variability should decrease (stability should increase), 
# but the rate of decrease in cv or increase in icv
# will be faster for lower z than the higher z
#======================================================
rm(list=ls()) 
library(ecofolio)
#z in the range of 0.6 to 2.8
calcsd2<- function (mu,z) sigmasq<- 0.01*mu^z # variance ~ constant * mean^z

#a_low<-rnorm(n=100, mean=50, sd=sqrt(calcsd2(50, 0.6)))
#a_high<-rnorm(n=100, mean=50, sd=sqrt(calcsd2(50, 1.5)))
#mean(a_high)
#var(a_high)
#plot(1:100, a_high, type="l")
#lines(1:100, a_low, col=2)

# generate a community
# nyr = number of years
# nsp = number of species in each community

gen_com_stab<-function(nyr, nsp, mu, z){
  
  eachmat<-matrix(NA, nrow=nyr, ncol=nsp)
  rownames(eachmat)<-paste("year",1:nyr,sep="")
  colnames(eachmat)<-paste("sp",1:nsp,sep="")
  for(i in 1:nsp){
    eachmat[,i]<- rnorm(n=nyr, mean=mu, sd=sqrt(calcsd2(mu,z)))
  }
  tot_ts<-apply(eachmat,MARGIN = 1, FUN=sum)
  res<-ecofolio::cv(tot_ts)
  res<-1/res
  return(res)
}

set.seed(seed=123)
# now generate data for CV vs. nsp for different z
sp<-c(2,4,8,16,32,64)
#zall<-c(0.6,1.0,1.5,2.0,2.8)
zall<-c(1.5,1.8,2.0,2.4,2.8)
icvmat<-matrix(NA,nrow=length(sp),ncol=length(zall))
colnames(icvmat)<-paste("z_",zall,sep="")
rownames(icvmat)<-paste("nsp_",sp,sep="")
for(j in 1:length(zall)){
  z<-zall[j]
  for(i in 1:length(sp)){
    nsp<-sp[i]
    res<-gen_com_stab(nyr=100,nsp=nsp,mu=50,z=z)
    icvmat[i,j]<-res
  }
}
matplot(sp,icvmat,type="b",xlab="Richness",ylab="Stability_com") # 1 to 5 higher z

#====================================================================
# Here, we will test for a range of z (Taylor's slope) when the 
# pe_avg_cv and pe_mv differ most
# Theoretically, I can see 
# for z=2: both would be same
# for z>2: pe_mv > pe_avg_cv   
# for z<2: pe_mv < pe_avg_cv
#====================================================================
# first we will check with real data if it's true?
sm_all<-readRDS("../../Results/gather_res/TaylorEstimate_alldata.RDS")

sm_all$pe_diff<-sm_all$pe_avg_cv - sm_all$pe_mv

df_1<-sm_all%>%filter(TLslope.z<2)
hist(df_1$pe_diff,100,xlab="pe_avg_cv - pe_mv",main="when z<2, 1676 obs.") # as expected for z<2: pe_mv < pe_avg_cv

df_2<-sm_all%>%filter(TLslope.z>=2)
hist(df_2$pe_diff,100,xlab="pe_avg_cv - pe_mv",main="when z>=2, 87 obs.") # as expected for z>2: pe_mv > pe_avg_cv
range(df_2$pe_diff)
sum(df_2$pe_diff>0)/nrow(df_2)*100 # nearly 5.7% data deviates, it's ok


# theoretically, when z=1
pe_mv_num<-(50*51*0.5)^(-0.5)
pe_avg_cv_num<-0
for(mu in 1:50){
  temp<-mu^(-0.5)
  pe_avg_cv_num<-pe_avg_cv_num+temp
}
pe_avg_cv_num<-pe_avg_cv_num/50
pe_avg_cv_num > pe_mv_num #checked



