# check rank abundance curve
#rm(list=ls())
library(tidyverse)
library(gridExtra)
library(BiodiversityR)
##############################
#resloc<-"../../Results/neutral_check/mySimHub1_S_30_j_100_cycles_1e+05/"
#resloc<-"../../Results/neutral_check/mySimHub1_S_30_j_500_cycles_1e+05/"
#resloc<-"../../Results/neutral_check/mySimHub1_S_30_j_500_cycles_1e+08/"
check_rank_abund<-function(resloc){
  nrep<-100
  p<-list() # to hold 100 ggplot as a list 
  pall<-list()# to hold 100 ggplot as a list for total abund
  for(i in 1:nrep){
    x<-readRDS(paste(resloc,"yr_by_sp_nrep_",i,".RDS",sep=""))
    x<-tail(x,26) # a year by sp matrix
    
    xr<-rankabundance(x) # computes each species rank based on abundance over the years
    # rank 1 means that species was abundant most across the years
    xr<-as.data.frame(xr)
    
    pall[[i]]<-ggplot(xr,aes(x=rank,y=abundance,label=rownames(xr)))+geom_line(col="gray",size=2)+
      xlab("Species rank")+ geom_text()+
      ylab("total abundance over the years")+
      ggtitle(paste("nrep = ",i,sep=""))+
      theme_classic()
    
    # or you can compute rank for each species for each year
    #xr2<-apply(x,MARGIN = 1, FUN=rank) # but here rank 1 means rarest
    # so flip the scale so that most abundant gets rank 1
    #xr2<-apply(x,MARGIN = 1, FUN=function(a){as.integer(rank(-a))})
    
    # check for each year all the species distribution 
    # it should be non-uniform
    x<-as.data.frame(x)
    S<-ncol(x) # number of species
    x$time<-rownames(x)
    df<-x %>% tidyr::gather("id", "value", 1:S)
    df$id<-factor(df$id, levels=c(paste("sp",1:S,sep=""))) # keep the sequence fixed
    # Define the number of colors you want
    
    p[[i]]<-ggplot(df,aes(x=id,y=value,group=time,col=time))+
      geom_smooth(se=F)+ggtitle(paste("nrep = ",i,sep=""))+
      xlab("id")+ylab("abundance, col=year")+  theme_classic()+
      theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  ########################
  ggsave(
    filename = paste(resloc,"check_rank_abundplots_eachyr.pdf",sep=""), 
    plot = marrangeGrob(p, layout_matrix = matrix(1:100,  nrow = 10, ncol=10, byrow=TRUE)), 
    width = 40, height = 40
  )
  ggsave(
    filename = paste(resloc,"check_rank_abundplots_acrossyr.pdf",sep=""), 
    plot = marrangeGrob(pall, layout_matrix = matrix(1:100,  nrow = 10, ncol=10, byrow=TRUE)), 
    width = 40, height = 40
  )
}
