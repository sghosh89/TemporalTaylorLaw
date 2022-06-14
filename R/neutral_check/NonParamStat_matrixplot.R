# This function is written to generate matrix plot for non-parametric stat results 
# Input : 
#   data 
#   resloc : location to save the results
#   posn_notneeded: a matrix of cell indices to not show in the plot
#   tot_target_sp: number of targetted sp (generally common sp)
#     tl.cex,cl.cex,line : numeric values for size of text in the labels and colorbar, line of mtext
#---------------------------
source("mycorrplot.R")
#---------------------------

NonParamStat_matrixplot<-function(data,nyr,tot_target_sp,resloc,tl.cex,cl.cex,line){
  
  #--------------------------Spearman plot---------------------------
  
  tempo<-data$spear
  indI<-data$posnI
  indI<-which(indI==1,arr.ind = T)
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  minval<-min(tempo,na.rm=T)
  maxval<-max(tempo,na.rm=T)
  
  cr<-max(abs(minval),abs(maxval))
  
  #pdf(paste(resloc,file="Spearman.pdf",sep=''),width=wd, height=ht)
  png(paste(resloc,file="Spearman.png",sep=''), width = 5000,height = 5000,  res = 300)
  posnI_ind<-which(data$posnI==1,arr.ind=T)
  posnN_ind<-which(data$posnN==1,arr.ind=T)
 
  mycorrplot(z=tempo,
             posnI_ind=posnI_ind,
             posnN_ind=posnN_ind,
             colrange=c(0,cr),tl.cex=tl.cex,cl.cex=cl.cex)
  dev.off()
  
  
  #--------------------------Kendall plot---------------------------
  
  tempo<-data$kend
  indI<-data$posnI
  indI<-which(indI==1,arr.ind = T)
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  minval<-min(tempo,na.rm=T)
  maxval<-max(tempo,na.rm=T)
  
  cr<-max(abs(minval),abs(maxval))
  
  
  #pdf(paste(resloc,file="Kendall.pdf",sep=''),width=wd, height=ht)
  png(paste(resloc,file="Kendall.png",sep=''), width = 5000,height = 5000,  res = 300)
  posnI_ind<-which(data$posnI==1,arr.ind=T)
  posnN_ind<-which(data$posnN==1,arr.ind=T)
  mycorrplot(z=tempo,
             posnI_ind=posnI_ind,
             posnN_ind=posnN_ind,
             colrange=c(0,cr),tl.cex=tl.cex,cl.cex=cl.cex)
  dev.off()
  
  
  #========================================= For cor npa stats ===============================================
  
  #--------------------------Corl plot---------------------------
  
  tempo<-data$Corl
  indI<-data$posnI
  indI<-which(indI==1,arr.ind = T)
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  
  minval<-min(tempo,na.rm=T) # Corl will be always >=0
  maxval<-max(tempo,na.rm=T)
  
  cr<-max(abs(minval),abs(maxval))
  
  #pdf(paste(resloc,file="Corl.pdf",sep=''),width=wd, height=ht)
  png(paste(resloc,file="Corl.png",sep=''), width = 5000,height = 5000,  res = 300)
  posnI_ind<-which(data$posnI==1,arr.ind=T)
  posnN_ind<-which(data$posnN==1,arr.ind=T)
  mycorrplot(z=tempo,
             posnI_ind=posnI_ind,
             posnN_ind=posnN_ind,
             colrange=c(-cr,cr),tl.cex=tl.cex,cl.cex=cl.cex) 
  # theoretically should be 0 to cr range 
  # but sometimes, you get slight negative correlation, if their are many ties
  dev.off()
  
  
  #--------------------------Coru plot---------------------------
  
  tempo<-data$Coru
  indI<-data$posnI
  indI<-which(indI==1,arr.ind = T)
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  
  minval<-min(tempo,na.rm=T)# Coru will be always >=0
  maxval<-max(tempo,na.rm=T)
  
  cr<-max(abs(minval),abs(maxval))
  
  #pdf(paste(resloc,file="Coru.pdf",sep=''),width=wd, height=ht)
  png(paste(resloc,file="Coru.png",sep=''), width = 5000,height = 5000,  res = 300)
  posnI_ind<-which(data$posnI==1,arr.ind=T)
  posnN_ind<-which(data$posnN==1,arr.ind=T)
  mycorrplot(z=tempo,
             posnI_ind=posnI_ind,
             posnN_ind=posnN_ind,
             colrange=c(-cr,cr),tl.cex=tl.cex,cl.cex=cl.cex)
  dev.off()
  
  #--------------------------Corl-Coru plot---------------------------
  
  tempo<-data$Corl - data$Coru
  indI<-data$posnI
  indI<-which(indI==1,arr.ind = T)
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  minval<-min(tempo,na.rm=T)
  maxval<-max(tempo,na.rm=T)
  
  cr<-max(abs(minval),abs(maxval))
  
  #pdf(paste(resloc,file="Corl-Coru.pdf",sep=''),width=wd, height=ht)
  png(paste(resloc,file="Corl-Coru.png",sep=''), width = 5000,height = 5000,  res = 300)
  z<-tempo
  posnI_ind<-which(data$posnI==1,arr.ind=T)
  posnN_ind<-which(data$posnN==1,arr.ind=T)
  mycorrplot(z=z,
             posnI_ind=posnI_ind,
             posnN_ind=posnN_ind,
             colrange=c(-cr,cr),tl.cex=tl.cex,cl.cex=cl.cex)
  
  
  #------------------- for target sp matrix ----------------------------------
  
  z[posnN_ind]<-NA # this line was added to exclude -vely correlated years from nL,nU,L,U 
  # calculation, but it does not matter as for -vely correlated cells [yr_i,yr_j] and 
  # [yr_j,yr_i] nL,nU both will increase by same number, whereas L+U remains same and
  # both L, U will change by same +, - factor
  
  
  z[upper.tri(z)]<-NA # symmetric matrix, so consider only lower triangular part
  
  nL<-sum(z>0,na.rm = T)
  nU<-sum(z<0,na.rm = T)
  L<-sum(z[which(z>0,arr.ind=T)])
  U<-sum(z[which(z<0,arr.ind=T)])
  
  
  npos<-nL+nU # number of positive correlation
  nint<-nrow(z)*(nrow(z)-1)/2 # number of all pairwise interaction
  nind<-nrow(posnI_ind) # number of indep. interaction
  
  nneg<-nrow(posnN_ind) # negative correlaion in lower.tri of a symmetric matrix
  
  
  summary_df<-data.frame(nsp=tot_target_sp,nyr=nyr,
                         nint=nint,nind=nind,npos=npos,nL=nL,nU=nU,
                         nneg=nneg,L=L,U=U)
  
  #------------------------------------------------------------------------------------------
  
  #mtext(paste0("nL =",nL,", nU =",nU),cex=5,side=1,adj=0.7)
  mtext((as.expression(bquote('N'['+']*' = '*.(npos)))),cex=tl.cex,line=line,side=1,adj=0.2,col="gray")# number of cells with pos. correlation
  mtext((as.expression(bquote(', N'['L']*' = '*.(nL)))),cex=tl.cex,line=line,side=1,adj=0.4,col="red")# number of cells with Lower or Left tail dep.
  mtext((as.expression(bquote(', N'['R']*' = '*.(nU)))),cex=tl.cex,line=line,side=1,adj=0.6,col="blue")# number of cells with Upper or Right tail dep.
  mtext((as.expression(bquote(', N'['-']*' = '*.(nneg)))),cex=tl.cex,line=line,side=1,adj=0.8,col="green")# number of cells with neg. correlation
  
  
  dev.off()
  
  res<-list(CorlmCoru=z,
            summary_df=summary_df)
  
  saveRDS(summary_df,paste(resloc,"summary_df.RDS",sep=""))
  
  #return(res)
}