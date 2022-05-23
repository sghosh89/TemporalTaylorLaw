# This script will do the copula analysis between years, will generate necessary plots
# for a given community

source("NonParamStat.R")
source("NonParamStat_matrixplot.R")

# Input:
# mat: a matrix or dataframe where each target species time series along each column, 
#         last column for all rest of the aggregated species, rows have name for sampling years
#resloc: path to save results
#nbin: 2 (default) to measure tail-dep.

# Output:
# a list and several plots to be saved in resloc path

tail_analysis<-function(mat, resloc, nbin=2){
  
  yrlist<-vector(mode="list",length=nrow(mat))
  names(yrlist)<-rownames(mat)
  mat<-as.data.frame(mat)# make sure to get a dataframe, otherwise in the 
  #next for loop Dat would be numeric and can't transpose
  for(i in 1:length(yrlist)){
    Dat<-mat[i,]
    Dat<-t(Dat)
    rownames(Dat)<-NULL
    yrlist[[i]]<-data.frame(sp=colnames(mat),Dat=Dat[,1])
  }
  
  z<-multcall(d_allyr = yrlist, resloc = resloc, nbin=nbin)
  
  tot_target_sp<-ncol(mat)
  nyr<-nrow(mat)# total number of years
  
  saveRDS(z,paste(resloc,"NonParamStat.RDS",sep=""))
  
  NonParamStat_matrixplot(data=z,
                          resloc=resloc,
                          tot_target_sp=tot_target_sp,
                          nyr=nyr,
                          tl.cex=1.2,cl.cex=2,line=1)
}