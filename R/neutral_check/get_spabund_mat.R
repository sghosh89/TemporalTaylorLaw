# function to get the year by species matrix 
# from the output of neutral model, mySimHub1

# Argument:
# spmat = a matrix, output from neutral model, 
# where rows are the total inidividuals present in the community
# and each column represents each timesteps - first 100 columns are 
# 1:100, then the timesteps are not 1 anymore.

# S = numeric value, integer, species richness to start with
library(dplyr)
get_spabund_mat<-function(spmat,S){
  
  # initiate year by sp matrix
  res<-matrix(NA,nrow=ncol(spmat),ncol=S)
  colnames(res)<-paste("sp",1:S,sep="")
  rownames(res)<-paste("t=",colnames(spmat),sep="")
  
  # fill year by sp matrix from spmat
  z<-apply(spmat, MARGIN=2, FUN= function(x){as.data.frame(table(x))})
  
  reft<-data.frame(x=as.factor(1:S),Freq=0) #reference table
  for(i in 1:length(z)){
    temp<-z[[i]]
    if(nrow(temp)<S){
      temp2<-right_join(temp,reft,by=c("x"))
      temp2$new<-coalesce(temp2$Freq.x,temp2$Freq.y)
      temp<-data.frame(x=temp2$x,Freq=temp2$new)
    }
    res[i,]<-temp$Freq
  }
  return(res)
}





