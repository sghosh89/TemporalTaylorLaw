# Select any two years [i,j] for a copula with all species' abundance/ biomass
# This function gives you a matrix with vi and vj as two columns

# Input :
# d_allyr : species abundance/biomass dataset in data[[yr]] format
# yr : year index
# i,j : pair of year indices
# level : significance level for BiCopIndepTest p-value
# ploton : (optional) logical, if T gives copula plot without transforming j-th 
# variable to it's -ve value

# Output :
# A list of 4 elements:
#                      mat : a matrix : copula of (vi,vj) with transforming j-th variable to it's -ve value for -ve corr.
#                      corval : Spearman's correlation
#                      pval   : pvalue of Kendall's cor.test
#                      IndepTestRes : BiCopIndepTest p-value

# and an optional plot of the copula

library(VineCopula)
vivj_matrix<-function(d_allyr,i,j,level=0.05,ploton){
  
  ds1<-d_allyr[[i]]
  ds2<-d_allyr[[j]]
  #----------------------------
  colnames(ds1)<-c("Sp","Dat")  # ensuring column names
  colnames(ds2)<-c("Sp","Dat")
  
  # Omitting the years and data containing NA in either d1 or d2 
  #from both d1 and d2
  if(anyNA(ds1$Dat)==T | anyNA(ds2$Dat)==T){
    dboth<-cbind(ds1,ds2$Dat)
    dboth<-na.omit(dboth)
    d1<-data.frame(Sp=dboth$Sp,Dat=dboth[,2])
    d2<-data.frame(Sp=dboth$Sp,Dat=dboth[,3])
  }else{
    d1<-ds1
    d2<-ds2
  }
  
  colnames(d1)[2]<-"Dat"  # ensuring column names
  colnames(d2)[2]<-"Dat"
  
  #get ranks modified now
  vi<-VineCopula::pobs(d1$Dat)
  vj<-VineCopula::pobs(d2$Dat)
  
  IndepTestRes<-VineCopula::BiCopIndTest(vi,vj)$p.value
  ct<-cor.test(vi,vj,alternative = "two.sided",method="spearman",exact=F)
  corval<-unname(ct$estimate)
  pval<-ct$p.value
  
  # should only return significant corval (either positive or negative)
  
  if(IndepTestRes<level && corval>0){ # for significant positive correlation
    corval<-corval
    if(ploton==T){
      plot(vi,vj,type='p',col=rgb(0,0,0,0.3),pch=19,xlim=c(0,1),ylim=c(0,1),
           xlab=names(d_allyr)[i],ylab=names(d_allyr)[j],cex.lab=1.5)
      mtext(paste0("(yr_x, yr_y) = (",i," , ",j,")"),
            side = 3, line=0.15, adj=0.5, col="black")
    }
    
  }else if(IndepTestRes<level && corval<0){ # for significant negative correlation
    corval<-corval
    if(ploton==T){
      plot(vi,vj,type='p',col=rgb(0,1,0,0.3),pch=19,xlim=c(0,1),ylim=c(0,1),
           xlab=names(d_allyr)[i],ylab=names(d_allyr)[j],cex.lab=1.5)
      mtext(paste0("(yr_x, yr_y) = (",i," , ",j,")"),
            side = 3, line=0.15, adj=0.5, col="black")
    }
    vj<-VineCopula::pobs(-(d2$Dat)) # reverse the variable
    
  }else{ # independent case
    corval<-NA # non-significant
    if(ploton==T){
      plot(-1,0,xlim=c(0,1),ylim=c(0,1),xlab=names(d_allyr)[i],ylab=names(d_allyr)[j],cex.lab=1.5)
      text(0.5,0.5,"Indep.",adj=c(0.5,.5),cex=2)
      mtext(paste0("(yr_x, yr_y) = (",i," , ",j,")"),
            side = 3, line=0.15, adj=0.5, col="black")
    }
  }
  
  
  Years<-d1$Year
  #-------------------------
  #n_datapt<-length(vi)
  #--------------------
  #plot(vi,vj,type="p")
  #-------------------------
  mat<-as.matrix(cbind(vi,vj))
  return(list(mat=mat,   # return reversed mat so that if you plot this mat you get +ve correlation 
              corval=corval,  # but return the actual -ve corr. value  
              pval=pval,
              IndepTestRes=IndepTestRes))  
}