# function to compute average positive and negative 
# pairwise Spearman correlations between columns (species)\

# Arg:
# mat = a matrix or a dataframe where species time series along each column
library(Hmisc)
compute_avg_cor<-function(mat){
  
  id<-which(colnames(mat)=="raresp")
  if(length(id)>0){
    mat<-mat[,-id]
  }
  
  # compute correlation between columns (i.e. among species)?
  mat<-as.matrix(mat)
  
  cor_btw_col<-Hmisc::rcorr(mat,type="spearman") # correlation btw columns
  cormat<-cor_btw_col$r
  cormat.p<-cor_btw_col$P
  cormat.p.nonsig<-which(cormat.p>=0.05,arr.ind=T) # non-sig correlation
  cormat[cormat.p.nonsig]<-NA # so keeping only sig +ve or -ve cor between columns
  cormat[upper.tri(cormat,diag=T)] <- NA # only lower triangular part
  
  avg_cor_btw_sp<-mean(cormat,na.rm=T)
  
  id<-which(cormat>0, arr.ind = TRUE)
  if(nrow(id)>0){
    avg_cor_positive<-mean(cormat[id],na.rm = T) 
  }else{
    avg_cor_positive<-NA # no positive pairwise correlation found
  }
  
  id<-which(cormat<0, arr.ind = TRUE)
  if(nrow(id)>0){
    avg_cor_negative<-mean(cormat[id],na.rm = T)
  }else{
    avg_cor_negative<-NA # no negative pairwise correlation found
  }
  
  # Now correlation between years (or rows)
  tmat<-t(mat) # species by year matrix
  if(nrow(tmat)>4){ # more than 4 species required to estimate correlations
    corln<-Hmisc::rcorr(tmat,type="spearman") # correlation btw columns
    cortmat<-corln$r
    cortmat.p<-corln$P
    cortmat.p.nonsig<-which(cortmat.p>=0.05,arr.ind=T) # non-sig correlation
    cortmat[cortmat.p.nonsig]<-NA # so keeping only sig +ve or -ve cor between columns
    cortmat[upper.tri(cortmat,diag=T)] <- NA # only lower triangular part
    avg_cor_btw_yr<-mean(cortmat,na.rm=T)
  }else{
    avg_cor_btw_yr<-NA # only 4 species were not enough to estimate correlation
  }
  
  res<-data.frame(avg_cor_btw_sp=avg_cor_btw_sp,
                  avg_cor_pos_btw_sp=avg_cor_positive,
            avg_cor_neg_btw_sp=avg_cor_negative,
            avg_cor_btw_yr=avg_cor_btw_yr)
  return(res)
}
