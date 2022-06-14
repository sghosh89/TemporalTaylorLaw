##############################################
# this code is tweked from original function
# SimHub1 of EcoVirtual package (random birth-death only)
# https://rdrr.io/cran/EcoVirtual/src/R/bioGeo.R
##############################################
library(EcoVirtual)#to load function "rich"

mySimHub1=function(S, j, D=1, cycles, m.weights=1, len,anima=TRUE){
  if(length(m.weights)>S){
    m.weigths <- m.weights[1:S]
    warning("length of m.weights truncated at end to match number of species (S)")
  }
  if(length(m.weights)<S){
    w1 <- rep(NA,S)
    w1 <- rep(m.weights, each=floor(S/length(m.weights)))
    w1[is.na(w1)] <- m.weights[length(m.weights)]
  }
  if(length(m.weights)==S) w1 <- m.weights
  if(cycles<200){cycles=200; cat("\n Minimum number of cycles: 200\n")}
  stepseq=round(seq(101, cycles+1, len=len))
  step=stepseq[2]- stepseq[1]
  
  ## Size of the community
  J <- S*j
  ## Matrices to save results
  ind.mat=matrix(nrow=J,ncol=100+length(stepseq)) 
  
  ## Initial conditions##
  ## All species start with the same number of individuals
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
  #################################################################
  ## incluindo 100 primeiros ciclos (including 100 first cycles)          ########
  #################################################################
  for(k in 2:100){
    ##Indice dos individuos que morrem
    # Index of individuals who die
    mortek <- sample(1:J,D, prob=w1[cod.sp])
    ##Indice dos individuos que produzem filhotes para substituir os mortos
    #Index of individuals that produce offspring to replace the dead
    novosk <- sample(1:J,D,replace=TRUE)
    ##Substituindo
    # replacing
    cod.sp[mortek]<-cod.sp[novosk]
    ind.mat[,k] <- cod.sp
  }
  ###########################
  cont=100
  tempo=0:99
  ##Aqui comecam as simulacoes
  # Here the simulations begin
  if(!is.null(stepseq)){
    for(i in 1:length(stepseq)){
      cont=cont+1
      for(j in 1:step){
        ##Indice dos individuos que morrem
        morte <- sample(1:J,D, prob=w1[cod.sp])
        ##Indice dos individuos que produzem filhotes para substituir os mortos
        novos <- sample(1:J,D,replace=TRUE)
        ##Substituindo
        cod.sp[morte]<-cod.sp[novos]
      }
      ## A cada step ciclos os resultados sao gravados
      # At each step cycles the results are recorded
      ind.mat[,cont] <- cod.sp
    }
    tempo <- c(tempo,stepseq)
  }
  colnames(ind.mat) <- tempo
  if(anima==TRUE){
    dev.new()
    animaHub(dadoHub=ind.mat)
  }
  #dev.new()
  #op = par(mar=c(6,5,5,2), las=1)
  #plot(as.numeric(colnames(ind.mat)),apply(ind.mat,2,rich), xlab="Time (cycles)", ylab="Number of species",ylim=c(0,S), cex.lab=1.2, type="l", col="red", lty=2,  main=paste("Neutral Model Without Colonization", "\n S=",S," J=",J),  sub=paste("Mean extintion=",(S-rich(ind.mat[,ncol(ind.mat)]))/cycles,"sp/cycle"), cex.sub=1, cex.axis=1.2, cex.lab=1.2, lwd=2) 
  #invisible(ind.mat)
  return(ind.mat)
}
#############################################
source("get_spabund_mat.R")
#############################################
# call the function to get 100 replicates

#nrep= number of replicates
#S= number of species
#j= number of starting individual per species
#cycles= cycles in neutral models, total timesteps
#len= number of timesteps you want to record the dynamics after initial 100 steps

call_mySimHub1<-function(resloc,nrep=100,S=30,j=500,cycles=1e5,len=1e2){
  
  J=S*j
  for(i in 1:nrep){
    
    # compute neutral dynamics
    res<-mySimHub1(S = S, j = j, D = 1, cycles = cycles, len=len,
                   m.weights = 1, anima = F)
    saveRDS(res,paste(resloc,"res_simhub1_nrep_",i,".RDS",sep=""))
    pdf(paste(resloc,"richness_with_t_nrep_",i,".pdf",sep=""),height=6,width=6)
    op = par(mar=c(6,5,5,2), las=1)
    plot(as.numeric(colnames(res)),apply(res,2,rich), 
         xlab="Time (cycles)", ylab="Number of species",
         ylim=c(0,S), cex.lab=1.2, type="l", col="red", lty=2,  
         main=paste("Neutral Model Without Colonization", 
                    "\n S=",S," J=",J),  sub=paste("Mean extintion=",(S-rich(res[,ncol(res)]))/cycles,"sp/cycle"), cex.sub=1, cex.axis=1.2, cex.lab=1.2, lwd=2) 
    par(op)
    dev.off()
    
    # compute year by species matrix
    res2<-get_spabund_mat(spmat=res,S)
    saveRDS(res2,paste(resloc,"yr_by_sp_nrep_",i,".RDS",sep=""))
    
    cat("-------- nrep = ",i," --------------------- \n")
  }
}
########################

