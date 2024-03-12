rm(list=ls())
#---------------------
path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

if(!dir.exists("../../Results/gather_res/")){
  dir.create("../../Results/gather_res/")
}

# for figures
source("conceptual_figs.R")
source("fig_with_interaction.R")
source("Taylor_all_visualization.R") # considering both realm altogether
source("test_empirical_z_taildep.R") # test with real data
source("Simulation_zmorethan2.R") # new mechanism for z>2

# raugh (conceptual fig 2)
# upper tail dep in conceptual figure
#a<-c(20:1)
#b<-c(20,18,19,17,16,15,14,13,6,12,3,5,11,9,1,10,7,8,4,2)
#op<-par(mar=c(0.5,0.5,0.5,0.5))
#plot(a,b,col="blue",pch=16,xlab="",ylab="",xlim=c(0,21),ylim=c(0,21),
#     xaxs="i",yaxs="i",xaxt="n", yaxt="n", asp=1, cex=2)
#par(op)
#text(a,b,20:1)

# lower tail dep in conceptual figure
#a<-c(20:1)
#b<-c(19,14,13,20,16,12,18,15,17,11:4,2,1,3)
#op<-par(mar=c(0.5,0.5,0.5,0.5))
#plot(a,b,col="red",pch=16,xlab="",ylab="",xlim=c(0,21),ylim=c(0,21),
#     xaxs="i",yaxs="i",xaxt="n", yaxt="n", asp=1, cex=2)
#par(op)
#text(a,b,20:1)

# no tail dep in conceptual figure
#a<-c(20:1)
#b<-c(16,14,13,20,17,12,18,15,19,5,4,2,10,11,1,8,9,6,7,3)
#plot(a,b,col="white")
#text(a,b,20:1)

