rm(list=ls())
#---------------------
path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

if(!dir.exists("../../Results/gather_res/")){
  dir.create("../../Results/gather_res/")
}

source("Taylor.R") # gather all TL estimates and 
# play with exploratory data visualization

source("generate_com.R") # simulation shows stability vs. nsp plot 
#for different z and check in data below hypothesis from theory
# for z>2: pe_mv > pe_avg_cv   
# for z<2: pe_mv < pe_avg_cv

source("Simulation_zmorethan2.R") # new mechanism for z>2
source("test_empirical_z_taildep.R") # test with real data

source("synchrony_and_z.R")
