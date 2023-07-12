# Load packages

library(doParallel)
library(doRNG)
library(rxode2)
library(dplyr)
library(tidyr)
library(adaptivetau)


# Parallellization

nodelist<-rep("localhost",30) # 32 cores
cl<-makePSOCKcluster(nodelist) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())

clusterEvalQ(cl, setwd("/gxfs_home/cau/sunzm666")) # set work dir correctly to each cluster instance


# Call function for main model structure
source("Model_function_Junqi_Adaptivetau.R")


###
# Scenario 83
# fast sequential treatment
###########

Scenario_n <- 83

#Create folder
dir.create(paste0("Scenario_",Scenario_n))

#define input 
t_half <- 2
n = 100
v_CS_A  <- 1
v_CS_B  <- 1
v_FIT <- 0.8
v_HILL_A <- c(0.5, 3)
v_HILL_B <- c(0.5, 3)
v_Gmin <- -6
v_RA <- 0
v_RB <- 0
v_eB <- 6
v_Css <- c(4, 8)
v_U_1 <- c(10^-6, 10^-8)
v_U_2 <- c(10^-6, 10^-8)

input_all <- expand_grid(v_HILL_A, v_HILL_B, v_FIT, v_Gmin, v_RA, v_RB, v_Css, v_U_1, v_U_2, v_eB) %>% 
  arrange(desc(v_FIT), v_HILL_A, v_HILL_B) %>% 
  mutate(ID = row_number()) %>% 
  as.data.frame()

saveRDS(input_all, file = paste0("Scenario_", Scenario_n, "/" ,"Scenario_", Scenario_n, "_input.rds" ))

select_ID <- 1:length(input_all$ID)

for(i in 1:length(select_ID )){
  
  id_i <- select_ID[i] 
  
  input <- input_all %>% 
    filter(ID %in% select_ID)
  
  sim <- CS_model(  v_Models = "1 day cycling",
                    t_half  =  2,
                    FIT     =  input$v_FIT[i],
                    HILL_A	=  input$v_HILL_A[i],
                    HILL_B	=  input$v_HILL_B[i],
                    Gmin_A	=  input$v_Gmin[i],
                    Gmin_B	=  input$v_Gmin[i],
                    MIC_S   =  1,
                    MIC_R   =  8,
                    u_1 = input$v_U_1[i], 
                    u_2 = input$v_U_2[i],
                    u_3 = 10^-5,
                    u_4 = 10^-5,
                    eB0  = input$v_eB[i],
                    F_Css_MIC = input$v_Css[i],
                    RA0 = input$v_RA[i],
                    RB0 = input$v_RB[i],
                    RAB0 = 0,
                    Bmax = 9,
                    V    = 1000,
                    n    = n,
                    ST   = 24*7,
                    DT = 1) %>% 
    mutate(SIM_ID = id_i,
           n = n)
  
  saveRDS(sim, file = paste0("Scenario_", Scenario_n, "/" ,"Scenario_", Scenario_n, "_sim_",id_i, "_n",n, "_", Sys.Date(), ".rds" ))
  
}

stopCluster(cl)

