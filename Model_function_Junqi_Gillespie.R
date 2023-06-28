
# PK-PD model for INC, SAME and SEQ treatments
# Refer to Linda B. S. Aulin, et al. 2021 (General structure and PK)
# Refer to W Wang, et al. 2015 (PK)
# Refer to Christin Nyhoegen, et al. 2023 (PD)

# Create variables, actual values were inputted in running script 
CS_model <- function(t_half = 2, F_Css_MIC = 1, 
                     v_Models = c("Mono A", 
                                  "Mono B",
                                  "3 day cycling",  
                                  "1 day cycling", 
                                  "Multi-step evolution"), 
                     FIT = 1, Gmin_A = -6, Gmin_B = -6, HILL_A  = 1, HILL_B = 1, u_1 = 10^-6, u_2 = 10^-6, u_3 = 10^-6, u_4 = 10^-6, eB0 = 6, RA0=0, RB0 =0 , RAB0 = 0,  Bmax = 9,  V = 1000, n = 100, DT = 1, ST = 24, MIC_S = 1, MIC_R =8, Gmax = 1.5)
  
# Mono: Single drug with constant concentrations (SAME evolution in main text)
# 3 day cycling: Switching drug every 3 days
# 1 day cycling: Switching drug every 1 day
# Multi-step evolution: Single drug with increasing concentrations (INC evolution in main text)
  
# argument explanation
# t_half; half life used for both drugs [h]
# F_Css_MIC scaling; factor to obtain Css average equal to X time MIC WT
# FIT; Fitness cost (0-1), 1 no cost, 0 no growth 
# Gmin_A; Maximal effect of drug A
# Gmin_B; Maximal effect of drug B
# HILL_A; Shape factor drug A
# HILL_B; Shape factor drug B
# u_1; Mutation rate of S against antibiotic A
# u_2; Mutation rate of S against antibiotic B
# u_3; Mutation rate of RB against antibiotic A
# u_4; Mutation rate of RA against antibiotic B
# eB0; Starting bacterial concentration [cfu/mL]
# RA0; Fraction of bacterial starting in RA
# RB0; Fraction of bacterial starting in RB
# RAB0; Fraction of bacterial starting in RAB
# Bmax; carrying capacity of the system [cfu/mL]
# V; volume of infection site
# n; number of simulations , 
# DT; time step used [h]
# ST; Duration of simulation [h]



{
  require(doParallel)
  require(doRNG)
  require(rxode2)
  require(dplyr)
  require(tidyr)
  require(adaptivetau)
  
  
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  # Used for integer check in output calculation 
  # (Such as when we want to calculate the results only by day instead of hour)
  
  start_t <- Sys.time()  #for checking run time
  
  ##############
  ## PK Model ##
  ##############
  
  #Define PK parameters and dosing, assume same PK for both drugs
  
  
  Vd = 5                  # [L] volume of distribution, based on plasma volume
  ke = log(2)/t_half      # [h^-1] elimination rate constant, calculated by half-life
  CL = ke*Vd              # [L/h] clearance 
  Tau = 12                # [h] dosing interval
  MIC = 1                 # [ug/ul] default MIC value
  
  Css = MIC*F_Css_MIC     # [mg/L] Target average steady-state concentration 
  
  Dose = Css*CL*Tau/Vd    # [mg/L] dose giving target Css
  
  
  PD_mod<- rxode2({
    d/dt(A) = -ke*A   
    d/dt(B) = -ke*B
    d/dt(S) = S * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax) > 0, 
                          ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax), 0) - 
                     (((1 - Gmin_A/Gmax)*(A_t/MIC_S)^HILL_A/((A_t/MIC_S)^HILL_A - (Gmin_A/Gmax)))   +
                        ((1 - Gmin_B/Gmax)*(B_t/MIC_S)^HILL_B/((B_t/MIC_S)^HILL_B - (Gmin_B/Gmax))))*Gmax)
    d/dt(RA) = RA * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax) > 0, 
                            ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax), 0) * FIT - 
                       (((1 - Gmin_A/Gmax)*(A_t/MIC_R)^HILL_A/((A_t/MIC_R)^HILL_A - (Gmin_A/Gmax)))   +
                          ((1 - Gmin_B/Gmax)*(B_t/MIC_S)^HILL_B/((B_t/MIC_S)^HILL_B - (Gmin_B/Gmax))))*Gmax*FIT)
    d/dt(RB) = RB * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax) > 0, 
                            ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax), 0) * FIT - 
                       (((1 - Gmin_A/Gmax)*(A_t/MIC_S)^HILL_A/((A_t/MIC_S)^HILL_A - (Gmin_A/Gmax)))   +
                          ((1 - Gmin_B/Gmax)*(B_t/MIC_R)^HILL_B/((B_t/MIC_R)^HILL_B - (Gmin_B/Gmax))))*Gmax*FIT)
    d/dt(RAB) = RAB * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax) > 0, 
                              ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax), 0) * FIT^2 - 
                         (((1 - Gmin_A/Gmax)*(A_t/MIC_R)^HILL_A/((A_t/MIC_R)^HILL_A - (Gmin_A/Gmax)))   +
                            ((1 - Gmin_B/Gmax)*(B_t/MIC_R)^HILL_B/((B_t/MIC_R)^HILL_B - (Gmin_B/Gmax))))*Gmax*FIT^2)
    
    
  });
  
  
  # Dosing regimens
  # specified using argument v_Models
  # PK_mod was record to initial drug concentration in PD model
  
  PK_mod<- rxode2({
    
    
    d/dt(A) = -ke*A   
    d/dt(B) = -ke*B
    
  });
  
  
  PK_inits <- c(A = 0, B = 0);
  
  PK_theta <- c(ke = ke)
  
  
  PK_model_list <- list()
  
  if( "Mono A" %in% v_Models){  
    
    ev_PK_A <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(dosing.to = 1, dose=Dose, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST, DT))  %>%
      as.tbl()
    
    PK_A      <- as.data.frame(PK_mod$run(PK_theta, ev_PK_A,     PK_inits)) %>% 
      mutate(Model= "Mono A")
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_A  }
  
  
  if( "Mono B" %in% v_Models){
    
    ev_PK_B <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(dosing.to = 2, dose=Dose, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST, DT))  %>%
      as.tbl()
    
    PK_B      <- as.data.frame(PK_mod$run(PK_theta, ev_PK_B, PK_inits)) %>% 
      mutate(Model = "Mono B") 
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_B }
  
  
  
  if( "3 day cycling" %in% v_Models){
    
    
    ev_PK_3day <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(start.time = 0 ,  dosing.to = 1, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 72,  dosing.to = 2, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 144, dosing.to = 1, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 216, dosing.to = 2, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 288, dosing.to = 1, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST,DT))%>%
      as.tbl()
    
    PK_3day   <- as.data.frame(PK_mod$run(PK_theta, ev_PK_3day,  PK_inits))%>% 
      mutate(Model = "3 day cycling")
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_3day
    
  }
  
  
  if( "1 day cycling" %in% v_Models){
    
    
    ev_PK_1day <- eventTable(amount.units="mg", time.units="hours") %>%   #2 dose cycling
      add.dosing(start.time = 0,     dosing.to =  1, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.dosing(start.time = Tau,   dosing.to =  1, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.dosing(start.time = Tau*2, dosing.to =  2, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.dosing(start.time = Tau*3, dosing.to =  2, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.sampling(seq(0,ST,DT))  %>%
      as.tbl()
    
    
    PK_1day   <- as.data.frame(PK_mod$run(PK_theta, ev_PK_1day,  PK_inits))%>% 
      mutate(Model = "1 day cycling")
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_1day
    
  }
  
  if( "Multi-step evolution" %in% v_Models){
    
    
    ev_PK_multi <- eventTable(amount.units="mg", time.units="hours") %>%   #Multi-step
      add.dosing(start.time = 0,   dosing.to =  1, dose=0.5*Dose, nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 24,  dosing.to =  1, dose=0.5*Dose, nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 48,  dosing.to =  1, dose=1*Dose,   nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 72,  dosing.to =  1, dose=2*Dose,   nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 96,  dosing.to =  1, dose=4*Dose,   nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 120, dosing.to =  1, dose=8*Dose,   nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 144, dosing.to =  1, dose=16*Dose,  nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 168, dosing.to =  1, dose=32*Dose,  nbr.doses=2, dosing.interval=Tau) %>%
      add.dosing(start.time = 192, dosing.to =  1, dose=64*Dose,  nbr.doses=2, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST,DT))  %>%
      as.tbl()
    
    
    PK_multi   <- as.data.frame(PK_mod$run(PK_theta, ev_PK_multi,  PK_inits))%>% 
      mutate(Model = "Multi-step evolution")
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_multi
    
  }
  
  
  
  all_models <- c("Mono A", "Mono B", "3 day cycling",  "1 day cycling", "Multi-step evolution")
  dose_reg <- all_models[all_models %in% v_Models]  # selecting only models in dosing regimens in v_Models
  
  
  ##############
  ## PD Model ##
  ##############
  
  
  S0  <- (10^eB0)*(1-RA0 - RB0 - RAB0)  # starting bacterial density [cfu/ml] of sensitive bacteria (S) = WT)
  U_1 <- u_1  # Mutation rate [mutation/bacteria/hour] against antibiotic A
  U_2 <- u_2  # Mutation rate [mutation/bacteria/hour] against antibiotic B
  U_3 <- u_3  # Mutation rate [mutation/bacteria/hour] of Resistance recover for A
  U_4 <- u_4  # Mutation rate [mutation/bacteria/hour] of Resistance recover for B
  
  n_obs  <- ST*DT+1   # number of observations ( simulation time(ST) * time step size (DT))
  n_models <- length(PK_model_list) #which models (i.e. dosing regimens) to simulate
  
  
  df_full_CS <-  foreach( ii = 1:n, .combine = "rbind",              # Parallelization
                          .packages = c("rxode2", "dplyr", "tidyr", "adaptivetau"),
                          .inorder = F,
                          # Only when the export order is not important
                          .options.RNG = 123) %dorng%  {            
                            # Setting seed to make sure we get identical values every time (pseudo random number)
                            
                            source("PK_list.R")
                            
                            # i_time <- Sys.time()
                            # Creating empty data frame to combined simulation data per realization ii
                            df_model_CS <- data.frame(time  = NA,
                                                      S     = NA,
                                                      RA    = NA,
                                                      RB    = NA,
                                                      RAB   = NA,
                                                      A     = NA,
                                                      B     = NA,
                                                      model = NA)
                            
                            df_CS <- data.frame(time = rep(x = NA, times = n_obs),
                                                S    = rep(x = NA, times = n_obs),
                                                RA   = rep(x = NA, times = n_obs),
                                                RB   = rep(x = NA, times = n_obs),
                                                RAB  = rep(x = NA, times = n_obs),
                                                A    = rep(x = NA, times = n_obs),
                                                B    = rep(x = NA, times = n_obs))
                            
                            
                            
                            
                            
                            
                            
                            
                            for(i_mod in 1:length(PK_model_list)) {
                              
                              PK_model = PK_model_list[[i_mod]]
                              
                              # Set up starting data (population t = 0)
                              df_CS[1,] <- data.frame(time = 0, S = S0, RA = RA0*S0, RB = RB0*S0, RAB = RAB0*S0,
                                                      A = PK_model$A[1], B = PK_model$B[1]);
                              
                              V = V
                              rS_RA <- u_1
                              rS_RB <- u_2
                              rRA_RAB <- u_2
                              rRB_RAB <- u_1
                              rRA_S <- u_3
                              rRB_S <- u_4
                              rRAB_RA <- u_4
                              rRAB_RB <- u_3
                              
                              # initial parameters in time 0
                              S_i = floor(S0*V)
                              RA_i = floor(RA0*V)
                              RB_i = floor(RB0*V)
                              RAB_i = floor(RAB0*V)
                              A_t = PK_model$A[1]
                              B_t = PK_model$B[1]
                              
                              # t is real-time time point
                              # s is maximal time step and used for data record
                              t = 0
                              s = 1 
                              lambda_s = 0
                              lambda_sup = 0
                              #########
                              # running the simulation one time step at the time (i) needed for the stochastic implementation for mutation rate
                              #-------------------
                              while (t < ST*DT) {
                                
                                
                                # assigning and checking previous time point for each subpopulation
                                # I still use this within while loop because I'm not sure whether the calculation leads to negative number 
                                # is S NA?, set to 0
                                if(is.na(S_i)){
                                  S_i    <- 0
                                  
                                  
                                  
                                  # is the total number of S bacteria larger than 1? , if not set to 0  (protects form running into small number errors)
                                }else if(S_i>=1) {
                                  S_i  <-  S_i
                                } else{
                                  S_i    <- 0}
                                
                                
                                # is RA NA?, set to 0
                                if(is.na(RA_i)){
                                  RA_i    <- 0
                                  
                                  
                                  # is the total number of RA bacteria larger than 1? , if not set to 0  (protects form running into small number errors)
                                }else if(RA_i>=1) {
                                  RA_i    <-  RA_i
                                } else{
                                  RA_i    <- 0}
                                
                                
                                # is RB NA?, set to 0
                                if(is.na(RB_i)){
                                  RB_i    <- 0
                                  
                                  # is the total number of RB bacteria larger than 1? , if not set to 0 (protects form running into small number errors)
                                } else if(RB_i>=1) {
                                  RB_i  <-  RB_i
                                } else{
                                  RB_i    <- 0}
                                
                                # is RAB NA?, set to 0
                                if(is.na(RAB_i)){
                                  RAB_i    <- 0
                                  
                                  # is the total number of RAB bacteria larger than 1?, if not set to 0 (protects form running into small number errors)
                                } else if(RAB_i>=1) {
                                  RAB_i  <-  RAB_i
                                } else {
                                  RAB_i    <- 0}
                                
                               

                                
                                
                                # Initialize lambda_max
                                lambda_max <- c(
                                  S_i * (Gmax + 0.01) * rS_RA + 
                                    S_i * (Gmax + 0.01) * rS_RB + 
                                    S_i * (Gmax + 0.01) * rS_RB * rS_RA,           
                                  RA_i * (Gmax*FIT + 0.01) * rRA_S + 
                                    (Gmax*FIT + 0.01) * rRA_RAB,      
                                  RB_i * (Gmax*FIT + 0.01) * rRB_S + 
                                    (Gmax*FIT + 0.01) * rRB_RAB,     
                                  RAB_i * (Gmax*FIT^2 + 0.01) * rRAB_RB + 
                                    (Gmax*FIT^2 + 0.01) * rRAB_RA)
                                
                              
                                
                                
                                
                                # Calculate lambda_sup

                                lambda_sup = sum(lambda_max)
                                
                                # Generate three random numbers r1,r2 and r3
                                # r1 is used to assess tau
                                # r2 is used to assess the first reaction
                                # r3 is used to do accept-rejection
                                
                                r1 = runif(1,0,1)
                                r2 = runif(1,0,1)
                                r3 = runif(1,0,1)
                                
                                # Calculate the time point tau 
                                
                                tau_PD <- (1/lambda_sup)*log(1/r1)
                                
                                # Doing Gillespie algorithm for mutation and ODE for growth and death
                                # if tau < 1
                                if(tau_PD < 1){
                                  
                                  t = t + tau_PD
                                  
                                  
                                  # ODE initialization
                                  t_inits <- c(S = S_i, RA = RA_i, RB = RB_i, RAB = RAB_i, A = A_t, B = B_t)
                                  theta <- c(
                                    ke = ke,
                                    Gmax = Gmax,
                                    V = V,
                                    FIT = FIT,
                                    Bmax = Bmax,
                                    Gmin_A = Gmin_A,
                                    Gmin_B = Gmin_B,
                                    HILL_A = HILL_A,
                                    HILL_B = HILL_B,
                                    MIC_S = MIC_S,
                                    MIC_R = MIC_R
                                  )
                                  
                                  ev <- PK_list(tau_PD = tau_PD,
                                                   t = t, 
                                                   Dose = Dose, 
                                                   ST = ST, 
                                                   Tau = Tau, 
                                                   dose_reg = dose_reg, 
                                                   i_mod = i_mod) 
                                  
                                  # ODE simulation for bacterial growth and death
                                  PD_ode <- tail(as.data.frame(rxSolve(PD_mod, 
                                                                       params = theta, 
                                                                       events = ev, 
                                                                       inits = t_inits,
                                                                       maxsteps = 1000000000)), n = 1)
                                  
                                  # Update population size
                                  
                                  S_i <- PD_ode$S
                                  RA_i <- PD_ode$RA
                                  RB_i <- PD_ode$RB
                                  RAB_i <- PD_ode$RAB
                                  A_t = PD_ode$A
                                  B_t = PD_ode$B

                                  
                                  
                                  # Update gS in a new time point
                                  gS <- ifelse(((1-((S_i+RA_i+RB_i+RAB_i)/(10^Bmax*V)))*Gmax) > 0, 
                                               ((1-((S_i+RA_i+RB_i+RAB_i)/(10^Bmax*V)))*Gmax), 0)
                                  gRA <- gS * FIT
                                  gRB <- gS * FIT
                                  gRAB <- gS * FIT^2
                                  
                                  
                                  
                                  # Initialize lambda_i
                                  lambda_i <- c(                         
                                    S_i * (gS + 0.01) * rS_RA,           # S to RA mutation
                                    S_i * (gS + 0.01) * rS_RB,           # S to RB mutation  
                                    S_i * (gS + 0.01) * rS_RB * rS_RA,   # S to RAB mutation
                                    RA_i * (gRA + 0.01) * rRA_S,         # RA to S recover
                                    RA_i * (gRA + 0.01) * rRA_RAB,       # RA to RAB mutation                                                         # RA populRAion deRAh
                                    RB_i * (gRB + 0.01) * rRB_S,         # RA to S recover
                                    RB_i * (gRB + 0.01) * rRB_RAB,       # RA to RAB mutation                                                      # RAB populRAion deRAh
                                    RAB_i * (gRAB + 0.01) * rRAB_RB,     # RAB to RB recover
                                    RAB_i * (gRAB + 0.01) * rRAB_RA      # RAB to RA recover
                                  )
                                  
                                  # Calculate lambada_s
                                  
                                  lambda_s = sum(lambda_i) 
                                  
                                  if (r3 <= lambda_s/lambda_sup){
                                    
                                    alpha_sum = cumsum(lambda_i)
                                    reaction_first <- as.numeric((which(alpha_sum >= r2*lambda_s))[[1]])
                                    
                                    if(reaction_first == 1){RA_i = RA_i + 1; S_i = S_i - 1}
                                    else if(reaction_first == 2){RB_i = RB_i + 1; S_i = S_i - 1}
                                    else if(reaction_first == 3){RAB_i = RAB_i + 1; S_i = S_i - 1}
                                    else if(reaction_first == 4){S_i = S_i + 1; RA_i = RA_i - 1}
                                    else if(reaction_first == 5){RAB_i = RAB_i + 1; RA_i = RA_i - 1}
                                    else if(reaction_first == 6){S_i = S_i + 1; RB_i = RB_i - 1}
                                    else if(reaction_first == 7){RAB_i = RAB_i + 1; RB_i = RB_i - 1}
                                    else if(reaction_first == 8){RB_i = RB_i + 1; RAB_i = RAB_i - 1}
                                    else if(reaction_first == 9){RA_i = RA_i + 1; RAB_i = RAB_i - 1}
                                    
                                    
                                  }
                                  
                                  
                                }
                                
                                # Doing ODE for growth and death
                                # if tau >= 1
                                if(tau_PD >= 1) {
                                  
                                  tau_PD = 1
                                  t = t + tau_PD
                                  
                                  
                                  
                                  
                                  
                                  # ODE initialization
                                  t_inits <- c(S = S_i, RA = RA_i, RB = RB_i, RAB = RAB_i, A = A_t, B = B_t)
                                  theta <- c(
                                    ke = ke,
                                    Gmax = Gmax,
                                    V = V,
                                    FIT = FIT,
                                    Bmax = Bmax,
                                    Gmin_A = Gmin_A,
                                    Gmin_B = Gmin_B,
                                    HILL_A = HILL_A,
                                    HILL_B = HILL_B,
                                    MIC_S = MIC_S,
                                    MIC_R = MIC_R
                                  )
                                  
                                  ev <- PK_list(tau_PD = 1,
                                                t = t, 
                                                Dose = Dose, 
                                                ST = ST, 
                                                Tau = Tau, 
                                                dose_reg = dose_reg, 
                                                i_mod = i_mod) 
                                  
                                  PD_ode <- tail(as.data.frame(rxSolve(PD_mod, 
                                                                       params = theta, 
                                                                       events = ev, 
                                                                       inits = t_inits,
                                                                       maxsteps = 1000000000)), n = 1)
                                  
                                  # Update population size
                                  
                                  S_i <- PD_ode$S
                                  RA_i <- PD_ode$RA
                                  RB_i <- PD_ode$RB
                                  RAB_i <- PD_ode$RAB
                                  A_t = PD_ode$A
                                  B_t = PD_ode$B
                                }
                                
                                

                                #---------------------------------------------------------
                                #---------------------------------------------------------

                                # record data every 1 h
                                if(t >= s){
                                  x_CS <- data.frame(time = s, 
                                                     S = S_i/V, 
                                                     RA = RA_i/V, 
                                                     RB = RB_i/V, 
                                                     RAB = RAB_i/V,
                                                     A = A_t,
                                                     B = B_t)
                                  
                                  df_CS[(s + 1),] <- x_CS
                                  s = s + 1
                                  
                                }
                                
                              }
                              
                              
                              
                              
                              
                              # Simulation for one drug regimen complete
                              
                              # Add name of dosing regimen to data frame
                              df_CS$model <- dose_reg[i_mod] 
                              
                              # bind all simulation for specific realization ii 
                              
                              df_model_CS <- df_model_CS %>%
                                bind_rows(df_CS)
                            }
                            
                            
                            # Simulations for different drug regimens complete
                            
                            
                            
                            # add realization identifier (Total numbers of simulation)
                            df_model_CS$index <- ii
                            #return simulated data for all dosing regimen per realization ii
                            return(df_model_CS)
                            
                            
                          }
  # Parallel process complete
  
  
  # Output and calculations for all simulations 
  
  mean_dat <-  df_full_CS %>% 
    gather(value = "CFU", key = "Population", -time, -A, -B, -index, -model) %>% #make into long format
    group_by(time, model, index) %>%
    mutate(total = sum(CFU)) %>% #calculate  total CFU per realization, time point and dosing regimen
    ungroup() %>%
    mutate(CFU_ratio = CFU/total) %>%
    #reorder levels of bacterial subpopulations 
    mutate(Order = ifelse(Population == "S", 1,
                          ifelse(Population == "RA", 2,
                                 ifelse(Population == "RB", 3, 4)))) %>%
    mutate(Population = reorder(Population, Order)) %>% 
    ungroup() %>%
    group_by(model, index, time) %>% 
    mutate(CFU_total = sum(CFU, na.rm = T)) %>%  #calculate  total CFU per realization, time point and dosing regimen
    ungroup() %>% 
    group_by(time, model, Population) %>% 
    #calculate metrics relating to CFU per time point and dosing regimen
    mutate(CFU_MEDIAN = median(CFU, na.rm = T), 
           CFU_SD   = sd(CFU, na.rm = T), 
           CFU_95   = quantile(CFU, .95, na.rm = T),
           CFU_05   = quantile(CFU, .05, na.rm = T),
           mean_ratio = mean(CFU_ratio, na.rm = T),
           CFU_T_MEDIAN = median(CFU_total, na.rm = T),
           CFU_T_SD   = sd(CFU_total, na.rm = T), 
           CFU_T_95   = quantile(CFU_total, .95, na.rm = T),
           CFU_T_05   = quantile(CFU_total, .05, na.rm = T),
           #including specific model parameters
           HILL_A = HILL_A,
           HILL_B = HILL_B,
           GMIN_A = Gmin_A,
           GMIN_B = Gmin_B,
           eB0 = eB0,
           U_1 = u_1,
           U_2 = u_2,
           U_3 = u_3,
           U_4 = u_4,
           MIC_S = MIC_S,
           MIC_R = MIC_R) %>% 
    ungroup() %>% 
    select(time, A, B, Population, model, CFU_MEDIAN, CFU_SD, CFU_95, CFU_05, CFU, CFU_total,
           CFU_T_MEDIAN, CFU_T_SD, CFU_T_95, CFU_T_05, mean_ratio , index, 
           HILL_A, HILL_B, GMIN_A, GMIN_B, U_1, U_2, U_3, U_4, eB0,
           MIC_S, MIC_R) %>% 
    filter(!is.na(model))
  
  run_t   <-  Sys.time() - start_t 
  print(run_t)
  
  return(mean_dat)
  
  
  
}
