
# PK-PD model for INC, SAME and SEQ treatments based on adaptive tau algorithm 
# Contact: Junqi Liao (jliao@zoologie.uni-kiel.de) 

# The general Model structure is taken from Aulin et al. 2021 as well as the PK Model (with slight changes) 
# all code directly used from Aulin et al 2021 and changes made to it are marked accordingly 
# For Information on the RxODE package check Wang et al. 2016
# For information on the adaptive tau-leaping algorithm check *please include*
# Refer to Linda B. S. Aulin, et al. 2021 (General structure and PK)
# Refer to W Wang, et al. 2015 (PK)
# Refer to Cao Y, et al. 2007 (PD)

# Create variables, actual values were inputted in running script 
# Model function as set up by Aulin et al., Variables for CS and CR were deleted
# we included variables for different MIC values of susceptible and resistant strains and parameter Gmax
# all parameter are set to default values, which can be overridden by the function input

CS_model <- function(t_half = 2, F_Css_MIC = 1, v_Models = c("Mono A", 
                                                             "Mono B", 
                                                             "3 day cycling",  
                                                             "1 day cycling", 
                                                             "Multi-step evolution"), 
                     FIT = 1, Gmin_A = -6, Gmin_B = -6, HILL_A  = 1, HILL_B = 1, u_1 = 10^-6, u_2 = 10^-6, eB0 = 6, RA0=0, RB0 =0 , RAB0 = 0,  Bmax = 9,  V = 1000, n = 100, DT = 0.1, ST = 24, MIC_S = 1, MIC_R =8, Gmax = 1.5)
  
  # We included following treatments from Aulin et al. 2021:
  # Mono: Single drug with constant concentrations (SAME evolution in main text)
  # 3 day cycling: Switching drug every 3 days
  # 1 day cycling: Switching drug every 1 day
  # We extended the list of treatments for INC evolution: 
  # Multi-step evolution: Single drug with increasing concentrations (INC evolution in main text)
  
  
  # argument explanation
  # Explanation of input parameters as given by Aulin et al. 2021 
# t_half; half life used for both drugs [h]
# F_Css_MIC scaling; factor to obtain Css average equal to X time MIC WT
# FIT; Fitness cost (0-1), 1 no cost, 0 no growth 
# Gmin_A; Maximal effect of drug A
# Gmin_B; Maximal effect of drug B
# HILL_A; Shape factor drug A
# HILL_B; Shape factor drug B
# u_1; Mutation rate of S against antibiotic A
# u_2; Mutation rate of S against antibiotic B
# eB0; Starting bacterial concentration [cfu/mL]
# RA0; Fraction of bacterial starting in RA
# RB0; Fraction of bacterial starting in RB
# RAB0; Fraction of bacterial starting in RAB
# Bmax; carrying capacity of the system [cfu/mL]
# V; volume of infection site
# n; number of simulations
# DT; time step used [h]
# ST; Duration of simulation [h]

# Additionally included Parameter 
# MIC_S; Minimum inhibitory concentration of the susceptible type
# MIC_R; Minimum inhibitory concentration of the resistant type  
# Gmax = maximum growth rate in the absence of antibiotics 


{
  # Check if all reqiuered packages are present 
  require(doParallel)
  require(doRNG)
  require(rxode2)
  require(dplyr)
  require(tidyr)
  require(adaptivetau)
  
  # Aulin et al. 2021: Function used for integer check in output calculation 
  # (Such as when we want to calculate the results only by day instead of hour)
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  # Used for integer check in output calculation 
  # (Such as when we want to calculate the results only by day instead of hour)
  
  
  
  
  start_t <- Sys.time()  #for checking run time
  
  ##############
  ## PK Model ##
  ##############
  
  #Define PK parameters and dosing, assume same PK for both drugs as given by Aulin et al. 2021
  # Our model excluded the Volume of distribution and assumes Dose to be a concentration 
  
  Vd = 1                  # [L] volume of distribution, based on plasma volume
  ke = log(2)/t_half      # [h^-1] elimination rate constant, calculated by half-life
  CL = ke*Vd              # [L/h] clearance 
  Tau = 12                # [h] dosing interval
  MIC = 1                 # [ug/ul] default MIC value
  
  Css = MIC*F_Css_MIC     # [mg/L] Target average steady-state concentration 
  
  Dose = Css*ke*Tau       # [mg/L] dose giving target Css
  
  
  
  # PK model for antibiotic A and B and definitions of treatments as given by Aulin et al 2021
  # ODE
  
  
  PK_mod<- rxode2({
    
    
    d/dt(A) = -ke*A   
    d/dt(B) = -ke*B
    
  });
  
  
  PK_inits <- c(A = 0, B = 0);
  
  PK_theta <- c(ke = ke)
  
  # Dosing regimens
  # specified using argument v_Models
  
  
  PK_model_list <- list()
  
  if( "Mono A" %in% v_Models){  
    
    ev_PK_A <- eventTable(amount.units="mg/L", time.units="hours") %>%
      add.dosing(dosing.to = 1, dose=Dose, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST, DT))  %>%
      as.tbl()
    
    PK_A      <- as.data.frame(PK_mod$run(PK_theta, ev_PK_A,     PK_inits)) %>% 
      mutate(Model= "Mono A")
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_A  }
  
  
  if( "Mono B" %in% v_Models){
    
    ev_PK_B <- eventTable(amount.units="mg/L", time.units="hours") %>%
      add.dosing(dosing.to = 2, dose=Dose, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST, DT))  %>%
      as.tbl()
    
    PK_B      <- as.data.frame(PK_mod$run(PK_theta, ev_PK_B, PK_inits)) %>% 
      mutate(Model = "Mono B") 
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_B }
  
  
  if( "3 day cycling" %in% v_Models){
    
    
    ev_PK_3day <- eventTable(amount.units="mg/L", time.units="hours") %>%
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
    
    
    ev_PK_1day <- eventTable(amount.units="mg/L", time.units="hours") %>%   #2 dose cycling
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
  # Additional treatment added 
  if( "Multi-step evolution" %in% v_Models){
    
    
    ev_PK_multi <- eventTable(amount.units="mg/L", time.units="hours") %>%   #Multi-step
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
  
  # The following code is not from the initial function implemented by Aulin et al. 2021
  # We used the adaptive tau leaping algorithm instead and incorporated a different PD Model 
  ##############
  ## PD Model ##
  ##############
  
  # 4 state model including S (= WT), RA, RB, RAB)
  # SSA Adaptive tau algorithm for PD 
  # The same order as mrate function following
  # List of all transitions (i.e. events) appearing in the simulation:
  transitions = list(c(S = +1),
                     c(S = -1),
                     c(S = -1, RA = +1),
                     c(S = -1, RB = +1),
                     c(S = -1, RAB = +1),
                     c(RA = +1),
                     c(RA = -1),
                     c(RA = -1, RAB = +1),
                     c(RB = +1),
                     c(RB = -1),
                     c(RB = -1, RAB = +1),
                     c(RAB = +1),
                     c(RAB = -1)
  )
  
  
  # Reactions for different mutation process, the same order as the previous transition list
  # Avoid replication rate less then 0 (i.e. when population size foes slightly above b_max)
  # Intrinsic death rate is set 0.01 and added to Gmax to get the replication rate r0 (needed to calculate
  # the number of mutants)
  mrates <- function(x, params, t) {
    S <- x["S"]
    RA <- x["RA"]
    RB <- x["RB"]
    RAB <- x["RAB"]
    
    return(c(
      S * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax) > 0, 
                  ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax), 0)),                                                      # S population growth (Including intrinsic death)
      S * (params$dS),                                                                                              # S population death (Drug effect)
      S * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)) > 0, 
                  ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)), 0)) * params$rS_RA * (1-params$rS_RB),      # S to A mutation
      S * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)) > 0, 
                  ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)), 0)) * params$rS_RB * (1-params$rS_RA),      # S to B mutation
      S * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)) > 0, 
                  ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)), 0)) * params$rS_RA * params$rS_RB,          # S to AB mutation
      RA * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax*FIT) > 0, 
                   ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax*FIT), 0)),                                                 # RA population growth
      RA * (params$dRA),                                                                                            # RA population death
      RA * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)*FIT) > 0, 
                   ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)*FIT), 0)) * params$rRA_RAB,                  # RA to RAB mutation
      RB * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax*FIT) > 0, 
                   ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax*FIT), 0)),                                                 # RB population growth
      RB * (params$dRB),                                                                                            # RB population death
      RB * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)*FIT) > 0, 
                   ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*(Gmax+0.01)*FIT), 0)) * params$rRB_RAB,                  # RB to RAB mutation
      RAB * (ifelse(((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax*FIT) > 0, 
                    ((1-((S+RA+RB+RAB)/(10^Bmax*V)))*Gmax*FIT), 0)),                                                # RAB population growth
      RAB * (params$dRAB)                                                                                           # RAB population death
      
    ))
  }
  
  # Set up bacterial growth and death as deterministic step (not stochastic activity)
  deterministic <- c(TRUE, TRUE, FALSE, FALSE, FALSE,
                     TRUE, TRUE, FALSE, FALSE,
                     TRUE, TRUE, FALSE, FALSE,
                     TRUE, TRUE, FALSE, FALSE)
  
  
  
  
  S0  <- (10^eB0)*(1-RA0 - RB0 - RAB0)  # starting bacterial density [cfu/ml] of sensitive bacteria (S) = WT)
  
  U_1 <- u_1  # Mutation rate [mutation/bacteria/hour] against antibiotic A
  U_2 <- u_2  # Mutation rate [mutation/bacteria/hour] against antibiotic B
  
  # Parallelization for running the different treatments as implemented by Aulin et al. 2021 
  # Structure of the data files analogue to Aulin et al. 2021
  n_obs  <- ST*DT+1   # number of observations ( simulation time(ST) * time step size (DT))
  n_models <- length(PK_model_list) #which models (i.e. dosing regimens) to simulate
  
  df_full_CS <-  foreach( ii = 1:n, .combine = "rbind",              # Parallelization
                          .packages = c("rxode2", "dplyr", "tidyr", "adaptivetau"),
                          .inorder = F,
                          # Only when the export order is not important
                          .options.RNG = 123) %dorng%  {            
                            # Setting seed to make sure we get identical values every time (pseudo random number)
                            
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
                            
                            for(i_mod in 1:length(PK_model_list)) {
                              #   #   
                              PK_mod = PK_model_list[[i_mod]]
                              # Creating empty data frame to  be populated with simulated data per dosing regimen i_mode for realization ii
                              
                              df_CS <- data.frame(time = rep(x = NA, times = n_obs),
                                                  S    = rep(x = NA, times = n_obs),
                                                  RA   = rep(x = NA, times = n_obs),
                                                  RB   = rep(x = NA, times = n_obs),
                                                  RAB  = rep(x = NA, times = n_obs),
                                                  A    = rep(x = NA, times = n_obs),
                                                  B    = rep(x = NA, times = n_obs))
                              
                              
                              
                              
                              # Set up starting data (population t = 0)
                              df_CS[1,] <- data.frame(time = 0, S = S0, RA = RA0*S0, RB = RB0*S0, RAB = RAB0*S0,
                                                      A = PK_mod$A[1], B = PK_mod$B[1]);
                              # Use absolute number of bacteria  
                              V = V
                              S_t = S0 * V
                              RA_t = RA0 * V
                              RB_t = RB0 * V
                              RAB_t = RAB0 * V
                              
                              #########
                              # running the simulation one time step at the time (i) needed for the stochastic implementation for mutation rate
                              # Run adaptive tau for discrete time steps of length 0.1 h: death rates are time dependent, which adaptive tau does account for i.e. every time step of 0.1, we update the death rates                              
                              # Time interval is 0.1 h based on comparison with Gillespie algorithm
                              
                              for (i in 1:(ST*10)) {
                                
                                
                                # assigning and checking previous time point for each subpopulation
                                # analogue to Aulin et al. 2021: check if the number of bacteria can be converted to integer (needed for random sampling)
                                # is S NA?, set to 0
                                if(is.na(S_t)){
                                  S_t    <- 0
                                  
                                  
                                  
                                  # is the total number of S bacteria larger than 1? , if not set to 0  (protects form running into small number errors)
                                }else if(S_t>=1) {
                                  S_t  <-  S_t
                                } else{
                                  S_t    <- 0}
                                
                                
                                # is RA NA?, set to 0
                                if(is.na(RA_t)){
                                  RA_t    <- 0
                                  
                                  
                                  # is the total number of RA bacteria larger than 1? , if not set to 0  (protects form running into small number errors)
                                }else if(RA_t>=1) {
                                  RA_t    <-  RA_t
                                } else{
                                  RA_t    <- 0}
                                
                                
                                # is RB NA?, set to 0
                                if(is.na(RB_t)){
                                  RB_t    <- 0
                                  
                                  # is the total number of RB bacteria larger than 1? , if not set to 0 (protects form running into small number errors)
                                } else if(RB_t>=1) {
                                  RB_t  <-  RB_t
                                } else{
                                  RB_t    <- 0}
                                
                                # is RAB NA?, set to 0
                                if(is.na(RAB_t)){
                                  RAB_t    <- 0
                                  
                                  
                                  
                                  # is the total number of RAB bacteria larger than 1?, if not set to 0 (protects form running into small number errors)
                                } else if(RAB_t>=1) {
                                  RAB_t  <-  RAB_t
                                } else {
                                  RAB_t    <- 0}
                                
                                # assign drug concentration from PK
                                A_t  <- PK_mod$A[i+1]
                                B_t  <- PK_mod$B[i+1]
                                
                                
                                #-------------------
                                # Stochastic implementation
                                x0 = c(S = floor(S_t), 
                                       RA = floor(RA_t),
                                       RB = floor(RB_t),
                                       RAB = floor(RAB_t))
                                
                                # Input parameters for growth rate and mutation rate in specific time 
                                # Updated for every time step 
                                params = list(
                                  
                                  dS = (((1 - Gmin_A/Gmax)*(A_t/MIC_S)^HILL_A/((A_t/MIC_S)^HILL_A - (Gmin_A/Gmax)))   +
                                          ((1 - Gmin_B/Gmax)*(B_t/MIC_S)^HILL_B/((B_t/MIC_S)^HILL_B - (Gmin_B/Gmax))))*Gmax,
                                  dRA = (((1 - Gmin_A/Gmax)*(A_t/MIC_R)^HILL_A/((A_t/MIC_R)^HILL_A - (Gmin_A/Gmax)))   +
                                           ((1 - Gmin_B/Gmax)*(B_t/MIC_S)^HILL_B/((B_t/MIC_S)^HILL_B - (Gmin_B/Gmax))))*Gmax*FIT,
                                  dRB = (((1 - Gmin_A/Gmax)*(A_t/MIC_S)^HILL_A/((A_t/MIC_S)^HILL_A - (Gmin_A/Gmax)))   +
                                           ((1 - Gmin_B/Gmax)*(B_t/MIC_R)^HILL_B/((B_t/MIC_R)^HILL_B - (Gmin_B/Gmax))))*Gmax*FIT,
                                  dRAB = (((1 - Gmin_A/Gmax)*(A_t/MIC_R)^HILL_A/((A_t/MIC_R)^HILL_A - (Gmin_A/Gmax)))   +
                                            ((1 - Gmin_B/Gmax)*(B_t/MIC_R)^HILL_B/((B_t/MIC_R)^HILL_B - (Gmin_B/Gmax))))*Gmax*FIT*FIT,
                                  
                                  rS_RA = u_1,
                                  rS_RB = u_2,
                                  rRA_RAB = u_2,
                                  rRB_RAB = u_1)
                                
                                # ssa.adaptive() was used to do stochastic implementation with Implicit-Explicit tau-leap method
                                r = ssa.adaptivetau(x0, transitions, mrates, params, deterministic = deterministic, tf = 0.1, maxtau = 0.1)
                                
                                # Input calculation results to data frame 
                                PD_mod <- tail(as.data.frame(r), n = 1)
                                S_t = PD_mod$S
                                RA_t = PD_mod$RA
                                RB_t = PD_mod$RB
                                RAB_t = PD_mod$RAB
                                
                                # Only collect data every hour 
                                if(is.wholenumber(i/10)){
                                  x_CS <- PD_mod %>%
                                    mutate(time = i/10,
                                           S = S/V,
                                           RA = RA/V,
                                           RB = RB/V,
                                           RAB = RAB/V,
                                           A = PK_mod$A[i+1],
                                           B = PK_mod$B[i+1])
                                  
                                  
                                  # Input x_CS data to two rows in df_CS
                                  df_CS[(i/10+1),] <- x_CS
                                }
                                
                                
                                
                                
                                
                                
                                
                                
                                
                              }
                              # Data is saved in same data structure as given in Aulin et al 2021
                              # Simulation for one drug regimen complete
                              
                              # Add name of dosing regimen to data frame
                              df_CS$model <- dose_reg[i_mod] 
                              
                              # bind all simulation for specific realization ii 
                              
                              df_model_CS <- df_model_CS %>%
                                bind_rows(df_CS)
                              
                            }
                            # Simulations for different drug regimens complate
                            
                            
                            
                            # add realization identifier (Total numbers of simulation)
                            df_model_CS$index <- ii
                            
                            
                            #return simulated data for all dosing regimen per realization ii
                            return(df_model_CS)
                            
                          }
  # Parallel process complete
  
  
  # Output and calculations for all simulations as implemented by Aulin et al. 2021 
  
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
           MIC_S = MIC_S,
           MIC_R = MIC_R) %>% 
    ungroup() %>% 
    select(time, A, B, Population, model, CFU_MEDIAN, CFU_SD, CFU_95, CFU_05, CFU, CFU_total,
           CFU_T_MEDIAN, CFU_T_SD, CFU_T_95, CFU_T_05, mean_ratio , index, 
           HILL_A, HILL_B, GMIN_A, GMIN_B, U_1, U_2, eB0,
           MIC_S, MIC_R) %>% 
    filter(!is.na(model))
  
  run_t   <-  Sys.time() - start_t 
  print(run_t)
  
  return(mean_dat)
  
  
  
}