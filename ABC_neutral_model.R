### Neutral model
### rates of cell division and turnover are constant
### so is the per capita rate of influx of cells into the compartment

## clearing the workspace and the memory
rm(list = ls()); gc()

library(tidyverse)
library(deSolve)
library(parallel)
library(mvtnorm)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

print(args)
num_cores =  detectCores() - 4

ModelName <- "neutral"
Population <- 'CD4' #args[1]
DataFile <- file.path("Datafiles", paste0(Population, "_ln.csv"))
OutputDir <- file.path("output", Population, ModelName)
Beta = 1/3.5

print(Population)


### importing data to be fitted 
#data_fit <- read.csv(DataFile) %>%
#  mutate(time_post_t0 = time - min(time)) %>% arrange(time_post_t0)

data_fit <- read.csv(DataFile) %>%
  mutate(age.at.S1K = time,
         time_post_t0 = age.at.S1K - min(age.at.S1K)) %>% arrange(age.at.S1K) %>% unique() %>%
  select("age.at.S1K", contains('counts'), contains('ki67'))


# unique time points in data for odes solver
unique_times <- data_fit %>% distinct(age.at.S1K, .keep_all = TRUE)

data_time <- data_fit$age.at.S1K                        # timepoints in observations 
solve_time <- unique_times$age.at.S1K                   #unique time points to solve odes  

time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time

#### data transformation functions
logit_trans <- function(x){
  
  log(x/(1-x))
}

expit_trans <- function(x){
  
  exp(x)/(1 + exp(x))
}

#####################################################################################################################################
#####################################################################################################################################

### model definition
## initial guesses of the parametrs to start the optimiser
start_vec <- c(psi = 0.4, delta = 0.05, rho = 0.002)#, N0_log = 14,kappa0_logit = 1.4)

### Estimating the total pool size of the compartment and the proportions of ki67+ cells

## function that defines influx of cell from the source compartment (SP cells in thymus)
theta_spline <- function(time, parms){
  psi = parms[1]  ## psi is the per capita rate constant of daily influx of SP4 cell
  
  t0 = 5    ### time zero in data
  
  ## parameters estimated from spline fit to the timecourse of counts of source compartment -- SP
  if (Population == "CD4") {
    spline_est_list <- list("basl"  = 4.3e5, 'theta' = 1.8e3, "n" = 2.1, "X" = 30, "q" = 3.7)
  } else if (Population == "CD8") {
    spline_est_list <- list("basl"  = 9e4, 'theta' = 68, "n" = 3, "X" = 25, "q" = 4.25)
  }
  
  ## spline for SP4 cell counts ~ time
  fit1 = spline_est_list$basl + (spline_est_list$theta * time^spline_est_list$n) *
         (1 - ((time^spline_est_list$q)/((spline_est_list$X^spline_est_list$q) + (time^spline_est_list$q))))
  
 
  ## output of this function is the number of cells entring naive CD4 T pool daily
  value = psi * fit1
  
  return(value)
}

## function to depict the changes in the proportions ki67+ cells within the precursor (SP4) popualtion
kiSource_spline <- function(time) {
  ## spline fitted separately to the proportions of ki67+ cells within thymic SP4 cells
  
  ## parameters estimated from spline fit to the timecourse of counts of 
  ## source compartment -- SP CD4
  if (Population == "CD4") {
    fit2 = exp(-0.03337899 * (time + 2.92554110)) + 0.13732103     ## best fitting spline for SP4 ki67 proportions
  } else if (Population == "CD8"){
    fit2 = exp(-0.01559996 * (time + 14.83715328)) + 0.24510453        ## best fitting spline for SP8 ki67 proportions
  }
  
  return (fit2)
}

#### fraction of ki67+ cells
ode_solver <- function(time_vec, parms_vec, return_obj){
  
  #fVRTE_logit  = parms_vec[4]  ## fraction of VRTE cells in ki67 high subset
  #N0_log       = parms[5]
  #kappa0_logit = parms[6]
  
  ## ode function
  ode_func <-  function(t, state, parms){
    with(as.list(c(state, parms)),{
      
      delta = parms[2]
      rho = parms[3]
      
      #Ki67 hi VRTE
      #dY1 <- theta_spline(t, parms) * kiSource_spline(t) - 
      #  (Beta_prime + rho + delta) * Y1
      
      #Ki67 hi mN
      dY1 <- theta_spline(t, parms) * kiSource_spline(t) + rho * (2 * Y2 + Y1) - (Beta + delta) * Y1
      
      #Ki67 lo mN
      dY2 <- theta_spline(t, parms) * (1 - kiSource_spline(t)) + Beta * Y1 - (rho + delta) * Y2
      
      #return the rate of change
      list(c(dY1, dY2))
      
    })  # end with(as.list ...
  }
  
  N0 = mean(data_fit$counts[1:2]) #exp(N0_log)  -- counts at t0
  kappa0 = mean(data_fit$ki67[1:2]) #expit_trans(kappa0_logit) -- ki67hi fraction at t0a
  #fVRTE = expit_trans(fVRTE_logit)
  
  #initial conditions
  state <- c(N0 * kappa0,
             N0 * (1 - kappa0))
  names(state) <- c("Y1", "Y2")
  
  #state <- c(Y1 = exp(14) * 0.5 * 0.8, Y2 = exp(14) * 0.8 * 0.5, Y3 = exp(14) * 0.2)
  
  # time points for which conc is reported
  # include the points where data is available
  times = time_vec
  
  sol_ode <- ode(y= state, times=times, func = ode_func, parms = parms_vec)
  
  sol_df <-  data.frame(sol_ode) %>%
    mutate(total_counts = Y1 + Y2,
           ki_prop = (Y1)/ total_counts) #%>% na.omit
  
  if (return_obj == "counts"){
    return(sol_df$total_counts)
  } else {
  return(sol_df$ki_prop)
  }
}

print(ode_solver(solve_time, start_vec, "counts"))
print(ode_solver(solve_time, start_vec, "ki67"))


##################################################################################################

### ABC sampler
 sample_one_from_pior <- function(Time){
   psi <- runif(1, 0, 1)
   delta <- rnorm(1, 0.05, 0.015)
   rho <- runif(1, 0, 1)
   
   theta = c(psi, delta, rho)
   
   res1 = ode_solver(Time, theta, 'counts')
   res2 = ode_solver(Time, theta, 'ki')
   
   return(list('res1' = res1, 'res2' = res2, 'theta' = theta))
 }

 ### define the distance between simulation and data
 dist_func <- function(sim_data, boot_data){
   CC = boot_data$counts
   KI = boot_data$ki67
   tdata = boot_data$age.at.S1K
   
   sim1 = sim_data$sim1[time_index]
   sim2 = sim_data$sim2[time_index]
   
   resi_cc = log(CC) - log(sim1)
   resi_ki67 = logit_trans(KI) - logit_trans(sim2)
   res = sum(resi_cc ** 2) + sum(resi_ki67 ** 2)
   
   return(res)
 }
 
 # Sample 1 theta from prior defined above and accept if |simulation - data| < pre-defined distance
 rej_samp_popt0 <- function(eps){
     # Sample theta from pi(theta) for population t=0, particle n=i
   sim <- sample_one_from_pior(solve_time)
   sim1 = sim$res1
   sim2 = sim$res2
   theta =  sim$theta
   
   sim_data <- data.frame(sim1, sim2)
   dis <- dist_func(sim_data, data_fit)
   #print(dis)
   
   while (dis > eps) {
     sim <- sample_one_from_pior(solve_time)
     sim1 = sim$res1
     sim2 = sim$res2
     theta =  sim$theta
     
     sim_data <- data.frame(sim1, sim2)
     dis <- dist_func(sim_data, data_fit)
     print(dis)
   }
   
   return(theta)
 }
 
## function that feeds in mapply for parallel computing
rej_samp_popt0_call <- function(i, eps){
  #Sample N particles from pi(theta)
  theta_n <- data.frame()
  theta_n[i, 1:3] <- rej_samp_popt0(eps)
  return(theta_n[i, 1:3])
}

rej_samp_popt0_call(10, 80)
 
rej_samp_Npart_popt0 <- function(n, eps){
    #Sample N particles from pi(theta)
  theta_n <- data.frame(mcmapply(rej_samp_popt0_call, i = seq(1, n, 1), eps= eps, mc.cores = num_cores),
                        row.names = NULL)# c("psi", "deltalog", "rholog"))
  #print(theta_n)
  return(t(data.matrix(theta_n, rownames.force = NA)))
}

theta_n_particles <- rej_samp_Npart_popt0(30, 50)
prob_wghts <-rep(1/nrow(theta_n_particles), nrow(theta_n_particles))

# Sample 1 theta from previous generation of posterior
sample_one_from_theta <- function(theta_n_particles, prob_wghts){
  j = sample(seq(1, length(prob_wghts)), 1, p = prob_wghts)
  theta = theta_n_particles[j, ]
  return(theta)
}

# Check if perturb theta is still within posterior bounds
# used also to calculate new weights
joint_pi <- function(theta, abound, bbound){
  pi = prod(1/(bbound - abound)) ## aboud and bbound are arrays
  for (i in 1:length(theta)){
    if (theta[i] < abound[i] | theta[i]>bbound[i]){
      pi = pi*0}
  }
  return(pi)
}

# used to calculate new weights with gaussian kernel
kernel_func <- function(theta_t, theta_tm1){
  cov_theta = 2 * cov(theta_tm1)
  pi_kernel = array()
  for (i in 1:nrow(theta_tm1)){
    pi_kernel[i] = dmvnorm(theta_t, mean = theta_tm1[i, ], sigma = cov_theta)
  }
  return(pi_kernel)
}

move_one_particle <- function(theta_prev_pop, prob_wghts, eps){
  #Move particle i for population t from theta(t-1) = theta_prev_pop
  
  # move particle with kernel
  # UPDATE BOUNDS FOR PARAMETERS HERE IF THEY ARE DIFFERENT
  abound = c(0, 0, 0)
  bbound = c(1, 1, 1)
  cov_theta = 2 * cov(theta_prev_pop)
  
  dis = eps+1
  pi_newtheta = 0
  NoOfSamp = 0
  while(dis > eps | pi_newtheta == 0){
    #Sample theta from previous theta with weights
    theta = sample_one_from_theta(theta_prev_pop, prob_wghts)
    # move particle with kernel
    newtheta = rmvnorm(1, mean = theta, sigma = cov_theta)
    # Check if you are in pi(theta)
    pi_newtheta = joint_pi(newtheta, abound, bbound)
    
    # Check the distance
    res1 <- ode_solver(solve_time, newtheta, 'counts')
    res2 <- ode_solver(solve_time, newtheta, 'ki67')
    
    sim_data <- data.frame('sim1' = res1, "sim2" = res2)
    dis =dist_func(sim_data, data_fit)
    NoOfSamp = NoOfSamp +1
  }
  NewWeights = joint_pi(newtheta, abound, bbound)/sum(prob_wghts/kernel_func(newtheta, theta_prev_pop))
  return(c("newtheta" = newtheta, 'NewWeights' = NewWeights, 'dis' = dis, 'NoOfSamp' = NoOfSamp))
}

move_one_particle_call <- function(i, theta_prev_pop, prob_wghts, eps){
  if (i %% 10 == 0) {
    print(paste0('particle: ', i))
  }
  res_n <- data.frame()
  res_n[i, 1:6] <-  move_one_particle(theta_prev_pop, prob_wghts, eps)
  #return(res_n[i, 1:6])
  return(c("theta" = res_n[i, 1:3], "weights" = res_n[i, 4], 'dis' = res_n[i, 5], 'NoOfSamp' = res_n[i, 6]))
}

# Move N particles
move_N_particles <- function(N, theta_prev, prob_wghts, eps){
    res <- pbmcmapply(move_one_particle_call, i = seq(1, N), MoreArgs = list(theta_prev, prob_wghts, eps), mc.cores = num_cores)
    return(data.frame(t(res), row.names = NULL))
}


abc_smc <- function(eps_max, eps_min, N){
  # eps = eps[0], ..., eps [T]  - T no of populations
  # N - no of particles
  
  t = 0
  eps = eps_max
  while (eps > eps_min) {
    print(paste0('Population ', t, ', epsilon: ', eps))
    if (t == 0) {
      # Sample N particles from pi(theta)
      thetaSS = rej_samp_Npart_popt0(N, eps)
      # Assign weights 1/N
      weightS = rep(1/N, N)
      # Theta(t = 0) = [thetaSS]
      thetaS = thetaSS
      dis = 0
      NoOfSamp = 0
    } else {
      res <- move_N_particles(N, thetaS, weightS, eps)
      weightSS <-  data.matrix(res[, 'weights'])
      dis <-  data.matrix(res[, 'dis'])
      NoOfSamp <-  data.matrix(res[, 'NoOfSamp'])
      weightS <- data.matrix(as.numeric(weightSS)/sum(as.numeric(weightSS)))
      thetaS <- data.matrix(res[, 1:3])
      
      # write.csv(res, paste0('post_neutral_eps_', eps, '.csv'))
      eps = median(as.numeric(dis))
    }
    t = t+1
  }
  write.csv(res, paste0('post_neutral_final.csv'))
  return(t)
}

abc_smc(eps_max = 50, eps_min = 10, N = 100)
