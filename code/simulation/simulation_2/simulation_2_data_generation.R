###Load packages

library(tidyverse)
library(rstan)
library(ggpubr)
library(shinystan)
library(cowplot)
library(devtools)
library(psych)

h_function_4para <- function(x, beta_1, beta_2, beta_3, beta_4) {
  beta_1 + beta_2/(1 + (x/beta_3)^(-beta_4))
}

# Number of simulation replicates
# For demonstration purposes, we provide code for a single replicate.
# In our full simulation study (reported in the manuscript), we used 500 replicates.
# To reproduce the full simulation, change this value to 500.
number_of_iteration <- 1

beta_1 <- 4
beta_2 <- 18000
beta_3 <- 8
beta_4 <- 1.5
sigma_y <- 0.1
standard_initial_concentration <- 125
###15 uncontaminated samples, 5 contaminated samples
number_of_unknown_samples <- 20
estimate_our_method_more_contamination <- matrix(0, nrow = number_of_iteration, ncol = number_of_unknown_samples)
estimate_2004_method_more_contamination <- matrix(0, nrow = number_of_iteration, ncol = number_of_unknown_samples)
for (m in 1:number_of_iteration){
  set.seed(m*10) # Need to change to random seed
  #Simulate Standards data
  d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                 0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048)
  y_standard <- c()
  for (i in 1:length(d_standard)){
    y_standard[i] <- rnorm(1, mean = log(h_function_4para(standard_initial_concentration*d_standard[i],
                                                          beta_1, 
                                                          beta_2,
                                                          beta_3,
                                                          beta_4)), sd = sigma_y)
  }
  #Simulate Unknown data 
  d_unknown = c(1/10, 1/100, 1/10000)
  unknown_concentration <- c(0, 0, 0, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16)
  y_unknown <- matrix(NA, nrow = length(unknown_concentration), ncol = length(d_unknown))
  for (i in 1:length(unknown_concentration)){
    for(j in 1:length(d_unknown)){
      y_unknown[i,j] <- rnorm(1, mean = log(h_function_4para(unknown_concentration[i]*d_unknown[j],
                                                             beta_1, 
                                                             beta_2,
                                                             beta_3,
                                                             beta_4)), sd = sigma_y)
    }
  }
  unknown_concentration_contaminated <- c(0, 1/16, 1/4, 2, 8)
  number_of_contaminated_samples <- 5
  delta_sigma_contaminated_more <- 1
  delta_beta_contaminated_more <- 1
  y_unknown_contaminated_more <- matrix(NA, nrow = length(unknown_concentration_contaminated), ncol = length(d_unknown))
  for (i in 1:number_of_contaminated_samples){
    for(j in 1:length(d_unknown)){
      y_unknown_contaminated_more[i,j] <- rnorm(1, mean = log(h_function_4para(unknown_concentration_contaminated[i]*d_unknown[j],
                                                                               beta_1, 
                                                                               beta_2*exp(delta_beta_contaminated_more),
                                                                               beta_3,
                                                                               beta_4)), sd = sigma_y*exp(delta_sigma_contaminated_more))
    }
  }
  ###severe contamination data
  combined_data_more_comtamination <- data.frame(
    id <- c(rep("standard", 26), rep("unknown_1", 3), rep("unknown_2", 3),rep("unknown_3", 3), rep("unknown_4", 3), rep("unknown_5", 3),
            rep("unknown_6", 3), rep("unknown_7", 3),rep("unknown_8", 3), rep("unknown_9", 3), rep("unknown_10", 3),
            rep("unknown_11", 3), rep("unknown_12", 3),rep("unknown_13", 3), rep("unknown_14", 3), rep("unknown_15", 3),
            rep("unknown_1_contaminated", 3), rep("unknown_2_contaminated", 3),rep("unknown_3_contaminated", 3), rep("unknown_4_contaminated", 3), rep("unknown_5_contaminated", 3)),
    dilution <- c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                  0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                  rep(c(1/10, 1/100, 1/10000), number_of_unknown_samples)),
    y <- c(y_standard, y_unknown[1,],y_unknown[2,],y_unknown[3,],y_unknown[4,],y_unknown[5,],
           y_unknown[6,],y_unknown[7,],y_unknown[8,],y_unknown[9,],y_unknown[10,],
           y_unknown[11,],y_unknown[12,],y_unknown[13,],y_unknown[14,],y_unknown[15,],
           y_unknown_contaminated_more[1,],y_unknown_contaminated_more[2,],
           y_unknown_contaminated_more[3,],y_unknown_contaminated_more[4,],
           y_unknown_contaminated_more[5,])
  )
  colnames(combined_data_more_comtamination) <- c("id", "dilution", "y")
  write.csv(combined_data_more_comtamination, file = paste("combined_data_more_comtamination", m, ".csv", sep=""))
  
  ###start simulation
  simulated_data_more <- list(
    N_stand_sample = 1,
    N_unknown_sample = number_of_unknown_samples,
    N_standard = 1*26,
    N_unknown = number_of_unknown_samples*3,
    y_standard = y_standard,
    d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                   0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048),
    ind_standard = c(rep(1,26)),
    y_unknown = c(y_unknown[1,],y_unknown[2,],y_unknown[3,],y_unknown[4,],y_unknown[5,],
                  y_unknown[6,],y_unknown[7,],y_unknown[8,],y_unknown[9,],y_unknown[10,],
                  y_unknown[11,],y_unknown[12,],y_unknown[13,],y_unknown[14,],y_unknown[15,],
                  y_unknown_contaminated_more[1,],y_unknown_contaminated_more[2,],
                  y_unknown_contaminated_more[3,],y_unknown_contaminated_more[4,],
                  y_unknown_contaminated_more[5,]),
    d_unknown = c(rep(c(1/10, 1/100, 1/10000), number_of_unknown_samples)),
    ind_unknown = c(rep(1:(number_of_unknown_samples),each = 3)),
    start_index = seq(1, 3*(number_of_unknown_samples - 1) + 1, by=3),
    end_index = seq(3, 3*(number_of_unknown_samples - 1) + 3, by=3),
    theta_0 = as.array(standard_initial_concentration)) 
  initfun_4para <- function(...) {
    list(beta =c(min(y_standard),max(y_standard) - min(y_standard), 1, 1))
  }
  
  simulated_model_more <- stan(file = "Bayesian_contamination_model.stan", data = simulated_data_more, init = initfun_4para)
  saveRDS(simulated_model_more, paste("simulation_more_", m, ".rds", sep=""))
  model_summmary_more <- summary(simulated_model_more, pars = c("phi"), probs = c(0.025, 0.5, 0.975))$summary


  ####base method  --- need to change log scale back to original scale
  ####base model is also named as 2004 model
  #Severe, base
  simulated_data_2004_more <- list(
    N_stand_sample = 1,
    N_unknown_sample = number_of_unknown_samples,
    N_standard = 1*26,
    N_unknown = number_of_unknown_samples*3,
    y_standard = exp(y_standard),
    d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                   0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048),
    ind_standard = c(rep(1,26)),
    y_unknown = exp(c(y_unknown[1,],y_unknown[2,],y_unknown[3,],y_unknown[4,],y_unknown[5,],
                      y_unknown[6,],y_unknown[7,],y_unknown[8,],y_unknown[9,],y_unknown[10,],
                      y_unknown[11,],y_unknown[12,],y_unknown[13,],y_unknown[14,],y_unknown[15,],
                      y_unknown_contaminated_more[1,],y_unknown_contaminated_more[2,],
                      y_unknown_contaminated_more[3,],y_unknown_contaminated_more[4,],
                      y_unknown_contaminated_more[5,])),
    d_unknown = c(rep(c(1/10, 1/100, 1/10000), number_of_unknown_samples)),
    ind_unknown = c(rep(1:(number_of_unknown_samples),each = 3)),
    start_index = seq(1, 3*(number_of_unknown_samples - 1) + 1, by=3),
    end_index = seq(3, 3*(number_of_unknown_samples - 1) + 3, by=3),
    A = mean(exp(y_standard)),
    theta_0 = as.array(standard_initial_concentration)) 
  initfun_4para <- function(...) {
    list(beta =c(min(y_standard),max(y_standard) - min(y_standard), 1, 1))
  }
 
  simulated_model_2004_more <- stan(file = "initial_model.stan", data = simulated_data_2004_more, init = initfun_4para)
  saveRDS(simulated_model_2004_more, paste("simulation_2004_more_", m, ".rds", sep=""))
  model_summmary_2004_more <- summary(simulated_model_2004_more, pars = c("phi"), probs = c(0.025, 0.5, 0.975))$summary
  
  combined_concentration <- c(unknown_concentration, unknown_concentration_contaminated)
  for (k in 1:number_of_unknown_samples) {
    estimate_our_method_more_contamination[m, k] <-  model_summmary_more[k,5]
    estimate_2004_method_more_contamination[m, k] <-  model_summmary_2004_more[k,5]
  }
}

write.csv(estimate_our_method_more_contamination, "estimate_our_method_more_contamination.csv")
write.csv(estimate_2004_method_more_contamination, "estimate_2004_method_more_contamination.csv")


