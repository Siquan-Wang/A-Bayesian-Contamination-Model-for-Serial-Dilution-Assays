####Load libraries

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
  unknown_concentration <- c(2.142117, 2.173073, 3.067431, 5.303637, 5.715181, 5.759132, 6.891102, 8.137805, 10.642271, 12.283897, 
                             19.371679, 20.043561, 23.413117, 28.621908, 32.093232)
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
  unknown_concentration_contaminated <- c(3.373859, 11.704103, 12.557783, 14.030920, 14.792289)
  number_of_contaminated_samples <- 5
  delta_sigma_contaminated_less <- 0.1
  delta_beta_contaminated_less <- 0.1
  delta_sigma_contaminated_more <- 1
  delta_beta_contaminated_more <- 1
  y_unknown_contaminated_less <- matrix(NA, nrow = length(unknown_concentration_contaminated), ncol = length(d_unknown))
  y_unknown_contaminated_more <- matrix(NA, nrow = length(unknown_concentration_contaminated), ncol = length(d_unknown))
  for (i in 1:number_of_contaminated_samples){
    for(j in 1:length(d_unknown)){
      y_unknown_contaminated_less[i,j] <- rnorm(1, mean = log(h_function_4para(unknown_concentration_contaminated[i]*d_unknown[j],
                                                                               beta_1, 
                                                                               beta_2*exp(delta_beta_contaminated_less),
                                                                               beta_3,
                                                                               beta_4)), sd = sigma_y*exp(delta_sigma_contaminated_less))
      y_unknown_contaminated_more[i,j] <- rnorm(1, mean = log(h_function_4para(unknown_concentration_contaminated[i]*d_unknown[j],
                                                                               beta_1, 
                                                                               beta_2*exp(delta_beta_contaminated_more),
                                                                               beta_3,
                                                                               beta_4)), sd = sigma_y*exp(delta_sigma_contaminated_more))
    }
  }
  ######
  combined_data_less_comtamination <- data.frame(
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
           y_unknown_contaminated_less[1,],y_unknown_contaminated_less[2,],
           y_unknown_contaminated_less[3,],y_unknown_contaminated_less[4,],
           y_unknown_contaminated_less[5,])
  )
  colnames(combined_data_less_comtamination) <- c("id", "dilution", "y")
  write.csv(combined_data_less_comtamination, file = paste("combined_data_less_comtamination", m, ".csv", sep=""))
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
  
  ###start simulation, less contamination, our proposed model
  simulated_data_less <- list(
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
                  y_unknown_contaminated_less[1,],y_unknown_contaminated_less[2,],
                  y_unknown_contaminated_less[3,],y_unknown_contaminated_less[4,],
                  y_unknown_contaminated_less[5,]),
    d_unknown = c(rep(c(1/10, 1/100, 1/10000), number_of_unknown_samples)),
    ind_unknown = c(rep(1:(number_of_unknown_samples),each = 3)),
    start_index = seq(1, 3*(number_of_unknown_samples - 1) + 1, by=3),
    end_index = seq(3, 3*(number_of_unknown_samples - 1) + 3, by=3),
    theta_0 = as.array(standard_initial_concentration)) 
  initfun_4para <- function(...) {
    list(beta =c(min(y_standard),max(y_standard) - min(y_standard), 1, 1))
  }
  
  simulated_model_less <- stan(file = "Bayesian_contamination_model.stan", data = simulated_data_less, init = initfun_4para)
  saveRDS(simulated_model_less, paste("simulation_2_less_exponential_generation", m, ".rds", sep=""))
  model_summmary_less <- summary(simulated_model_less, pars = c("phi"), probs = c(0.025, 0.5, 0.975))$summary

# severe contamination, our proposed model
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
  saveRDS(simulated_model_more, paste("simulation_2_more_exponential_generation", m, ".rds", sep=""))
  model_summmary_more <- summary(simulated_model_more, pars = c("phi"), probs = c(0.025, 0.5, 0.975))$summary
  
  ####base method  --- need to change log scale back to original scale
  # less contamination, base model
  simulated_data_2004_less <- list(
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
                      y_unknown_contaminated_less[1,],y_unknown_contaminated_less[2,],
                      y_unknown_contaminated_less[3,],y_unknown_contaminated_less[4,],
                      y_unknown_contaminated_less[5,])),
    d_unknown = c(rep(c(1/10, 1/100, 1/10000), number_of_unknown_samples)),
    ind_unknown = c(rep(1:(number_of_unknown_samples),each = 3)),
    start_index = seq(1, 3*(number_of_unknown_samples - 1) + 1, by=3),
    end_index = seq(3, 3*(number_of_unknown_samples - 1) + 3, by=3),
    A = mean(exp(y_standard)),
    theta_0 = as.array(standard_initial_concentration)) 
  initfun_4para <- function(...) {
    list(beta =c(min(y_standard),max(y_standard) - min(y_standard), 1, 1))
  }
  
  simulated_model_2004_less <- stan(file = "initial_model.stan", data = simulated_data_2004_less, init = initfun_4para)
  saveRDS(simulated_model_2004_less, paste("simulation_2_2004_less_exponential_generation", m, ".rds", sep=""))
  
  #severe contamination, base model
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
  saveRDS(simulated_model_2004_more, paste("simulation_2_2004_more_exponential_generation", m, ".rds", sep=""))
}
