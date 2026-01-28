############### Load library and data
library(tidyverse)
library(rstan)
library(ggpubr)
library(shinystan)
library(psych)
library(loo)

h_function_4para <- function(x, beta_1, beta_2, beta_3, beta_4) {
  beta_1 + beta_2/(1 + (x/beta_3)^(-beta_4))
}
derf <- read_csv("derf.csv")
initial_concentration <- 125
y_standard = c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13])

############### Fit our Bayesian computation model
derf_standard_unknown_data <- list(
  N_stand_sample = 1,
  N_unknown_sample = 23,
  N_standard = 1*26,
  N_unknown = 23*3,
  y_standard = log(c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13])),
  d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                 0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048),
  ind_standard = c(rep(1,26)),
  y_unknown = log(c(derf$FI[14:82])),
  d_unknown = c(derf$Dilution_1[14:82]),
  ind_unknown = c(rep(1:23,each = 3)),
  start_index = seq(1, 3*22+1, by=3),
  end_index = seq(3, 3*22+3, by=3),
  theta_0 = as.array(initial_concentration)) 
initfun_4para <- function(...) {
  list(beta =c(min(y_standard),max(y_standard) - min(y_standard),1, 1))
}

##### Standard + Unknown, model 1
model_fit <- stan(file = "Bayesian_contamination_model.stan", data = derf_standard_unknown_data, init = initfun_4para)
saveRDS(model_fit, "model_fit.rds")

