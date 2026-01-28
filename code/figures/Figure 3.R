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

derf <- read_csv("derf.csv")
y_standard = log(c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13]))
d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
               0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048)
y_unknown = log(c(derf$FI[14:82]))
d_unknown = c(derf$Dilution_1[14:82])
standard_concentration = 125

d_zm = as.vector(seq(0,0.1,0.0001))
d_zm_standard = as.vector(seq(0,1,0.0001)) 

standard_unknown <- readRDS("model_fit.rds")  
matrix_of_draws <- as.data.frame(standard_unknown, par = c("lambda", "beta", "beta2_unknown", "phi"))
matrix_of_draws_logprob_latent1 <- as.data.frame(standard_unknown, par = c("lp_mixture_latent_1"))
matrix_of_draws_logprob_latent2 <- as.data.frame(standard_unknown, par = c("lp_mixture_latent_2"))
set.seed(1234)
number_of_draws_more <- 4000
random_sample_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
start_index = seq(1, 3*22+1, by=3)
end_index = seq(3, 3*22+3, by=3)
colnames(matrix_of_draws) <- c("lambda","beta_1","beta_2","beta_3","beta_4","beta2_unknown_1","beta2_unknown_2","beta2_unknown_3",
                               "beta2_unknown_4","beta2_unknown_5","beta2_unknown_6","beta2_unknown_7","beta2_unknown_8","beta2_unknown_9","beta2_unknown_10","beta2_unknown_11",
                               "beta2_unknown_12","beta2_unknown_13","beta2_unknown_14","beta2_unknown_15","beta2_unknown_16","beta2_unknown_17","beta2_unknown_18","beta2_unknown_19",
                               "beta2_unknown_20","beta2_unknown_21","beta2_unknown_22","beta2_unknown_23","phi_1","phi_2","phi_3","phi_4",
                               "phi_5","phi_6","phi_7","phi_8","phi_9","phi_10","phi_11","phi_12",
                               "phi_13","phi_14","phi_15","phi_16","phi_17","phi_18","phi_19","phi_20",
                               "phi_21","phi_22","phi_23")
colnames(matrix_of_draws_logprob_latent1) <- c("lp_mixture_latent_1_1","lp_mixture_latent_1_2","lp_mixture_latent_1_3","lp_mixture_latent_1_4",
                               "lp_mixture_latent_1_5","lp_mixture_latent_1_6","lp_mixture_latent_1_7","lp_mixture_latent_1_8","lp_mixture_latent_1_9","lp_mixture_latent_1_10","lp_mixture_latent_1_11","lp_mixture_latent_1_12",
                               "lp_mixture_latent_1_13","lp_mixture_latent_1_14","lp_mixture_latent_1_15","lp_mixture_latent_1_16","lp_mixture_latent_1_17","lp_mixture_latent_1_18","lp_mixture_latent_1_19","lp_mixture_latent_1_20",
                               "lp_mixture_latent_1_21","lp_mixture_latent_1_22","lp_mixture_latent_1_23")
colnames(matrix_of_draws_logprob_latent2) <- c("lp_mixture_latent_2_1","lp_mixture_latent_2_2","lp_mixture_latent_2_3","lp_mixture_latent_2_4",
                                       "lp_mixture_latent_2_5","lp_mixture_latent_2_6","lp_mixture_latent_2_7","lp_mixture_latent_2_8","lp_mixture_latent_2_9","lp_mixture_latent_2_10","lp_mixture_latent_2_11","lp_mixture_latent_2_12",
                                       "lp_mixture_latent_2_13","lp_mixture_latent_2_14","lp_mixture_latent_2_15","lp_mixture_latent_2_16","lp_mixture_latent_2_17","lp_mixture_latent_2_18","lp_mixture_latent_2_19","lp_mixture_latent_2_20",
                                       "lp_mixture_latent_2_21","lp_mixture_latent_2_22","lp_mixture_latent_2_23" )
### Generate results based on posterior draws
random_sample_lambda_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
mixture_probability_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
mixture_probability_latent1_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
mixture_probability_latent2_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
for (i in 1:23){
  random_sample_more[i,] <- sample(seq(1:4000), number_of_draws_more)
  for (j in 1:number_of_draws_more){
    mixture_probability_latent1_more[i,j] <- prod(as.vector(exp(matrix_of_draws_logprob_latent1[,start_index[i]:end_index[i]])[random_sample_more[i, j],]))
    mixture_probability_latent2_more[i,j] <- prod(as.vector(exp(matrix_of_draws_logprob_latent2[,start_index[i]:end_index[i]])[random_sample_more[i, j],]))
    mixture_probability_more[i,j] <- matrix_of_draws$lambda[random_sample_more[i, j]]*mixture_probability_latent2_more[i,j]/(matrix_of_draws$lambda[random_sample_more[i, j]]*mixture_probability_latent2_more[i,j] + (1-matrix_of_draws$lambda[random_sample_more[i, j]])*mixture_probability_latent1_more[i,j])
    random_sample_lambda_more[i,j] <- rbinom(prob = mixture_probability_more[i,j], size = 1, n = 1)
  }
}
d_y1_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y2_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y3_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y4_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y5_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y6_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y7_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y8_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y9_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y10_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y11_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y12_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y13_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y14_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y15_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y16_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y17_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y18_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y19_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y20_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y21_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y22_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y23_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y_standard_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm_standard))

for (i in 1:number_of_draws_more){
  d_y1_more[i,] <- (1 -random_sample_lambda_more[1,i])*log(h_function_4para(matrix_of_draws$phi_1[random_sample_more[1,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[1,i]],matrix_of_draws$beta_2[random_sample_more[1,i]],matrix_of_draws$beta_3[random_sample_more[1,i]],matrix_of_draws$beta_4[random_sample_more[1,i]])) + 
    random_sample_lambda_more[1,i]*log(h_function_4para(matrix_of_draws$phi_1[random_sample_more[1,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[1,i]],matrix_of_draws$beta2_unknown_1[random_sample_more[1,i]],matrix_of_draws$beta_3[random_sample_more[1,i]],matrix_of_draws$beta_4[random_sample_more[1,i]]))
  d_y2_more[i,] <- (1 -random_sample_lambda_more[2,i])*log(h_function_4para(matrix_of_draws$phi_2[random_sample_more[2,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[2,i]],matrix_of_draws$beta_2[random_sample_more[2,i]],matrix_of_draws$beta_3[random_sample_more[2,i]],matrix_of_draws$beta_4[random_sample_more[2,i]])) + 
    random_sample_lambda_more[2,i]*log(h_function_4para(matrix_of_draws$phi_2[random_sample_more[2,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[2,i]],matrix_of_draws$beta2_unknown_2[random_sample_more[2,i]],matrix_of_draws$beta_3[random_sample_more[2,i]],matrix_of_draws$beta_4[random_sample_more[2,i]]))
  d_y3_more[i,] <- (1 -random_sample_lambda_more[3,i])*log(h_function_4para(matrix_of_draws$phi_3[random_sample_more[3,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[3,i]],matrix_of_draws$beta_2[random_sample_more[3,i]],matrix_of_draws$beta_3[random_sample_more[3,i]],matrix_of_draws$beta_4[random_sample_more[3,i]])) + 
    random_sample_lambda_more[3,i]*log(h_function_4para(matrix_of_draws$phi_3[random_sample_more[3,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[3,i]],matrix_of_draws$beta2_unknown_3[random_sample_more[3,i]],matrix_of_draws$beta_3[random_sample_more[3,i]],matrix_of_draws$beta_4[random_sample_more[3,i]]))
  d_y4_more[i,] <- (1 -random_sample_lambda_more[4,i])*log(h_function_4para(matrix_of_draws$phi_4[random_sample_more[4,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[4,i]],matrix_of_draws$beta_2[random_sample_more[4,i]],matrix_of_draws$beta_3[random_sample_more[4,i]],matrix_of_draws$beta_4[random_sample_more[4,i]])) + 
    random_sample_lambda_more[4,i]*log(h_function_4para(matrix_of_draws$phi_4[random_sample_more[4,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[4,i]],matrix_of_draws$beta2_unknown_4[random_sample_more[4,i]],matrix_of_draws$beta_3[random_sample_more[4,i]],matrix_of_draws$beta_4[random_sample_more[4,i]]))
  d_y5_more[i,] <- (1 -random_sample_lambda_more[5,i])*log(h_function_4para(matrix_of_draws$phi_5[random_sample_more[5,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[5,i]],matrix_of_draws$beta_2[random_sample_more[5,i]],matrix_of_draws$beta_3[random_sample_more[5,i]],matrix_of_draws$beta_4[random_sample_more[5,i]])) + 
    random_sample_lambda_more[5,i]*log(h_function_4para(matrix_of_draws$phi_5[random_sample_more[5,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[5,i]],matrix_of_draws$beta2_unknown_5[random_sample_more[5,i]],matrix_of_draws$beta_3[random_sample_more[5,i]],matrix_of_draws$beta_4[random_sample_more[5,i]]))
  d_y6_more[i,] <- (1 -random_sample_lambda_more[6,i])*log(h_function_4para(matrix_of_draws$phi_6[random_sample_more[6,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[6,i]],matrix_of_draws$beta_2[random_sample_more[6,i]],matrix_of_draws$beta_3[random_sample_more[6,i]],matrix_of_draws$beta_4[random_sample_more[6,i]])) + 
    random_sample_lambda_more[6,i]*log(h_function_4para(matrix_of_draws$phi_6[random_sample_more[6,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[6,i]],matrix_of_draws$beta2_unknown_6[random_sample_more[6,i]],matrix_of_draws$beta_3[random_sample_more[6,i]],matrix_of_draws$beta_4[random_sample_more[6,i]]))
  d_y7_more[i,] <- (1 -random_sample_lambda_more[7,i])*log(h_function_4para(matrix_of_draws$phi_7[random_sample_more[7,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[7,i]],matrix_of_draws$beta_2[random_sample_more[7,i]],matrix_of_draws$beta_3[random_sample_more[7,i]],matrix_of_draws$beta_4[random_sample_more[7,i]])) + 
    random_sample_lambda_more[7,i]*log(h_function_4para(matrix_of_draws$phi_7[random_sample_more[7,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[7,i]],matrix_of_draws$beta2_unknown_7[random_sample_more[7,i]],matrix_of_draws$beta_3[random_sample_more[7,i]],matrix_of_draws$beta_4[random_sample_more[7,i]]))
  d_y8_more[i,] <- (1 -random_sample_lambda_more[8,i])*log(h_function_4para(matrix_of_draws$phi_8[random_sample_more[8,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[8,i]],matrix_of_draws$beta_2[random_sample_more[8,i]],matrix_of_draws$beta_3[random_sample_more[8,i]],matrix_of_draws$beta_4[random_sample_more[8,i]])) + 
    random_sample_lambda_more[8,i]*log(h_function_4para(matrix_of_draws$phi_8[random_sample_more[8,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[8,i]],matrix_of_draws$beta2_unknown_8[random_sample_more[8,i]],matrix_of_draws$beta_3[random_sample_more[8,i]],matrix_of_draws$beta_4[random_sample_more[8,i]]))
  d_y9_more[i,] <- (1 -random_sample_lambda_more[9,i])*log(h_function_4para(matrix_of_draws$phi_9[random_sample_more[9,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[9,i]],matrix_of_draws$beta_2[random_sample_more[9,i]],matrix_of_draws$beta_3[random_sample_more[9,i]],matrix_of_draws$beta_4[random_sample_more[9,i]])) + 
    random_sample_lambda_more[9,i]*log(h_function_4para(matrix_of_draws$phi_9[random_sample_more[9,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[9,i]],matrix_of_draws$beta2_unknown_9[random_sample_more[9,i]],matrix_of_draws$beta_3[random_sample_more[9,i]],matrix_of_draws$beta_4[random_sample_more[9,i]]))
  d_y10_more[i,] <- (1 -random_sample_lambda_more[10,i])*log(h_function_4para(matrix_of_draws$phi_10[random_sample_more[10,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[10,i]],matrix_of_draws$beta_2[random_sample_more[10,i]],matrix_of_draws$beta_3[random_sample_more[10,i]],matrix_of_draws$beta_4[random_sample_more[10,i]])) + 
    random_sample_lambda_more[10,i]*log(h_function_4para(matrix_of_draws$phi_10[random_sample_more[10,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[10,i]],matrix_of_draws$beta2_unknown_10[random_sample_more[10,i]],matrix_of_draws$beta_3[random_sample_more[10,i]],matrix_of_draws$beta_4[random_sample_more[10,i]]))
  d_y11_more[i,] <- (1 -random_sample_lambda_more[11,i])*log(h_function_4para(matrix_of_draws$phi_11[random_sample_more[11,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[11,i]],matrix_of_draws$beta_2[random_sample_more[11,i]],matrix_of_draws$beta_3[random_sample_more[11,i]],matrix_of_draws$beta_4[random_sample_more[11,i]])) + 
    random_sample_lambda_more[11,i]*log(h_function_4para(matrix_of_draws$phi_11[random_sample_more[11,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[11,i]],matrix_of_draws$beta2_unknown_11[random_sample_more[11,i]],matrix_of_draws$beta_3[random_sample_more[11,i]],matrix_of_draws$beta_4[random_sample_more[11,i]]))
  d_y12_more[i,] <- (1 -random_sample_lambda_more[12,i])*log(h_function_4para(matrix_of_draws$phi_12[random_sample_more[12,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[12,i]],matrix_of_draws$beta_2[random_sample_more[12,i]],matrix_of_draws$beta_3[random_sample_more[12,i]],matrix_of_draws$beta_4[random_sample_more[12,i]])) + 
    random_sample_lambda_more[12,i]*log(h_function_4para(matrix_of_draws$phi_12[random_sample_more[12,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[12,i]],matrix_of_draws$beta2_unknown_12[random_sample_more[12,i]],matrix_of_draws$beta_3[random_sample_more[12,i]],matrix_of_draws$beta_4[random_sample_more[12,i]]))
  d_y13_more[i,] <- (1 -random_sample_lambda_more[13,i])*log(h_function_4para(matrix_of_draws$phi_13[random_sample_more[13,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[13,i]],matrix_of_draws$beta_2[random_sample_more[13,i]],matrix_of_draws$beta_3[random_sample_more[13,i]],matrix_of_draws$beta_4[random_sample_more[13,i]])) + 
    random_sample_lambda_more[13,i]*log(h_function_4para(matrix_of_draws$phi_13[random_sample_more[13,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[13,i]],matrix_of_draws$beta2_unknown_13[random_sample_more[13,i]],matrix_of_draws$beta_3[random_sample_more[13,i]],matrix_of_draws$beta_4[random_sample_more[13,i]]))
  d_y14_more[i,] <- (1 -random_sample_lambda_more[14,i])*log(h_function_4para(matrix_of_draws$phi_14[random_sample_more[14,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[14,i]],matrix_of_draws$beta_2[random_sample_more[14,i]],matrix_of_draws$beta_3[random_sample_more[14,i]],matrix_of_draws$beta_4[random_sample_more[14,i]])) + 
    random_sample_lambda_more[14,i]*log(h_function_4para(matrix_of_draws$phi_14[random_sample_more[14,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[14,i]],matrix_of_draws$beta2_unknown_14[random_sample_more[14,i]],matrix_of_draws$beta_3[random_sample_more[14,i]],matrix_of_draws$beta_4[random_sample_more[14,i]]))
  d_y15_more[i,] <- (1 -random_sample_lambda_more[15,i])*log(h_function_4para(matrix_of_draws$phi_15[random_sample_more[15,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[15,i]],matrix_of_draws$beta_2[random_sample_more[15,i]],matrix_of_draws$beta_3[random_sample_more[15,i]],matrix_of_draws$beta_4[random_sample_more[15,i]])) + 
    random_sample_lambda_more[15,i]*log(h_function_4para(matrix_of_draws$phi_15[random_sample_more[15,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[15,i]],matrix_of_draws$beta2_unknown_15[random_sample_more[15,i]],matrix_of_draws$beta_3[random_sample_more[15,i]],matrix_of_draws$beta_4[random_sample_more[15,i]]))
  d_y16_more[i,] <- (1 -random_sample_lambda_more[16,i])*log(h_function_4para(matrix_of_draws$phi_16[random_sample_more[16,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[16,i]],matrix_of_draws$beta_2[random_sample_more[16,i]],matrix_of_draws$beta_3[random_sample_more[16,i]],matrix_of_draws$beta_4[random_sample_more[16,i]])) + 
    random_sample_lambda_more[16,i]*log(h_function_4para(matrix_of_draws$phi_16[random_sample_more[16,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[16,i]],matrix_of_draws$beta2_unknown_16[random_sample_more[16,i]],matrix_of_draws$beta_3[random_sample_more[16,i]],matrix_of_draws$beta_4[random_sample_more[16,i]]))
  d_y17_more[i,] <- (1 -random_sample_lambda_more[17,i])*log(h_function_4para(matrix_of_draws$phi_17[random_sample_more[17,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[17,i]],matrix_of_draws$beta_2[random_sample_more[17,i]],matrix_of_draws$beta_3[random_sample_more[17,i]],matrix_of_draws$beta_4[random_sample_more[17,i]])) + 
    random_sample_lambda_more[17,i]*log(h_function_4para(matrix_of_draws$phi_17[random_sample_more[17,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[17,i]],matrix_of_draws$beta2_unknown_17[random_sample_more[17,i]],matrix_of_draws$beta_3[random_sample_more[17,i]],matrix_of_draws$beta_4[random_sample_more[17,i]]))
  d_y18_more[i,] <- (1 -random_sample_lambda_more[18,i])*log(h_function_4para(matrix_of_draws$phi_18[random_sample_more[18,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[18,i]],matrix_of_draws$beta_2[random_sample_more[18,i]],matrix_of_draws$beta_3[random_sample_more[18,i]],matrix_of_draws$beta_4[random_sample_more[18,i]])) + 
    random_sample_lambda_more[18,i]*log(h_function_4para(matrix_of_draws$phi_18[random_sample_more[18,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[18,i]],matrix_of_draws$beta2_unknown_18[random_sample_more[18,i]],matrix_of_draws$beta_3[random_sample_more[18,i]],matrix_of_draws$beta_4[random_sample_more[18,i]]))
  d_y19_more[i,] <- (1 -random_sample_lambda_more[19,i])*log(h_function_4para(matrix_of_draws$phi_19[random_sample_more[19,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[19,i]],matrix_of_draws$beta_2[random_sample_more[19,i]],matrix_of_draws$beta_3[random_sample_more[19,i]],matrix_of_draws$beta_4[random_sample_more[19,i]])) + 
    random_sample_lambda_more[19,i]*log(h_function_4para(matrix_of_draws$phi_19[random_sample_more[19,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[19,i]],matrix_of_draws$beta2_unknown_19[random_sample_more[19,i]],matrix_of_draws$beta_3[random_sample_more[19,i]],matrix_of_draws$beta_4[random_sample_more[19,i]]))
  d_y20_more[i,] <- (1 -random_sample_lambda_more[20,i])*log(h_function_4para(matrix_of_draws$phi_20[random_sample_more[20,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[20,i]],matrix_of_draws$beta_2[random_sample_more[20,i]],matrix_of_draws$beta_3[random_sample_more[20,i]],matrix_of_draws$beta_4[random_sample_more[20,i]])) + 
    random_sample_lambda_more[20,i]*log(h_function_4para(matrix_of_draws$phi_20[random_sample_more[20,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[20,i]],matrix_of_draws$beta2_unknown_20[random_sample_more[20,i]],matrix_of_draws$beta_3[random_sample_more[20,i]],matrix_of_draws$beta_4[random_sample_more[20,i]]))
  d_y21_more[i,] <- (1 -random_sample_lambda_more[21,i])*log(h_function_4para(matrix_of_draws$phi_21[random_sample_more[21,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[21,i]],matrix_of_draws$beta_2[random_sample_more[21,i]],matrix_of_draws$beta_3[random_sample_more[21,i]],matrix_of_draws$beta_4[random_sample_more[21,i]])) + 
    random_sample_lambda_more[21,i]*log(h_function_4para(matrix_of_draws$phi_21[random_sample_more[21,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[21,i]],matrix_of_draws$beta2_unknown_21[random_sample_more[21,i]],matrix_of_draws$beta_3[random_sample_more[21,i]],matrix_of_draws$beta_4[random_sample_more[21,i]]))
  d_y22_more[i,] <- (1 -random_sample_lambda_more[22,i])*log(h_function_4para(matrix_of_draws$phi_22[random_sample_more[22,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[22,i]],matrix_of_draws$beta_2[random_sample_more[22,i]],matrix_of_draws$beta_3[random_sample_more[22,i]],matrix_of_draws$beta_4[random_sample_more[22,i]])) + 
    random_sample_lambda_more[22,i]*log(h_function_4para(matrix_of_draws$phi_22[random_sample_more[22,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[22,i]],matrix_of_draws$beta2_unknown_22[random_sample_more[22,i]],matrix_of_draws$beta_3[random_sample_more[22,i]],matrix_of_draws$beta_4[random_sample_more[22,i]]))
  d_y23_more[i,] <- (1 -random_sample_lambda_more[23,i])*log(h_function_4para(matrix_of_draws$phi_23[random_sample_more[23,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[23,i]],matrix_of_draws$beta_2[random_sample_more[23,i]],matrix_of_draws$beta_3[random_sample_more[23,i]],matrix_of_draws$beta_4[random_sample_more[23,i]])) + 
    random_sample_lambda_more[23,i]*log(h_function_4para(matrix_of_draws$phi_23[random_sample_more[23,i]]*d_zm,matrix_of_draws$beta_1[random_sample_more[23,i]],matrix_of_draws$beta2_unknown_23[random_sample_more[23,i]],matrix_of_draws$beta_3[random_sample_more[23,i]],matrix_of_draws$beta_4[random_sample_more[23,i]]))
  d_y_standard_more[i,] <- log(h_function_4para(125*d_zm_standard,matrix_of_draws$beta_1[random_sample_more[23,i]],matrix_of_draws$beta_2[random_sample_more[23,i]],matrix_of_draws$beta_3[random_sample_more[23,i]],matrix_of_draws$beta_4[random_sample_more[23,i]]))
}
d_y_standard_median <- apply(d_y_standard_more,2, median)
d_y_standard_up <- apply(d_y_standard_more,2, quantile, prob = 0.975)
d_y_standard_down <- apply(d_y_standard_more,2, quantile, prob = 0.025)
d_y1_median <- apply(d_y1_more,2, median)
d_y1_up <- apply(d_y1_more,2, quantile, prob = 0.975)
d_y1_down <- apply(d_y1_more,2, quantile, prob = 0.025)
d_y2_median <- apply(d_y2_more,2, median)
d_y2_up <- apply(d_y2_more,2, quantile, prob = 0.975)
d_y2_down <- apply(d_y2_more,2, quantile, prob = 0.025)
d_y3_median <- apply(d_y3_more,2, median)
d_y3_up <- apply(d_y3_more,2, quantile, prob = 0.975)
d_y3_down <- apply(d_y3_more,2, quantile, prob = 0.025)
d_y4_median <- apply(d_y4_more,2, median)
d_y4_up <- apply(d_y4_more,2, quantile, prob = 0.975)
d_y4_down <- apply(d_y4_more,2, quantile, prob = 0.025)
d_y5_median <- apply(d_y5_more,2, median)
d_y5_up <- apply(d_y5_more,2, quantile, prob = 0.975)
d_y5_down <- apply(d_y5_more,2, quantile, prob = 0.025)
d_y6_median <- apply(d_y6_more,2, median)
d_y6_up <- apply(d_y6_more,2, quantile, prob = 0.975)
d_y6_down <- apply(d_y6_more,2, quantile, prob = 0.025)
d_y7_median <- apply(d_y7_more,2, median)
d_y7_up <- apply(d_y7_more,2, quantile, prob = 0.975)
d_y7_down <- apply(d_y7_more,2, quantile, prob = 0.025)
d_y8_median <- apply(d_y8_more,2, median)
d_y8_up <- apply(d_y8_more,2, quantile, prob = 0.975)
d_y8_down <- apply(d_y8_more,2, quantile, prob = 0.025)
d_y9_median <- apply(d_y9_more,2, median)
d_y9_up <- apply(d_y9_more,2, quantile, prob = 0.975)
d_y9_down <- apply(d_y9_more,2, quantile, prob = 0.025)
d_y10_median <- apply(d_y10_more,2, median)
d_y10_up <- apply(d_y10_more,2, quantile, prob = 0.975)
d_y10_down <- apply(d_y10_more,2, quantile, prob = 0.025)
d_y11_median <- apply(d_y11_more,2, median)
d_y11_up <- apply(d_y11_more,2, quantile, prob = 0.975)
d_y11_down <- apply(d_y11_more,2, quantile, prob = 0.025)
d_y12_median <- apply(d_y12_more,2, median)
d_y12_up <- apply(d_y12_more,2, quantile, prob = 0.975)
d_y12_down <- apply(d_y12_more,2, quantile, prob = 0.025)
d_y13_median <- apply(d_y13_more,2, median)
d_y13_up <- apply(d_y13_more,2, quantile, prob = 0.975)
d_y13_down <- apply(d_y13_more,2, quantile, prob = 0.025)
d_y14_median <- apply(d_y14_more,2, median)
d_y14_up <- apply(d_y14_more,2, quantile, prob = 0.975)
d_y14_down <- apply(d_y14_more,2, quantile, prob = 0.025)
d_y15_median <- apply(d_y15_more,2, median)
d_y15_up <- apply(d_y15_more,2, quantile, prob = 0.975)
d_y15_down <- apply(d_y15_more,2, quantile, prob = 0.025)
d_y16_median <- apply(d_y16_more,2, median)
d_y16_up <- apply(d_y16_more,2, quantile, prob = 0.975)
d_y16_down <- apply(d_y16_more,2, quantile, prob = 0.025)
d_y17_median <- apply(d_y17_more,2, median)
d_y17_up <- apply(d_y17_more,2, quantile, prob = 0.975)
d_y17_down <- apply(d_y17_more,2, quantile, prob = 0.025)
d_y18_median <- apply(d_y18_more,2, median)
d_y18_up <- apply(d_y18_more,2, quantile, prob = 0.975)
d_y18_down <- apply(d_y18_more,2, quantile, prob = 0.025)
d_y19_median <- apply(d_y19_more,2, median)
d_y19_up <- apply(d_y19_more,2, quantile, prob = 0.975)
d_y19_down <- apply(d_y19_more,2, quantile, prob = 0.025)
d_y20_median <- apply(d_y20_more,2, median)
d_y20_up <- apply(d_y20_more,2, quantile, prob = 0.975)
d_y20_down <- apply(d_y20_more,2, quantile, prob = 0.025)
d_y21_median <- apply(d_y21_more,2, median)
d_y21_up <- apply(d_y21_more,2, quantile, prob = 0.975)
d_y21_down <- apply(d_y21_more,2, quantile, prob = 0.025)
d_y22_median <- apply(d_y22_more,2, median)
d_y22_up <- apply(d_y22_more,2, quantile, prob = 0.975)
d_y22_down <- apply(d_y22_more,2, quantile, prob = 0.025)
d_y23_median <- apply(d_y23_more,2, median)
d_y23_up <- apply(d_y23_more,2, quantile, prob = 0.975)
d_y23_down <- apply(d_y23_more,2, quantile, prob = 0.025)

### Plots for standard data
zoom_1 <-  ggplot()+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_median), size = 1,color = "black", linetype = "solid")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_up), size = 1,color = "gray", linetype = "dashed")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_down), size = 1,color = "gray", linetype = "dashed")+
  geom_point(aes(x = d_standard, y = y_standard), size = 1.5)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Dust Mite Plate Standard Sample")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(0,1.05), ylim = c(0, max(y_standard) *1.1))

zoom_2 <-  ggplot()+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_median), size = 1,color = "black", linetype = "solid")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_up), size = 1,color = "gray", linetype = "dashed")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_down), size = 1,color = "gray", linetype = "dashed")+
  geom_point(aes(x = d_standard, y = y_standard), size = 1.5)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Dust Mite Plate Standard Sample")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(0,0.1))

zoom_3 <-  ggplot()+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_median), size = 1,color = "black", linetype = "solid")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_up), size = 1,color = "gray", linetype = "dashed")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_down), size = 1,color = "gray", linetype = "dashed")+
  geom_point(aes(x = d_standard, y = y_standard), size = 1.5)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Dust Mite Plate Standard Sample")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(0,1.05*d_zm_standard[which(d_y_standard_up > 5)[1]]), ylim = c(0, 5))


contamination_summary <- summary(standard_unknown, pars = c("p_tilde"), probs = c(0.25, 0.5, 0.75))$summary
print(contamination_summary)
phi_summary <- summary(standard_unknown, pars = c("phi"), probs = c(0.25, 0.5, 0.75))$summary
#write.csv(phi_summary, "phi_summary.csv")
### 2004 plots for unknown data (fitting with standard + unknown)
library(scales)
show_col(hue_pal()(4))
legend_title <- "line type"
plot1 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y1_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y1_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y1_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[14:16])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.98, 50% interval [0.18,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ # change linetype
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 1")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    #axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    #axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 4))

plot2 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y2_median, linetype = "posterior median"), size = 1,color = "black")+ 
  geom_line(aes(x = d_zm, y = d_y2_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y2_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[17:19])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.75,label = "0.001, 50% interval [0,0.10]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 2")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 8))

plot3 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y3_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y3_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y3_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[20:22])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "1, 50% interval [1,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 3")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))

plot4 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y4_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y4_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y4_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[23:25])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.75,label = "0.002, 50% interval [0,0.18]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 4")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 8))

plot5 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y5_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y5_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y5_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[26:28])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.75,label = "0.001, 50% interval [0,0.09]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 5")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 8))

plot6 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y6_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y6_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y6_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[29:31])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.59, 50% interval [0.01,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 6")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    #axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 5))

plot7 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y7_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y7_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y7_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[32:34])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.005, 50% interval [0,0.25]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 7")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 6))

plot8 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y8_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y8_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y8_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[35:37])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.001, 50% interval [0,0.06]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 8")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    #axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 6))

plot9 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y9_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y9_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y9_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[38:40])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "1, 50% interval [1,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 9")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))

plot10 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y10_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y10_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y10_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[41:43])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "0.16, 50% interval [0,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 10")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    #axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))

plot11 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y11_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y11_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y11_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[44:46])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.07, 50% interval [0,0.98]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 11")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 5))

plot12 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y12_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y12_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y12_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[47:49])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.04, 50% interval [0,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 12")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 6))

plot13 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y13_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y13_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y13_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[50:52])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.75,label = "0.0005, 50% interval [0,0.05]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+ 
  ylab("logy")+
  labs(title = "Unknown sample 13")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 8))

plot14 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y14_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y14_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y14_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[53:55])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.002, 50% interval [0,0.16]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 14")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0,6))

plot15 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y15_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y15_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y15_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[56:58])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.06, 50% interval [0,0.69]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 15")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    #axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 4))

plot16 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y16_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y16_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y16_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[59:61])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "1, 50% interval [1,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 16")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    #axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0,4))

plot17 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y17_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y17_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y17_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[62:64])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "1, 50% interval [1,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 17")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))

plot18 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y18_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y18_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y18_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[65:67])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.004, 50% interval [0,0.12]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 18")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 5))

plot19 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y19_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y19_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y19_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[68:70])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.04, 50% interval [0,0.99]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 19")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 6))

plot20 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y20_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y20_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y20_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[71:73])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "0.006, 50% interval [0,0.21]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 20")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 5))

plot21 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y21_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y21_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y21_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[74:76])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.75,label = "0.003, 50% interval [0,0.23]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 21")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    #axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 8))

plot22 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y22_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y22_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y22_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[77:79])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 0.5,label = "1, 50% interval [0.74,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 22")+
  theme( legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0,5))

plot23 <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y23_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y23_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y23_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[80:82])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "1, 50% interval [1,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 23")+
  theme(legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    #axis.ticks.x=element_blank(),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y = element_text(color = "grey20", size = 10),
    axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))

plot23_legend <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y23_median, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y23_up, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y23_down, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[80:82])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "50% posterior probability of contamination: [1,1]"), size = 2.5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unknown sample 23")+
  theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major  = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = 8, angle = 25),
        #axis.ticks.x=element_blank(),
        axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 10),
        axis.title.y = element_blank() #element_text(size = 14.5)
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))
legend <- get_legend(plot23_legend)
A_plot <-ggarrange(zoom_1,zoom_2,zoom_3, legend, ncol = 4)

col_1 <- ggarrange(plot10,plot21,plot8,plot6, nrow = 5, align = "v" )#,common.legend = TRUE)
col_2 <- ggarrange(plot17,plot5,plot7,plot20,plot1, nrow = 5, align = "v" )#,common.legend = TRUE)
col_3 <- ggarrange(plot23,plot4,plot12,plot22,plot16, nrow = 5, align = "v" )#,common.legend = TRUE)
col_4 <- ggarrange(plot9,plot2,plot14,plot18,plot15, nrow = 5, align = "v" )#,common.legend = TRUE)
col_5 <- ggarrange(plot3,plot13,plot19,plot11, nrow = 5, align = "v" )#,common.legend = TRUE)

ggdraw() +
  draw_plot(A_plot, x = 0.1, y = 1/6*5, width = 0.8, height = 1/6) +
  draw_plot(col_1, x = 0, y = 0, width = 0.2, height = 5/6) +
  draw_plot(col_2, x = 1/5, y = 0, width = 0.2, height = 5/6)+
  draw_plot(col_3, x = 2/5, y = 0, width = 0.2, height = 5/6)+
  draw_plot(col_4, x = 3/5, y = 0, width = 0.2, height = 5/6)+
  draw_plot(col_5, x = 4/5, y = 0, width = 0.2, height = 5/6)
ggsave("Figure3.png", width = 18, height = 12)

