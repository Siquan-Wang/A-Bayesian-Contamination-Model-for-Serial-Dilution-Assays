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
initial_concentration = 125
y_standard = c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13])
d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
               0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048)
y_unknown = c(derf$FI[14:82])
d_unknown = c(derf$Dilution_1[14:82])
### Generate the results for the initial model
data_initial <- list(
  N_stand_sample = 1,
  N_unknown_sample = 23,
  N_standard = 1*26,
  N_unknown = 23*3,
  y_standard = c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13]),
  d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
                 0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048),
  ind_standard = c(rep(1,26)),
  y_unknown = c(derf$FI[14:82]),
  d_unknown = c(derf$Dilution_1[14:82]),
  ind_unknown = c(rep(1:23,each = 3)),
  start_index = seq(1, 3*22+1, by=3),
  end_index = seq(3, 3*22+3, by=3),
  A = mean(y_standard),
  theta_0 = as.array(initial_concentration))
initfun_4para <- function(...) {
  list(beta =c(min(y_standard),max(y_standard) - min(y_standard), 1, 1))
}

initial_model <- stan(file = "initial_model.stan", data = data_initial, init = initfun_4para)

### Generate the results for the intermediate model
derf_intermediate <- list(
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
  theta_0 = as.array(initial_concentration)) ##Be careful that you need to change this number
initfun_4para <- function(...) {
  list(beta =c(min(y_standard),max(y_standard) - min(y_standard),1, 1))
}

intermediate_model <- stan(file = "intermediate_model.stan", data = derf_intermediate, init = initfun_4para)

### Generating Plots
final_model <- readRDS("model_fit.rds")
matrix_of_draws_initial <- as.data.frame(initial_model, par = c("beta", "phi"))
matrix_of_draws_intermediate <- as.data.frame(intermediate_model, par = c("beta", "phi"))
matrix_of_draws_final <- as.data.frame(standard_unknown, par = c("lambda", "beta", "beta2_unknown", "phi"))
set.seed(1234)
number_of_draws_more <- 4000
random_sample_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
start_index = seq(1, 3*22+1, by=3)
end_index = seq(3, 3*22+3, by=3)
d_zm_standard = as.vector(seq(0,1,0.0001)) 
d_zm = as.vector(seq(0,0.1,0.0001))
colnames(matrix_of_draws_initial) <- c("beta_1","beta_2","beta_3","beta_4","phi_1","phi_2","phi_3","phi_4",
                               "phi_5","phi_6","phi_7","phi_8","phi_9","phi_10","phi_11","phi_12",
                               "phi_13","phi_14","phi_15","phi_16","phi_17","phi_18","phi_19","phi_20",
                               "phi_21","phi_22","phi_23")
colnames(matrix_of_draws_intermediate) <- c("beta_1","beta_2","beta_3","beta_4","phi_1","phi_2","phi_3","phi_4",
                                       "phi_5","phi_6","phi_7","phi_8","phi_9","phi_10","phi_11","phi_12",
                                       "phi_13","phi_14","phi_15","phi_16","phi_17","phi_18","phi_19","phi_20",
                                       "phi_21","phi_22","phi_23")
colnames(matrix_of_draws_final) <- c("lambda","beta_1","beta_2","beta_3","beta_4","beta2_unknown_1","beta2_unknown_2","beta2_unknown_3",
                                     "beta2_unknown_4","beta2_unknown_5","beta2_unknown_6","beta2_unknown_7","beta2_unknown_8","beta2_unknown_9","beta2_unknown_10","beta2_unknown_11",
                                     "beta2_unknown_12","beta2_unknown_13","beta2_unknown_14","beta2_unknown_15","beta2_unknown_16","beta2_unknown_17","beta2_unknown_18","beta2_unknown_19",
                                     "beta2_unknown_20","beta2_unknown_21","beta2_unknown_22","beta2_unknown_23","phi_1","phi_2","phi_3","phi_4",
                                     "phi_5","phi_6","phi_7","phi_8","phi_9","phi_10","phi_11","phi_12",
                                     "phi_13","phi_14","phi_15","phi_16","phi_17","phi_18","phi_19","phi_20",
                                     "phi_21","phi_22","phi_23")
###Generate results from posterior distribution draws
random_sample_lambda_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
mixture_probability_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
mixture_probability_latent1_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
mixture_probability_latent2_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
matrix_of_draws_logprob_latent1 <- as.data.frame(final_model, par = c("lp_mixture_latent_1"))
matrix_of_draws_logprob_latent2 <- as.data.frame(final_model, par = c("lp_mixture_latent_2"))
for (i in 1:23){
  random_sample_more[i,] <- sample(seq(1:4000), number_of_draws_more)
  for (j in 1:number_of_draws_more){
    mixture_probability_latent1_more[i,j] <- prod(as.vector(exp(matrix_of_draws_logprob_latent1[,start_index[i]:end_index[i]])[random_sample_more[i, j],]))
    mixture_probability_latent2_more[i,j] <- prod(as.vector(exp(matrix_of_draws_logprob_latent2[,start_index[i]:end_index[i]])[random_sample_more[i, j],]))
    mixture_probability_more[i,j] <- matrix_of_draws_final$lambda[random_sample_more[i, j]]*mixture_probability_latent2_more[i,j]/(matrix_of_draws_final$lambda[random_sample_more[i, j]]*mixture_probability_latent2_more[i,j] + (1-matrix_of_draws_final$lambda[random_sample_more[i, j]])*mixture_probability_latent1_more[i,j])
    random_sample_lambda_more[i,j] <- rbinom(prob = mixture_probability_more[i,j], size = 1, n = 1)
  }
}
d_y23_initial <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y23_intermediate <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y23_final <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm))
d_y_standard_initial <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm_standard))
d_y_standard_intermediate <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm_standard))
d_y_standard_final <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm_standard))

for (i in 1:number_of_draws_more){
  d_y23_initial[i,] <- h_function_4para(matrix_of_draws_initial$phi_23[random_sample_more[23,i]]*d_zm,matrix_of_draws_initial$beta_1[random_sample_more[23,i]],matrix_of_draws_initial$beta_2[random_sample_more[23,i]],matrix_of_draws_initial$beta_3[random_sample_more[23,i]],matrix_of_draws_initial$beta_4[random_sample_more[23,i]])  
  d_y_standard_initial[i,] <- h_function_4para(125*d_zm_standard,matrix_of_draws_initial$beta_1[random_sample_more[23,i]],matrix_of_draws_initial$beta_2[random_sample_more[23,i]],matrix_of_draws_initial$beta_3[random_sample_more[23,i]],matrix_of_draws_initial$beta_4[random_sample_more[23,i]])
  d_y23_intermediate[i,] <- log(h_function_4para(matrix_of_draws_intermediate$phi_23[random_sample_more[23,i]]*d_zm,matrix_of_draws_intermediate$beta_1[random_sample_more[23,i]],matrix_of_draws_intermediate$beta_2[random_sample_more[23,i]],matrix_of_draws_intermediate$beta_3[random_sample_more[23,i]],matrix_of_draws_intermediate$beta_4[random_sample_more[23,i]]))  
  d_y_standard_intermediate[i,] <- log(h_function_4para(125*d_zm_standard,matrix_of_draws_intermediate$beta_1[random_sample_more[23,i]],matrix_of_draws_intermediate$beta_2[random_sample_more[23,i]],matrix_of_draws_intermediate$beta_3[random_sample_more[23,i]],matrix_of_draws_intermediate$beta_4[random_sample_more[23,i]]))
  d_y23_final[i,] <- (1-random_sample_lambda_more[23,i])*log(h_function_4para(matrix_of_draws_final$phi_23[random_sample_more[23,i]]*d_zm,matrix_of_draws_final$beta_1[random_sample_more[23,i]],matrix_of_draws_final$beta_2[random_sample_more[23,i]],matrix_of_draws_final$beta_3[random_sample_more[23,i]],matrix_of_draws_final$beta_4[random_sample_more[23,i]])) + 
    random_sample_lambda_more[23,i]*log(h_function_4para(matrix_of_draws_final$phi_23[random_sample_more[23,i]]*d_zm,matrix_of_draws_final$beta_1[random_sample_more[23,i]],matrix_of_draws_final$beta2_unknown_23[random_sample_more[23,i]],matrix_of_draws_final$beta_3[random_sample_more[23,i]],matrix_of_draws_final$beta_4[random_sample_more[23,i]]))
  d_y_standard_final[i,] <- log(h_function_4para(125*d_zm_standard,matrix_of_draws_final$beta_1[random_sample_more[23,i]],matrix_of_draws_final$beta_2[random_sample_more[23,i]],matrix_of_draws_final$beta_3[random_sample_more[23,i]],matrix_of_draws_final$beta_4[random_sample_more[23,i]]))
}

#Prepare lines in GGplot
d_y_standard_median_initial <- apply(d_y_standard_initial,2, median)
d_y_standard_up_initial <- apply(d_y_standard_initial,2, quantile, prob = 0.975)
d_y_standard_down_initial <- apply(d_y_standard_initial,2, quantile, prob = 0.025)
d_y23_median_initial <- apply(d_y23_initial,2, median)
d_y23_up_initial <- apply(d_y23_initial,2, quantile, prob = 0.975)
d_y23_down_initial <- apply(d_y23_initial,2, quantile, prob = 0.025)

d_y_standard_median_intermediate <- apply(d_y_standard_intermediate,2, median)
d_y_standard_up_intermediate <- apply(d_y_standard_intermediate,2, quantile, prob = 0.975)
d_y_standard_down_intermediate <- apply(d_y_standard_intermediate,2, quantile, prob = 0.025)
d_y23_median_intermediate <- apply(d_y23_intermediate,2, median)
d_y23_up_intermediate <- apply(d_y23_intermediate,2, quantile, prob = 0.975)
d_y23_down_intermediate <- apply(d_y23_intermediate,2, quantile, prob = 0.025)

d_y_standard_median_final <- apply(d_y_standard_intermediate,2, median)
d_y_standard_up_final <- apply(d_y_standard_intermediate,2, quantile, prob = 0.975)
d_y_standard_down_final <- apply(d_y_standard_final,2, quantile, prob = 0.025)
d_y23_median_final <- apply(d_y23_final,2, median)
d_y23_up_final <- apply(d_y23_final,2, quantile, prob = 0.975)
d_y23_down_final <- apply(d_y23_final,2, quantile, prob = 0.025)

### Generate plots
zoom_1_initial <-  ggplot()+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_median_initial), size = 1,color = "black", linetype = "solid")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_up_initial), size = 1,color = "gray", linetype = "dashed")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_down_initial), size = 1,color = "gray", linetype = "dashed")+
  geom_point(aes(x = d_standard, y = c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13])), size = 1.5)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("dilution")+
  ylab("y")+
  labs(title = "Standard sample using initial model")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size = 8, angle = 25),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(0,1.05), ylim = c(0, max(c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13])) *1.1))

zoom_1_intermediate <-  ggplot()+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_median_intermediate), size = 1,color = "black", linetype = "solid")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_up_intermediate), size = 1,color = "gray", linetype = "dashed")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_down_intermediate), size = 1,color = "gray", linetype = "dashed")+
  geom_point(aes(x = d_standard, y = y_standard), size = 1.5)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Standard sample using intermediate model")+
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

zoom_1_final <-  ggplot()+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_median_final), size = 1,color = "black", linetype = "solid")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_up_final), size = 1,color = "gray", linetype = "dashed")+
  geom_line(aes(x = d_zm_standard, y = d_y_standard_down_final), size = 1,color = "gray", linetype = "dashed")+
  geom_point(aes(x = d_standard, y = y_standard), size = 1.5)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Standard sample using final model")+
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

plot23_initial<- ggplot()+
  geom_line(aes(x = d_zm, y = d_y23_median_initial, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y23_up_initial, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y23_down_initial, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = derf$FI[80:82]),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
 scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("y")+
  labs(title = "Unk 23 using initial model")+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major  = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = 8, angle = 25),
        axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
        axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 2500))

plot23_intermediate <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y23_median_intermediate, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y23_up_intermediate, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y23_down_intermediate, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[80:82])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unk 23 using intermediate model")+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major  = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = 8, angle = 25),
        axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
        axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))


plot23_final <- ggplot()+
  geom_line(aes(x = d_zm, y = d_y23_median_final, linetype = "posterior median"), size = 1,color = "black")+ #A is just a label, not a type
  geom_line(aes(x = d_zm, y = d_y23_up_final, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_line(aes(x = d_zm, y = d_y23_down_final, linetype = "95% credible interval"), size = 1,color = "gray")+
  geom_point(aes(x = c(1/10, 1/100, 1/10000), y = log(derf$FI[80:82])),color = "black",size = 2)+
  scale_x_continuous(expand = c(0, 0),breaks = c(1/10000, 1/100, 1/10), labels = c("0.0001", 1/100, 1/10))+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(x = 0.05, y = 1.25,label = "1, 50% interval [1,1]"), size = 5)+
  scale_linetype_manual(legend_title, values = c("posterior median" = "solid", "95% credible interval" = "dashed"))+ 
  xlab("dilution")+
  ylab("logy")+
  labs(title = "Unk 23 using final model")+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major  = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = 8, angle = 25),
        axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
        axis.text.y = element_text(color = "grey20", size = 10),
  )+
  coord_cartesian(xlim = c(-0.001,0.11), ylim = c(0, 10))

ggarrange(zoom_1_initial, zoom_1_intermediate, zoom_1_final, plot23_initial, plot23_intermediate, plot23_final, ncol = 3, nrow = 2)
ggsave("work_flow_plot.png", width = 14, height = 8)