library(tidyverse) 
library(ggpubr)
library(cowplot)
derf <- read_csv("derf.csv")
standard_concentration <- 125
y_standard = c(derf$FI[1:13] - derf$SD[1:13], derf$FI[1:13] + derf$SD[1:13])
d_standard = c(0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048,
               0, 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048)
#standard_df <- cbind(y_standard, d_standard*standard_concentration)
#colnames(standard_df) <- c("obs", "concentration")
standard_id <- c(rep("Standard", length(y_standard)))
standard_df <- cbind(y_standard, d_standard*standard_concentration, standard_id)
colnames(standard_df) <- c("obs", "concentration","Sample")
standard_df <- as.data.frame(standard_df)
ggplot()+
  geom_point(aes(x = concentration, y = obs), data = standard_df)+
  geom_smooth(method = "lm",aes(x = concentration, y = obs), data = standard_df)
#if directly apply geom_smooth using lm, loess or spline, will not produce the required curve
#so we manually generate the 4-para logistic mean curve
h_function_4para <- function(x, beta_1, beta_2, beta_3, beta_4) {
  beta_1 + beta_2/(1 + (x/beta_3)^(-beta_4))
}
d_zm_standard = as.vector(seq(0,1,0.0001)) 
standard_unknown <- readRDS("standard_unknown_4para_model1_mixture_proportion_more_diffuse.rds")  
matrix_of_draws <- as.data.frame(standard_unknown, par = c("beta"))
set.seed(1234)
number_of_draws_more <- 4000
random_sample_more <- matrix(NA, nrow = 23, ncol = number_of_draws_more)
colnames(matrix_of_draws) <- c("beta_1","beta_2","beta_3","beta_4")
for (i in 1:23){
  random_sample_more[i,] <- sample(seq(1:4000), number_of_draws_more)
}
d_y_standard_more <- matrix(NA, nrow = number_of_draws_more, ncol = length(d_zm_standard))

for (i in 1:number_of_draws_more){
  d_y_standard_more[i,] <- log(h_function_4para(standard_concentration*d_zm_standard,matrix_of_draws$beta_1[random_sample_more[23,i]],matrix_of_draws$beta_2[random_sample_more[23,i]],matrix_of_draws$beta_3[random_sample_more[23,i]],matrix_of_draws$beta_4[random_sample_more[23,i]]))
}
d_y_standard_median <- apply(d_y_standard_more,2, median)
d_y_standard_up <- apply(d_y_standard_more,2, quantile, prob = 0.975)
d_y_standard_down <- apply(d_y_standard_more,2, quantile, prob = 0.025)

standard_df$obs <- as.numeric(standard_df$obs)
standard_df$concentration <- as.numeric(standard_df$concentration)
#Regenerate plot on  log scale
ggplot()+
  geom_point(aes(x = concentration, y = log(obs)), data = standard_df)+
  geom_line(aes(x = d_zm_standard*standard_concentration, y = d_y_standard_median), size = 1,color = "black", linetype = "solid")

#add one unknown sample
unknown_3_concnetration <- (6.92 + 68.04 + 1233.87)/3
unknown_3_y <- c(660, 646, 74)
unknown_3_x <- unknown_3_concnetration*c(1/10, 1/100, 1/10000)
unknown_3_df <- as.data.frame(cbind(unknown_3_y, unknown_3_x))
colnames(unknown_3_df) <- c("obs", "concentration")

ggplot()+
  geom_point(aes(x = concentration, y = log(obs)), data = standard_df)+
  geom_line(aes(x = d_zm_standard*standard_concentration, y = d_y_standard_median), size = 1,color = "black", linetype = "solid")+
  geom_point(aes(x = concentration, y = log(obs)), data = unknown_3_df, shape = 24)

#add multiple unknown samples
unknown_23_concnetration <- (15.87 + 78.78 + 647.51)/3
unknown_23_y <- c(1904, 783, 39)
unknown_23_x <- unknown_23_concnetration*c(1/10, 1/100, 1/10000)
unknown_10_concnetration <- (32.47 + 32.21)/2
unknown_10_y <- c(4259.5, 241, 6)
unknown_10_x <- unknown_10_concnetration*c(1/10, 1/100, 1/10000)

sample_id <- c(rep("Unk 3", 3), rep("Unk 23", 3), rep("Unk 10", 3))
unknown_data_combined_df <- cbind(c(unknown_3_y, unknown_23_y, unknown_10_y), c(unknown_3_x, unknown_23_x, unknown_10_x), sample_id)
unknown_data_combined_df <- as.data.frame(unknown_data_combined_df)
colnames(unknown_data_combined_df) <- c("obs", "concentration", "Sample")
str(unknown_data_combined_df)

Figure1_nologx <-ggplot()+
  geom_line(aes(x = d_zm_standard*standard_concentration, y = d_y_standard_median), size = 0.3,color = "black", linetype = "solid")+
  geom_point(aes(x = as.numeric(concentration), y = log(as.numeric(obs)), shape = Sample, color = Sample), data = standard_df)+
  geom_point(aes(x = as.numeric(concentration), y = log(as.numeric(obs)), shape = Sample, color = Sample), data = unknown_data_combined_df)+
  geom_line(aes(x = as.numeric(concentration), y = log(as.numeric(obs)), color = Sample), data = unknown_data_combined_df)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("Concentration")+
  ylab("Log Observed Signal")+
  scale_color_manual(values=c("black","#FF0000", "#E69F00", "#56B4E9"))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.text.x=element_text(size = 15),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 15),
    text = element_text(size = 15),
    legend.position = "none"
  )+
  coord_cartesian(xlim = c(0,75), ylim = c(0,10.5))

Figure1_logx <- ggplot()+
  geom_line(aes(x = log(d_zm_standard*standard_concentration), y = d_y_standard_median), size = 0.3,color = "black", linetype = "solid")+
  geom_point(aes(x = log(as.numeric(concentration)), y = log(as.numeric(obs)), shape = Sample, color = Sample), data = standard_df)+
  geom_point(aes(x = log(as.numeric(concentration)), y = log(as.numeric(obs)), shape = Sample, color = Sample), data = unknown_data_combined_df)+
  geom_line(aes(x = log(as.numeric(concentration)), y = log(as.numeric(obs)), color = Sample), data = unknown_data_combined_df)+
  scale_x_continuous(expand = c(0, 0), breaks = c(-4, 0, 4, 8))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("Log Concentration")+
  ylab("Log Observed Signal")+
  scale_color_manual(values=c("black","#FF0000", "#E69F00", "#56B4E9"))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.text.x=element_text(size = 15),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 15),
    text = element_text(size = 15),
    legend.position = c(0.8, 0.4)
  )+
  coord_cartesian(xlim = c(-4,10), ylim = c(0,10.5))


ggarrange(Figure1_nologx, Figure1_logx)
ggsave("Figure2_updated_reviewer_commented.png", width = 8, height = 3)
#in the log concentration plot, if we choose some 'normal' sample, their concentration/10000 will be no way above exp(-5), so there
#will be one missing point
#we could always change the color scale maunally
