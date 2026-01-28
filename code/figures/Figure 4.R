library(tidyverse) 
library(ggpubr)
library(cowplot)
library(janitor)
figure4_data <- read_csv("Figure 4 Data.csv")

method_current_data <- figure4_data[-1, 1:4] %>% 
  setNames(c("id","q1","median","q3")) %>% 
  #row_to_names(row_number = 1) %>% 
  mutate(Method = "Bayesian contamination model",
         y_method = 1,
         id = str_remove(id,"Unk "),
         id = factor(id, levels = c(10,17,5,4,2,21,23, 9, 3,13,8,7,12,14,19,6,
                                    20,22,18,11,1,16,15)))

method_classical_data <- figure4_data[-1, c(1,8)] %>% 
  setNames(c("id","median")) %>% 
  #row_to_names(row_number = 1) %>% 
  mutate(Method = "Classical calibration",
         median = as.numeric(median),
         y_method = 3,
         id = str_remove(id,"Unk "),
         id = factor(id, levels = c(10,17,5,4,2,21, 23, 9, 3,13,8,7,12,14,19,6,
                                    20,22,18,11,1,16,15)))

### Add unknown 23 and 3 back
method_classical_data$median[method_classical_data$id == "3"]
method_classical_data$median[method_classical_data$id == "23"]
method_classical_data_new <-method_classical_data %>%
  mutate(median_new = ifelse(median<= 35, median, 35))
View(method_classical_data_new)
method_classical_data_new_part1 <- method_classical_data_new[1:10,]
method_classical_data_new_part2 <- method_classical_data_new[11:23,]
method_current_data_part1 <- method_current_data[1:10,]
method_current_data_part2 <- method_current_data[11:23,]

plot_part1<- ggplot()+
  geom_point(aes(x = id, y = as.numeric(median), shape = Method, color = Method),  data = method_current_data_part1)+
  geom_errorbar(aes(x = id, ymin = as.numeric(q1), ymax = as.numeric(q3)), width = 0.2,  data = method_current_data_part1)+
  geom_point(aes(x = id, y = as.numeric(median_new), shape = Method, color = Method),data = method_classical_data_new_part1)+
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 5, 10, 15, 20, 25, 30, 35))+
  scale_x_discrete(labels=c("10", "17*", "5", "4", "2", "21", "23*", "9*", "3*", "13"))+
  annotate(geom = "text", x = 7, y = 33, label = "247", hjust = "left")+
  annotate(geom = "text", x = 9, y = 33, label = "436", hjust = "left")+
  xlab("Sample ID")+
  ylab("Estimated Concentration ")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.text.x=element_text(size = 13),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 13),
    text = element_text(size = 13),
    legend.position="bottom"
  )+
  coord_cartesian( ylim = c(0,35.5))

plot_part2<- ggplot()+
  geom_point(aes(x = id, y = as.numeric(median), shape = Method, color = Method),  data = method_current_data_part2)+
  geom_errorbar(aes(x = id, ymin = as.numeric(q1), ymax = as.numeric(q3)), width = 0.2,  data = method_current_data_part2)+
  geom_point(aes(x = id, y = as.numeric(median_new), shape = Method, color = Method),data = method_classical_data_new_part2)+
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 1, 2, 3,4, 5))+
  scale_x_discrete(labels=c("8", "7", "12", "14", "19", "6*", "20", "22*", "18", "11", "1*", "16*", "15"))+
  xlab("Sample ID")+
  ylab("Estimated Concentration ")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major  = element_blank(),
    panel.border = element_blank(),
    axis.text.x=element_text(size = 13),
    axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
    axis.text.y = element_text(color = "grey20", size = 13),
    text = element_text(size = 13),
    legend.position="bottom"
  )+
  coord_cartesian( ylim = c(0,5))

ggarrange(plot_part1, plot_part2, widths = c(10, 13) , common.legend= TRUE, legend = "bottom",
                              labels = "AUTO")

ggsave("Figure4.png", width = 10, height = 4)



