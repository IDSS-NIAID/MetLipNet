#### Required Libraries ####
library(tidyverse)

#### Read in Data ####

dat <-read.csv("Gwen_GiGi_CleanedData.csv")

#### Data Scaling ####

dat_scaled <- dat %>% 
  tidyr::pivot_longer(cols = 4:ncol(dat), names_to = "condition", values_to = "Value_raw") %>% #column that's being correlated needs to be named condition
  #filter(Class != "LM") %>% 
  group_by(Class, Metabolite) %>% 
  mutate(Value_scaled = as.numeric(scale(Value_raw, center=TRUE, scale=TRUE))) %>% 
  ungroup() 



##### Data Density Plots #####

###### Pre Scaling #######
dat_density <- dat_scaled 


density_plot_pre <- ggplot(dat_density, aes(x = Value_raw, fill = Class))+
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "YlGnBu") +
  
  scale_x_continuous(limits = c(0, 1e+02))+
  scale_y_continuous(expand = c(0,0))+
  
  labs(title = "Un-adjusted Data", 
       subtitle = paste("Sample type: Cell"),
       y = "density", 
       fill = "Experiment")+
  
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA), 
        panel.grid = element_blank())

plot(density_plot_pre)

ggsave(plot = density_plot_pre, 
       filename = paste(Sys.Date(), "pre_scaling_density_plot.svg", sep = "_"), 
       width = 3, height = 3, units = "in", dpi = 300)


###### Post Scaling #######

density_plot_post <- ggplot(dat_density, aes(x = Value_scaled, fill = Class))+
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "YlGnBu") +
  
  scale_x_continuous(limits = c(-4, 5))+
  scale_y_continuous(expand = c(0,0))+
  
  labs(title = "Scaled Data", 
       subtitle = paste("Sample type: Cell"),
       y = "density", 
       fill = "Experiment")+
  
  
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5), 
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA), 
    panel.grid = element_blank())

plot(density_plot_post)


ggsave(plot = density_plot_post, 
       filename = paste(Sys.Date(), "post_scaling_density_plot.svg", sep = "_"), 
       width = 4.5, height = 3, units = "in", dpi = 300)



