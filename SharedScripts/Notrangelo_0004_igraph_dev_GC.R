
#### Package Installation ####

# run this line every time there is an update to the MetLipNet package
devtools::install_github("IDSS-NIAID/MetLipNet", force = TRUE)



# install dependent packages
install.packages("igraph")
install.packages("dplyr")
install.packages("corrr")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("rstatix")
install.packages("corrr")



#### Required Libraries ####

# run this chunk to open libraries containing functions to execute the script
library(MetLipNet)
library(future)
library(ggraph)
library(openxlsx)
library(dplyr)
library(future)
library(furrr)
library(RColorBrewer)
library(tidyverse)
library(corrr)
library(igraph)
library(ggplot2)
library(rstatix)
library(dplyr)


#### Read in Data ####

dat <-read.csv("Gwen_GiGi_CleanedData.csv")

dat_longer <- dat %>% 
  tidyr::pivot_longer(cols = 4:ncol(dat), names_to = "condition", values_to = "Value_raw") %>% 
  filter(Class != "LM") %>% 
  select(-Shape, -Class) 

#### Remove other omics ####
remove_me <- c("Bile", "Clinical", "Estrogen")

dat_longer_trim <- dat %>% 
  tidyr::pivot_longer(cols = 4:ncol(dat), names_to = "condition", values_to = "Value_raw") %>% 
  filter(!Class %in% remove_me) %>% 
  select(-Shape, -Class) 


#### Data Scaling ####

dat_scaled <- dat %>% 
  tidyr::pivot_longer(cols = 4:ncol(dat), names_to = "condition", values_to = "Value_raw") %>% #column that's being correlated needs to be named condition
  filter(Class != "LM") %>% 
  group_by(Class, Metabolite) %>% 
  mutate(Value_scaled = as.numeric(scale(Value_raw, center=TRUE, scale=TRUE))) %>% 
  ungroup() %>% 
  select(-Shape, -Class, -Value_raw) 


#### Correlations ####
# This function calculates the pairwise Pearson correlations for metabolite intensities 
# and provides the correlation estimate and p-value for each pair.

# view the help page
?cal_met_cor

# The returned data frame will be ready for the igraph function.
dat_cor_result <- MetLipNet::cal_met_cor(dat_scaled, intensity_col = "Value_scaled", identifier_col = "Metabolite", max_workers = 13) #have to change to dat_scaled if want to use this data!!! Value raw to value scaled

#### Correlation Data Export ####
#write.csv(R_Data_Corr, file="Notarangelo0001 R vals pearson all v all.csv")#Option to write the R matrix to a csv. You can do the same thing for the p values for handling elsewhere if needed.
#write.csv(p_Data_Corr, file="Notarangelo0001 p vals pearson all v all.csv")



#### Data Preparation for igraph ####

# modify the correlation result data frame for attributes that apply to metabolite pairs
dat_igraph <- dat_cor_result %>% 
  filter(p_value < 0.01) %>% 
  filter(estimate > 0.7) %>%
  mutate(edge_color = case_when(estimate < 0 ~ "estimate<0",
                                estimate > 0 ~ "estimate>0",
                                TRUE ~ "grey"))


#### igraph Network Generation ####

# This creates the network object which can be modified from here
g <- igraph::graph_from_data_frame(dat_igraph, directed = FALSE)


# Remove nodes with less than 3 connections 
gRemoveNodes <- which(degree(g)<3) 
trim_g <- delete.vertices(g, gRemoveNodes)

isoNodes <- which(degree(trim_g)==0)
trim_g1 <- delete.vertices(trim_g, isoNodes)


#### Optimize the network layout ####

# function that relies on bootsrap sampling and XGBoost to predict the best layout
best_layout <- MetLipNet::optimize_network_layout(trim_g1, n_samples = 50, n_bootstrap = 100)


#### Data Preparation for igraph network visualization ####

attributes <- dat %>% 
  select(Metabolite, Class, Shape)

# modify node attributes
nodes <- tibble(Metabolite = V(trim_g1)$name) %>% 
  left_join(attributes, join_by(Metabolite))

V(trim_g1)$Class <- nodes$Class
V(trim_g1)$Shape <- nodes$Shape


#### Plot the network ####

# generate the plot
p <- ggraph(trim_g1, layout = best_layout) +
  geom_edge_link(aes(color = edge_color), alpha = 0.8) +
  geom_node_point(size = 6, aes(color = Class, shape = Shape)) +
  geom_node_text(aes(label = name), repel = TRUE) + #comment this out to remove node names 
  scale_edge_color_manual(values = c("black", "red"))+
  scale_color_brewer(palette = "Accent") # can change the color pallet here 
#ggtitle(plot_title)+
#labs(edge_color = "Legend"))

# view the plot
plot(p)

# save the plot
ggsave(plot = p, "Test_Network_Plot_p01.svg", width = 6000, height = 4000, units = "px") # can open this in illustrator and make clusters the size you want 



