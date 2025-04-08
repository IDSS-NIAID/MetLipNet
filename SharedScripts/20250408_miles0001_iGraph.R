
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

dat <-read.csv("20250407_HighGroups_Only_For_iGraph.csv")

meta_data_names <- dat %>% 
  select(Class, Shape) %>% # specify meta data columns 
  names()

#### Correlations ####
# This function calculates the pairwise Pearson correlations for metabolite intensities 
# and provides the correlation estimate and p-value for each pair.

# view the help page
?cal_met_cor

# The returned data frame will be ready for the igraph function.
dat_cor_result <- MetLipNet::cal_met_cor(data = dat, 
                                         meta_cols = meta_data_names, 
                                         identifier_col = "Metabolite", 
                                         method = "pearson",
                                         cor_threshold = 0.5,
                                         p_threshold = 0.01)
#have to change to dat_scaled if want to use this data!!! Value raw to value scaled

#### Correlation Data Export ####
#write.csv(R_Data_Corr, file="Notarangelo0001 R vals pearson all v all.csv")#Option to write the R matrix to a csv. You can do the same thing for the p values for handling elsewhere if needed.
#write.csv(p_Data_Corr, file="Notarangelo0001 p vals pearson all v all.csv")



#### Data Preparation for igraph ####

clinical_covar<- c("bodyweightlbs (AVG)","bodyweightkg", "bmi (AVG)", "waistcircm (AVG)", "fatmasslbs (AVG)", "fatmasspercent", "fatfreemasslb", "fatfreemasspercent",
"wholebodymuscle AVG","torsomuscle", "visceraladipose (AVG)", "restingee", "Systolic", "Diastolic", "CHOL_0", "HDL_0", "TG_0", "ALT_0", "AST_0", "GLU_0", "LDL_0", "VLDL_0",
"INS_0", "VO2_pred_max", "TNF_0")

lipids<-as.vector(unlist(dat$Metabolite[1:867]))




# modify the correlation result data frame for attributes that apply to metabolite pairs
dat_igraph <- dat_cor_result %>% 
  mutate(edge_color = case_when(cor < 0 ~ "estimate<0",
                                cor > 0 ~ "estimate>0",
                                TRUE ~ "grey")) %>% 
  filter(from %in% clinical_covar| to %in% clinical_covar | to %in% lipids)

#### igraph Network Generation ####

# This creates the network object which can be modified from here
g <- igraph::graph_from_data_frame(dat_igraph, directed = FALSE)


# Remove nodes with less than 3 connections 
#gRemoveNodes <- which(degree(g)<3) 
#trim_g <- delete.vertices(g, gRemoveNodes)

#isoNodes <- which(degree(trim_g)==0)
#trim_g1 <- delete.vertices(trim_g, isoNodes)

center_node <- "TAG(54:5_FA18:2)"  # change this to a clinical variable, whatever you're interested 
radius <- 2 # degrees of separation (2 for 2nd-degree neighbors), anything beyonmd 2 might get confusing 

# Extract the ego graph
ego_nodes <- unlist(ego(g, order = radius, nodes = center_node, mode = "all"))
subg <- induced_subgraph(g, ego_nodes)



# set the seed for reproducibility
# defining the layout external to the plot function with the seed set ensures the resulting plot 
# will look the same every time it is regenerated
set.seed(3.14)
layout <- layout_with_fr(subg)

# Plot the ego graph
plot(subg,
     layout = layout,
     vertex.label = V(subg)$name,
     vertex.size = 15,
     vertex.color = ifelse(V(subg)$name == center_node, "red", "lightblue"),
     vertex.label.cex = 0.8,
     main = paste("Ego Network (radius =", radius, ") for:", center_node)) #use this one to trouble shoot, this will be far faster, only use gg graph to make nice plot



#### Optimize the network layout ####

# function that relies on bootsrap sampling and XGBoost to predict the best layout
#best_layout <- MetLipNet::optimize_network_layout(trim_g1, n_samples = 40, n_bootstrap = 50)


#### Data Preparation for igraph network visualization ####

attributes <- dat %>% 
  select(Metabolite, Class, Shape)

# modify node attributes
nodes <- tibble(Metabolite = V(subg)$name) %>% # change to subg// check for consistency in naming!!! 
  left_join(attributes, join_by(Metabolite))

V(subg)$Class <- nodes$Class
V(subg)$Shape <- nodes$Shape

best_layout <- MetLipNet::optimize_network_layout(subg, n_samples = 80, n_bootstrap = 200)

#### Plot the network ####

# generate the plot
p <- ggraph(subg, layout = best_layout) +
  geom_edge_link(aes(color = edge_color), alpha = 0.8) +
  geom_node_point(size = 6, aes(color = Class, shape = Shape)) +
  geom_node_text(aes(label = name), repel = TRUE) + #comment this out to remove node names 
  scale_edge_color_manual(values = c("black", "red"))+
  scale_color_brewer(palette = "Accent") # can change the color pallet here 
#ggtitle(plot_title)+
#labs(edge_color = "Legend"))

# view the plot
plot(p)

#install.packages("svglite", dependencies = TRUE, INSTALL_opts = '--no-lock')
# save the plot
#ggsave(plot = p, "milesm0001_test.svg", width = 6000, height = 4000, units = "px") # can open this in illustrator and make clusters the size you want 



