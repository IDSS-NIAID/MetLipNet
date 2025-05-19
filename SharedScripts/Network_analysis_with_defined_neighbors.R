library(igraph)

# Collect data from mtcars
# The mtcars dataset is a built-in dataset in R that contains measurements on 11 different attributes for 32 different cars.
cars <- mtcars
car_names <- rownames(cars)
threshold <- 2

# Calculate similarity based on mpg
adj_matrix <- outer(cars$mpg, cars$mpg, function(x, y) abs(x - y) < threshold)
diag(adj_matrix) <- 0
colnames(adj_matrix) <- rownames(adj_matrix) <- car_names

# Make the similarity network graph
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Choose node(s) of interest
center_node <- "Mazda RX4"  # change this to your desired car
radius <- 2  # degrees of separation (2 for 2nd-degree neighbors)

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
<<<<<<< Updated upstream
     vertex.label = V(subg)$name,
=======
          vertex.label = V(subg)$name,
>>>>>>> Stashed changes
     vertex.size = 15,
     vertex.color = ifelse(V(subg)$name == center_node, "red", "lightblue"),
     vertex.label.cex = 0.8,
     main = paste("Ego Network (radius =", radius, ") for:", center_node))
