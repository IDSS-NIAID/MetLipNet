#' Generate, Save, and Plot Networks for Each RT Window
#'
#' This function generates network layouts for each RT window, saves the network data to an Excel file,
#' and produces network plots that are saved as JPEG files.
#'
#' @param data A data frame containing network data with at least the columns: `Window`, `mz1`, `mz2`, and `molecule_match`.
#' @param output_xlsx Character. The name of the output Excel file where network data will be saved. Default is "network_data.xlsx".
#' @param output_dir Character. The directory where network plots will be saved. Default is "network_plots".
#'
#' @return No return value. Outputs an Excel file with network data and saves network plots as JPEG images.
#'
#' @examples
#' library(igraph)
#' library(ggraph)
#' sample_data <- data.frame(
#'   Window = rep(1:3, each = 5),
#'   mz1 = runif(15, 100, 500),
#'   mz2 = runif(15, 100, 500),
#'   molecule_match = sample(c("A", "B", "C"), 15, replace = TRUE)
#' )
#' net_viz(sample_data)
#'
#' @export
net_viz <- function(data, output_xlsx = "network_data.xlsx", output_dir = "network_plots") {
  dir.create(output_dir, showWarnings = FALSE)
  wb <- createWorkbook()
  
  unique_windows <- unique(data$Window)
  
  for (window in unique_windows) {
    window_data <- filter(data, Window == window)
    
    # Save network data to Excel
    addWorksheet(wb, paste0("RT_Window_", window))
    writeData(wb, sheet = paste0("RT_Window_", window), window_data)
    
    # Prepare data for graph
    graph_data <- data.frame(from = round(window_data$mz1, 0), to = round(window_data$mz2, 0), edge_color = window_data$molecule_match)
    g <- graph_from_data_frame(graph_data, directed = TRUE)
    
    # Optimize layout for each window
    best_layout <- optimize_network_layout(g)
    
    plot_title <- paste("Network for RT Window", window)
    
    p <- ggraph(g, layout = best_layout) +
      geom_edge_link(aes(color = edge_color), alpha = 0.8) +
      geom_node_point(size = 5, color = "black") +
      geom_node_text(aes(label = name), repel = TRUE) +
      theme_minimal() +
      ggtitle(plot_title)
    
    plot(p)
    
    # Save plot as JPEG
    ggsave(filename = file.path(output_dir, paste0("network_RT_Window_", window, ".jpeg")), plot = p, width = 30, height = 20)
  }
  
  openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)
}
