#' Optimize Network Layout using XGBoost with Bootstrapping and Validation
#'
#' This function optimizes the layout of a given igraph using multiple layout algorithms.
#' It employs bootstrapping for stress value estimation and uses XGBoost to predict the best layout.
#'
#' @param g An igraph object representing the network.
#' @param n_samples Integer. Number of layout samples to test. Default is 50.
#' @param n_bootstrap Integer. Number of bootstrap iterations for stress calculation. Default is 10.
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#' @param seed An optional seed for reproducibility. Set to `NULL` to disable (default: 31415926)
#' 
#' @return A matrix representing the optimized layout of the graph.
#'
#' @details 
#' This function evaluates different network layouts (`"fr"`, `"kk"`, `"lgl"`) in parallel
#' and selects the one with the lowest stress value, calculated using bootstrapped sampling.
#' XGBoost is trained to predict the best layout based on stress values.
#'
#' Users should set up parallel processing manually by calling:
#' \code{future::plan(multisession, workers = availableCores())} before executing this function.
#' 
#' @examples
#' \dontrun{
#'   library(igraph)
#'   g <- erdos.renyi.game(100, p = 0.05)  # Generate a random network
#'   best_layout <- optimize_network_layout(g)
#'   plot(g, layout = best_layout)
#' }
#'
#' @export
optimize_network_layout <- function(g, n_samples = 50, n_bootstrap = 10, verbose = TRUE, seed = 31415926) {
  layouts <- list(
    "fr" = igraph::layout_with_fr,
    "kk" = igraph::layout_with_kk,
    "lgl" = function(g) {
      if (!igraph::is_connected(g)) {
        message("Graph is disconnected. Falling back to Kamada-Kawai layout.")
        return(igraph::layout_with_kk(g))
      } else {
        return(igraph::layout_with_lgl(g, maxiter = 100))
      }
    }
  )
  
  parameters <- tibble(
    layout = sample(names(layouts), n_samples, replace = TRUE)
  ) %>% split(1:n_samples)
  
  total_samples <- length(parameters)
  plan(multisession)
  results <- future_map(seq_along(parameters), function(i) {
    param <- parameters[[i]]
    layout_fn <- layouts[[param$layout]]
    
    if (verbose) {
      message(sprintf("[Optimization %d/%d] Running %s layout", i, total_samples, param$layout))
    }
    
    g_layout <- layout_fn(g)
    
    # Bootstrapping stress values
    stress_values <- replicate(n_bootstrap, {
      sampled_nodes <- sample(1:nrow(g_layout), size = round(0.8 * nrow(g_layout)), replace = TRUE)
      sum(dist(g_layout[sampled_nodes, ])^2)
    })
    
    stress_value <- mean(stress_values, na.rm = TRUE)
    if (is.nan(stress_value) || is.infinite(stress_value)) stress_value <- 0
    
    # Validate stress value before training XGBoost
    if (length(unique(parameters$layout)) < 2 || stress_value == 0) {
      return(list(layout = g_layout, stress = stress_value, predicted_layout = param$layout))
    }
    
    # Train XGBoost model only on valid values
    train_data <- data.frame(layout = as.numeric(factor(param$layout)), stress = stress_value)
    dtrain <- xgb.DMatrix(data = as.matrix(train_data$stress), label = train_data$layout)
    
    xgb_model <- xgboost(data = dtrain, max.depth = 2, eta = 1, nrounds = 10, objective = "multi:softmax", num_class = length(layouts), verbose = 0)
    predicted_layout <- predict(xgb_model, dtrain)
    
    if (verbose) {
      message(sprintf("[Optimization %d/%d] Completed %s layout | Stress: %.2f | Predicted Best: %s\n", i, total_samples, param$layout, stress_value, names(layouts)[predicted_layout + 1]))
    }
    
    list(layout = g_layout, stress = stress_value, predicted_layout = names(layouts)[predicted_layout + 1])
  }, .options = furrr_options(seed = 31415926))
  plan(sequential)
  
  results_df <- bind_rows(results)
  best_result <- results[[which.min(sapply(results, function(x) x$stress))]]
  return(best_result$layout)
}
