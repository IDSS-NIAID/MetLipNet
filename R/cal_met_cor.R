#' Calculate Pairwise Correlations for Metabolite Intensities
#'
#' This function calculates the metabolite pairwise correlations 
#' and provides the correlation estimate and p-value for each pair.
#'
#' @param data A data frame containing metabolite names, intensity values, and metadata columns.
#' @param intensity_col Character. The column name in `data` that contains intensity values. Default is "intensity".
#' @param identifier_col Character. The column name that uniquely identifies metabolites. Default is "metabolite".
#' @param method Character. The correlation method to use. Options: "pearson", "spearman", "kendall". Default is "pearson".
#' @param remove_na Logical. Whether to drop rows with NA values before correlation calculation. Default is TRUE.
#' @param p_threshold Numeric. The p-value threshold for filtering results. Default is NULL (no filtering).
#' @param max_workers Integer. The maximum number of parallel workers to use. Default is 4.
#'
#' @return A tibble with columns:
#'   \item{from}{Character. Name of the first metabolite in the correlation pair.}
#'   \item{to}{Character. Name of the second metabolite in the correlation pair.}
#'   \item{estimate}{Numeric. Correlation coefficient.}
#'   \item{p_value}{Numeric. P-value of the correlation test.}
#'
#' @examples
#' library(dplyr)
#' library(rstatix)
#' library(corrr)
#' set.seed(123)
#' df <- expand.grid(
#'   metabolite = paste0("M", 1:50),
#'   condition = paste0("C", 1:3)
#' ) %>%
#'   mutate(intensity = rnorm(n(), mean = 1000, sd = 200))
#' 
#' cal_met_cor(df)
#'
#' @export
cal_met_cor <- function(data, intensity_col = "intensity", identifier_col = "metabolite", method = "pearson", remove_na = TRUE, p_threshold = NULL, max_workers = 4) {
  if (nrow(data) < 2) {
    return(tibble(from = character(), to = character(), estimate = numeric(), p_value = numeric()))
  }
  
  # Identify metadata columns dynamically
  metadata_cols <- setdiff(names(data), c(intensity_col, identifier_col))
  
  # Ensure identifier_col exists before transformation
  if (!(identifier_col %in% names(data))) {
    stop(paste("Identifier column", identifier_col, "not found in dataset"))
  }
  
  # Pivot data so that metabolites are in columns and samples are in rows
  intensity_wide <- data %>%
    pivot_wider(names_from = identifier_col, values_from = intensity_col)
  
  # Remove rows with insufficient data for correlation if specified
  if (remove_na) {
    intensity_wide <- intensity_wide %>% drop_na()
  }
  
  # Ensure at least two metabolites exist for correlation
  if (ncol(intensity_wide) < 3) {  # One column will be sample IDs, so we need at least two metabolite columns
    warning("Not enough metabolites for correlation. Returning empty data frame.")
    return(tibble(from = character(), to = character(), estimate = numeric(), p_value = numeric()))
  }
  
  # Filter only numeric columns for correlation
  numeric_cols <- intensity_wide %>% select(where(is.numeric))
  
  available_workers <- min(availableCores() - 1, max_workers)
  if (available_workers < 1) available_workers <- 1
  
  plan(multisession, workers = available_workers)
  
  col_pairs <- combn(ncol(numeric_cols), 2, simplify = FALSE)
  
  cor_results <- future_map_dfr(col_pairs, function(pair) {
    cor_test_result <- suppressMessages(suppressWarnings(
      cor_test(data = numeric_cols, vars = colnames(numeric_cols)[pair[1]], vars2 = colnames(numeric_cols)[pair[2]], method = method)
    ))
  
  
  # Extract p-values separately using rstatix::cor_test()
  p_matrix <- rstatix::cor_test(numeric_cols, method = method)
  if (!"var1" %in% names(p_matrix) || !"var2" %in% names(p_matrix)) {
    stop("Unexpected column names in cor_test output")
  }
  
  p_matrix <- p_matrix %>%
    select(var1, var2, p) %>%
    rename(from = var1, to = var2, p_value = p)
  
  # Convert correlation matrix into long format
  cor_data <- cor_matrix %>%
    corrr::stretch() %>%
    rename(from = x, to = y, estimate = r) %>%
    filter(from > to)
  
  # Join correlation estimates and p-values
  result <- left_join(cor_data, p_matrix, by = c("from", "to"))
  
  # Apply p-value threshold if specified
  if (!is.null(p_threshold)) {
    result <- result %>% filter(p_value <= p_threshold)
  }
  
  return(result)
}