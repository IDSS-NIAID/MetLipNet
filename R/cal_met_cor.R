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
#'
#' @return A data frame containing pairwise Pearson correlation estimates and p-values.
#'
#' @examples
#' library(dplyr)
#' library(rstatix)
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
cal_met_cor <- function(data, intensity_col = "intensity", identifier_col = "metabolite", method = "pearson", remove_na = TRUE, p_threshold = NULL) {
  if (nrow(data) < 2) {
    return(tibble(metabolite1 = character(), metabolite2 = character(), estimate = numeric(), p_value = numeric()))
  }
  
  # Identify metadata columns dynamically
  metadata_cols <- setdiff(names(data), c(intensity_col, identifier_col))
  
  # Prepare intensity data in wide format for correlation
  intensity_wide <- data %>%
    pivot_wider(names_from = all_of(metadata_cols), values_from = all_of(intensity_col))
  
  # Remove rows with insufficient data for correlation if specified
  if (remove_na) {
    intensity_wide <- intensity_wide %>% drop_na()
  }
  
  if (ncol(intensity_wide) <= 2) {
    return(tibble(metabolite1 = character(), metabolite2 = character(), estimate = numeric(), p_value = numeric()))
  }
  
  # Compute correlation matrix
  cor_matrix <- cor_mat(intensity_wide %>% select(-all_of(identifier_col)), method = method)
  
  # Extract p-values separately
  p_matrix <- cor_test(intensity_wide %>% select(-all_of(identifier_col)), method = method) %>%
    select(var1, var2, p) %>%
    rename(metabolite1 = var1, metabolite2 = var2, p_value = p)
  
  # Convert correlation matrix into long format
  cor_data <- cor_matrix %>%
    cor_gather() %>%
    rename(metabolite1 = var1, metabolite2 = var2, estimate = cor) %>%
    filter(metabolite1 > metabolite2)
  
  # Join correlation estimates and p-values
  result <- left_join(cor_data, p_matrix, by = c("metabolite1", "metabolite2")) %>% 
    select(-p_value)
  
  # Apply p-value threshold if specified
  if (!is.null(p_threshold)) {
    result <- result %>% filter(p_value <= p_threshold)
  }
  
  return(result)
}