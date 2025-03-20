#' Calculate Pairwise Correlations for Metabolite Intensities
#'
#' This function calculates the metabolite pairwise correlations 
#' and provides the correlation estimate and p-value for each pair.
#'
#' @param data A data frame containing metabolite names, intensity values, and metadata columns.
#' @param identifier_col Character. The column name that uniquely identifies metabolites. Default is "metabolite".
#' @param method Character. The correlation method to use. Options: "pearson", "spearman", "kendall". Default is "pearson".
#' @param cor_threshold Numeric. The absolute correlation estimate threshold for filtering results. Default is NULL (no filtering).
#' @param p_threshold Numeric. The p-value threshold for filtering results. Default is NULL (no filtering).
#'
#' @return A tibble with columns:
#'   \item{from}{Character. Name of the first metabolite in the correlation pair.}
#'   \item{to}{Character. Name of the second metabolite in the correlation pair.}
#'   \item{cor}{Numeric. Correlation coefficient.}
#'   \item{p}{Numeric. P-value of the correlation test.}
#'
#' @examples
#'
#' @export
cal_met_cor <- function(data, meta_cols = NULL, identifier_col = "Metabolite", method = "pearson", cor_threshold = NULL, p_threshold = NULL) {

  # Ensure identifier_col exists before transformation
  if (!(identifier_col %in% names(data))) {
    stop(paste("Identifier column", identifier_col, "not found in dataset"))
  }
  
  # Pre-process data
  data_clean <-data %>% 
    select(-any_of(meta_data_names)) %>% 
    t() %>% 
    janitor::row_to_names(row_number = 1)
  
  # Run the correlation
  data_corr_result <- rcorr(data_clean, type = method)
  
  # Flatten the correlation results
  data_corr_flat <- cor_gather(data_corr_result) %>% 
    rename(from = var1, to = var2)
  
  # Apply correlation estimate and p-value thresholds
  isl_out <- data_corr_flat %>% 
    filter(abs(cor)>cor_threshold & p<p_threshold)
  
  print(paste("Snap"))
  
  return(isl_out)
}
