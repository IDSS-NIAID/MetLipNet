#' Categorize Data by Retention Time (RT) Windows with Overlapping Intervals
#'
#' This function applies a sliding window approach to categorize data based on 
#' retention time (RT) with overlapping windows. Each data point can belong to 
#' multiple windows based on the defined overlap.
#'
#' @param data A data frame containing retention time (RT) and mass-to-charge ratio (mz) columns.
#' @param rt_col A character string specifying the column name for retention time (RT). Default is `"RT"`.
#' @param mz_col A character string specifying the column name for mass-to-charge ratio (mz). Default is `"mz"`.
#' @param window_width A numeric value specifying the width of each RT window. Default is `2.5`.
#' @param overlap A numeric value specifying the overlap between consecutive windows. Default is `0.75`.
#'
#' @return A data frame with an additional column `Window`, indicating the assigned window(s) for each observation.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#' 
#' # Example data
#' data <- data.frame(
#'   RT = c(1.0, 2.3, 3.7, 5.1, 6.8, 8.0),
#'   mz = c(100, 150, 200, 250, 300, 350)
#' )
#' 
#' # Apply categorization
#' categorized_data <- categorize_by_rt_window(data, rt_col = "RT", window_width = 2.5, overlap = 0.75)
#' print(categorized_data)
#'
categorize_by_rt_window <- function(data, rt_col = "RT", mz_col = "mz", window_width = 2.5, overlap = 0.75) {
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  
  min_rt <- min(data[[rt_col]], na.rm = TRUE)
  max_rt <- max(data[[rt_col]], na.rm = TRUE)
  
  # Define the window breaks with overlap
  window_starts <- seq(min_rt, max_rt - window_width, by = window_width - overlap)
  window_ends <- window_starts + window_width
  
  # Create a mapping of each RT value to its respective windows
  expanded_data <- lapply(seq_along(window_starts), function(i) {
    window_data <- data %>%
      dplyr::filter(.data[[rt_col]] >= window_starts[i] & .data[[rt_col]] < window_ends[i]) %>%
      dplyr::mutate(Window = i)
    return(window_data)
  })
  
  # Combine results into a single data frame
  categorized_data <- dplyr::bind_rows(expanded_data)
  return(categorized_data)

}
