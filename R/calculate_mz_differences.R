#' Calculate m/z Differences within RT Windows and Detect Known Molecular Differences
#'
#' This function calculates the pairwise m/z differences within retention time (RT) windows
#' and detects known molecular differences based on predefined values.
#'
#' @param data A data frame containing m/z values and RT windows.
#' @param mz_col Character. The column name in `data` that contains m/z values. Default is "mz".
#' @param window_col Character. The column name in `data` that defines RT windows. Default is "Window".
#' @param selected_groups Character vector. Specific molecular differences to detect. Default is NULL (uses all known differences).
#'
#' @return A data frame containing pairwise m/z differences and detected molecular matches.
#'
#' @export
#' @import dplyr
calculate_mz_differences <- function(data, mz_col = "mz", window_col = "Window", selected_groups = NULL) {
  known_differences <- c(
    "17" = "Ammonia", "14" = "CH2", "28" = "Carbonyl", "58" = "Aldehyde", "60" = "Alcohol or Ether",
    "15" = "Methyl (-CH3)", "16" = "Hydroxyl or Amine", "45" = "Carboxyl (-COOH)",
    "32" = "Sulfur (-S)", "80" = "Phosphate (-HPO3)", "96" = "Sulfate (-SO4)",
    "91" = "Benzyl (-C6H5CH2)", "77" = "Phenyl (-C6H5)", "42" = "Acetyl (-COCH3)",
    "44" = "CO2",  "19" = "KYN-KYN Acid", "57" = "3OHKYN-Q Acid",
    "87" = "NAD-Q Acid", "19" = "3OHKYN-Xan Acid", "8" = "TRP-KYN"
  )
  
  if (is.null(selected_groups)) {
    selected_groups <- names(known_differences)
  }
  
  known_differences <- known_differences[names(known_differences) %in% as.character(selected_groups)]
  
  data %>%
    group_by(!!sym(window_col)) %>%
    group_modify(~ {
      window_data <- .x
      pairwise_diffs <- expand.grid(mz1 = window_data[[mz_col]], mz2 = window_data[[mz_col]]) %>%
        filter(mz1 > mz2) %>%
        mutate(
          difference = round(mz1 - mz2),
          molecule_match = case_when(
            as.character(difference) %in% names(known_differences) ~ known_differences[as.character(difference)],
            TRUE ~ NA_character_
          )
        ) %>%
        filter(!is.na(molecule_match))
      
      return(pairwise_diffs)
    })
}
