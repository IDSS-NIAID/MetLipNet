# Untargeted-Metabolomics-Network-Analysis
 Includes functions to handle the preparation of untargeted metabolomics data for network visualization using the established igraph package.


# RT Window Categorization with Overlapping Intervals

This repository provides a function to categorize data based on retention time (RT) using a sliding window approach with overlapping intervals.

## Installation

To install and use this function in your R environment:

```r
# Install dependencies if not already installed
install.packages(c("dplyr", "tidyr"))

# Source the function if not in an R package
source("path/to/categorize_by_rt_window.R")  # Adjust path as needed


library(dplyr)
library(tidyr)

# Generate example data
df <- data.frame(
  RT = runif(1000, min = 0.026, max = 13.9),
  mz = rnorm(1000, mean = 200, sd = 50) + rnorm(1000, mean = 0, sd = 10) * runif(1000, min = 0.5, max = 1.5)
)

# Apply categorization function
df_categorized <- categorize_by_rt_window(df, overlap = 0.1)

# Display first 5 rows per RT window
df_categorized %>%
  group_by(Window) %>%
  slice_head(n = 5) %>%
  print()
