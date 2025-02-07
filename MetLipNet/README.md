# Untargeted-Metabolomics-Network-Analysis
 Includes functions to handle the preparation of untargeted metabolomics data for network visualization using the established igraph package.


## RT Window Categorization with Overlapping Intervals

This function is used to categorize data based on retention time (RT) using a sliding window approach with overlapping intervals.

## Installation

To install and use this function in your R environment:

```r
# Install dependencies if not already installed
install.packages(c("dplyr", "tidyr"))


library(dplyr)
library(tidyr)

# Generate example data
df <- data.frame(
  RT = runif(1000, min = 0.026, max = 13.9),
  mz = rnorm(1000, mean = 200, sd = 50) + rnorm(1000, mean = 0, sd = 10) * runif(1000, min = 0.5, max = 1.5)
)

# Apply categorization function with adjusted window overlap
df_categorized <- categorize_by_rt_window(df, overlap = 0.5)

# Display first 3 rows per RT window
df_categorized %>%
  group_by(Window) %>%
  slice_head(n = 3) %>%
  print()


## Expected Output

The function will categorize data into overlapping retention time (RT) windows. Below is an example of the output:

# A tibble: 3 × 3 per Window
     RT    mz Window
  <dbl> <dbl>  <int>
1  0.50  201.     1
2  1.23  189.     1
3  1.87  210.     1

4  2.78  205.     2
5  3.12  198.     2
6  3.98  215.     2

7  4.87  220.     3
8  5.33  190.     3
9  5.79  195.     3
...

