# MetLipNet
 Functions to handle the preparation of untargeted metabolomics data for network analysis.


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

# Apply function to calculate m/z differences and match to pre-determined chemical groups
df_mz_differences_filtered <- calculate_mz_differences(df_categorized, selected_groups = c("17", "14", "28"))

head(df_mz_differences_filtered)


## Expected Output

The `categorize_by_rt_window` function will categorize data into overlapping retention time (RT) windows. 
Below is an example of the output:

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



The `calculate_mz_differences` function calculates pairwise m/z differences within retention time windows,
and matches the differences to selected chemical groups. 
Below is an example of the output:

# A tibble: 6 × 5
# Groups:   Window [1]
  Window   mz1   mz2 difference molecule_match
   <int> <dbl> <dbl>      <dbl> <chr>         
1      1  225.  211.         14 CH2           
2      1  239.  211.         28 Carbonyl      
3      1  164.  147.         17 Ammonia       
4      1  165.  147.         17 Ammonia       
5      1  161.  147.         14 CH2           
6      1  183.  166.         17 Ammonia   

