

library(ppcor)
library(reshape2)

# define covariates
co_vars <- c("gestational age", "Maternal age", "ethnicity", "gravida")

# define all other variables to correlate
corr_vars <- setdiff(names(dat), co_vars)

# compute partial correlation matrix

            #\- correlation varibles -/ \- adjustment variables -/
pcor_out <- pcor(dat[ , corr_vars], dat[ , co_vars])

# extract partial correlation values
pcor_matrix <- pcor_out$estimate

# convert to tidy format
pcor_long <- melt(pcor_matrix, varnames = c("from", "to"), value.name = "partial_cor")

# remove duplicates and self-correlations
pcor_long <- pcor_long[pcor_long$from != pcor_long$to, ]
pcor_long <- pcor_long[!duplicated(t(apply(pcor_long[ ,1:2], 1, sort))), ]

# add p-values
pval_matrix <- pcor_out$p.value
pval_long <- melt(pval_matrix, varnames = c("from", "to"), value.name = "p_value")

# merge p-values and partial correlations
pcor_long <- merge(pcor_long, pval_long, by = c("from", "to"))