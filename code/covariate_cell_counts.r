##### Cell counts  #####
# LA Salas et al. (2022) Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling. Nat Commun 13, 761. doi: https://doi.org/10.1038/s41467-021-27864-7

library(EpiDISH)
data(cent12CT.m)

cc <- epidish(beta.m = beta_val, ref.m = cent12CT.m, method = c("RPC")) # beta_val is the matrix of beta values for all CpG sites in your cohort, rather than just the target CpG sites
cell_counts <- as.data.frame(cc$estF)
cell_counts <- data.frame(Sample_Name = row.names(cell_counts), cell_counts)

head(cell_counts, n = 1)

# check for NAs
cc_na <- table(apply(cell_counts, 1, anyNA))
message(sprintf("Cell counts were successfully predicted for %s individuals.", cc_na[1]))
if (length(cc_na) > 1) {
    message(sprintf("Cell counts contain NAs for for %s individuals. Please check the input data.", cc_na[2]))
}

# save
write.csv(cell_counts, file = "cells.csv", row.names = FALSE, quote = FALSE)
