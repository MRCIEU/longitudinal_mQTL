library(tidyverse)

##########################################################
# To identify participants with methylation data across multiple time points
##########################################################
time_num <- samples %>%
    group_by(ind_id) %>%
    summarise(num = n_distinct(time_point)) # count how many time point identifiers (samples) each individual has
samples$time_num <- time_num$num[match(samples$ind_id, time_num$ind_id)]
samples <- samples[!(samples$time_num == 1), ] # rule out individuals with only one sample


##########################################################
# Based on the available quality control data within the cohort, individuals will be excluded if discrepancies are identified, such as incorrect designated sex, or if their genotype data are missing or show mismatches.
##########################################################
# add cohort-specific quality control code


##########################################################
# To derive 12 cell counts from beta-values of complete DNA methylation data
##########################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("EpiDISH")
library(EpiDISH) # infer the fractions of a priori known cell subtypes
data(cent12CT.m)
cc <- epidish(beta.m = beta_val, ref.m = cent12CT.m, method = c("RPC")) # beta_val is a data matrix with rows labeling the CpG sites (should use same ID as in ref.m) and columns labeling samples
cell_counts <- as.data.frame(cc$estF)
cell_counts <- data.frame(samp_id = row.names(cell_counts), cell_counts)

# check
cc_na <- table(apply(cell_counts, 1, anyNA))
message(sprintf("Cell counts were successfully predicted for %s individuals.", cc_na[1]))
if (length(cc_na) > 1) {
    message(sprintf("Cell counts contain NAs for for %s individuals. Please check the input data.", cc_na[2]))
}

# merge 12 cell types with other variables
samples <- merge(samples, cell_counts, by = "samp_id")

##########################################################
# To merge principal components with other variables
##########################################################
# add code


##########################################################
# To derive samp_ids from ind_ids
##########################################################
samples <- samples %>%
    mutate(
        samp_id = case_when(
            time_point == 0 ~ paste0(ind_id, "1"),
            time_point == 5 ~ paste0(ind_id, "2"),
            time_point == 9 ~ paste0(ind_id, "3")
        )
    )


##########################################################
# Save
##########################################################
write.table(samples, "data/samp_matrix.txt", row.names = F, quote = F, sep = "\t")



# To curate genotype data later, extract the ind_id and samp_id of participants who possess methylation data for each time point
p1_id <- samples %>%
    filter(time_point == 0) %>%
    select(FID = ind_id, IID = samp_id) %>%
    distinct()
write.table(p1_id, "data/time/gsamples_t1.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")

p2_id <- samples %>%
    filter(time_point == 5) %>%
    select(FID = ind_id, IID = samp_id) %>%
    distinct()
write.table(p2_id, "data/time/gsamples_t2.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")

p3_id <- samples %>%
    filter(time_point == 9) %>%
    select(FID = ind_id, IID = samp_id) %>%
    distinct()
write.table(p3_id, "data/time/gsamples_t3.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
