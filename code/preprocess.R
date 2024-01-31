#!/usr/bin/Rscript

###################################################################
##### Set-up #####
###################################################################

##### Libraries #####
library(argparse)
library(foreach)
library(tidyverse)

##### Arguments #####
parser <- ArgumentParser()
parser$add_argument("--meth", type = "character")
parser$add_argument("--geno", type = "character")
parser$add_argument("--moutput", type = "character")
parser$add_argument("--goutput", type = "character")
Args <- parser$parse_args()


###################################################################
##### Dataset #####
###################################################################

##### Methylation #####
cat(paste0("Reading methylation: ", Args$meth, "\n"))
M <- readRDS(Args$meth)
cat("Dimensions of methylation: ", dim(M), "\n")

##### Genotypes #####
lines <- readLines("data/dataset.raw")

# Process the header
header <- strsplit(lines[1], " ")[[1]]
data_names <- sapply(header[7:length(header)], function(x) strsplit(x, "_")[[1]][1])
geno <- data.frame(matrix(ncol = length(data_names) + 1, nrow = 0))
colnames(geno) <- c("IID", data_names)

# Process and write the data lines
for (i in 2:length(lines)) {
    split_line <- strsplit(lines[i], " ")[[1]]
    ALN <- split_line[2]
    data <- split_line[7:length(split_line)]
    geno <- rbind(geno, c(ALN, data))
}
colnames(geno) <- c("IID", data_names)
rownames(geno) <- geno$IID
geno$IID <- NULL
cat("Dimensions of genotype matrix: ", dim(geno), "\n")

##### IDs #####
id <- colnames(M)
id <- intersect(id, rownames(geno))
cat(paste0("Number of samples in final dataset: ", length(id), "\n"))

###################################################################
##### Process #####
###################################################################

##### Methylation #####
M <- M[, id, drop = FALSE]
cat("Saving methylation data ...", dim(M), "\n")
saveRDS(M, "data/meth_matrix.rds")

##### Genotypes row-sample col-snp -> row-snp col-sample #####
df <- t(geno[id, , drop = FALSE]) # transposes the matrix
colnames(df) <- id
df <- as.data.frame(df)
cat("Saving genotype data: \n")
saveRDS(df, Args$goutput)


q("no")
