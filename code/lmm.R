#!/usr/bin/Rscript

###################################################################
##### Set-up #####
###################################################################
##### Libraries #####
suppressMessages(library(argparse))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(lmerTest))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(ewaff))
library(tidyverse)
options(stringsAsFactors = F)


##### Arguments #####
# registerDoMC(detectCores())
parser <- ArgumentParser()
parser$add_argument("--samp", type = "character")
parser$add_argument("--geno", type = "character")
parser$add_argument("--meth", type = "character")
parser$add_argument("--anal", type = "character")
parser$add_argument("--sout", type = "character")
parser$add_argument("--out", type = "character")
Args <- parser$parse_args()


###################################################################
##### load in datasets #####
###################################################################
##### Genotypes #####
cat(paste0("Reading genotypes: ", Args$geno, "\n"))
geno <- readRDS(Args$geno)

##### Methylation #####
cat(paste0("Reading methylation: ", Args$meth, "\n"))
meth0 <- readRDS(Args$meth)

##### Covariates #####
cat(paste0("Reading samples: ", Args$samp, "\n"))
samples <- read_tsv(file = Args$samp, col_names = TRUE)

##### Samples #####
covar <- samples[match(colnames(geno), samples$samp_id), ]
variable_names <- c("ind_id", "samp_id", "time_point", "batch", "age", "sex", "CD4Tnv", "Baso", "CD4Tmem", "Bmem", "Bnv", "Treg", "CD8Tmem", "CD8Tnv", "Eos", "NK", "Neu", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10")

cat(paste0("Methylation and Genotype IDs equal: ", all(colnames(geno) == colnames(meth0)), "\n"))
cat(paste0("Sample and Genotype IDs equal: ", all(colnames(geno) == covar$samp_id), "\n"))
cat(paste0("The names of the variables are the same: ", all(variables_names %in% colnames(covar)), "\n"))

#### Outlier exclusions for each time point #####
message("Generating summary stats of methylation")
message("Eliminating methylation outliers: β-values 3 SDs away from the time-point specific mean")
summariseMeth <- function(X)
{
  norm.beta <- X
	norm.beta.copy<-ewaff.handle.outliers(norm.beta, method=c("iqr"), iqr.limit=3, winsorize.pct=0.05)
  message("Counting outliers in methylation matrix")
	outliers <-norm.beta.copy[[2]]

	message("Estimating means, SDs, medians")
	means <- rowMeans(norm.beta.copy[[1]], na.rm=T)
	sds <- rowSds(norm.beta.copy[[1]], na.rm=T)
	medians <- rowMedians(norm.beta.copy[[1]], na.rm=T)
  
	dat <- data.frame(cpg=rownames(norm.beta.copy[[1]]), mean=means, median=medians, sd=sds, outlier=rowSums(outliers[, 1:2], na.rm = TRUE))
	return(list(norm.beta = norm.beta.copy[[1]], dat = dat))
}

meth <- meth0
stats <- tibble(cpg=rownames(meth))
time_lev <- as.character(levels(as.factor(covar$time_point)))
for (lev in time_lev) {
ret <- summariseMeth(meth0[, covar[covar$time_point == lev, ]$samp_id])
meth[, covar[covar$time_point == lev, ]$samp_id] <- ret$norm.beta
stats <- stats %>% mutate(
  !!paste0("mean_", lev) := ret$dat$mean[match(cpg,ret$dat$cpg)], 
  !!paste0("sd_", lev) := ret$dat$sd[match(cpg,ret$dat$cpg)], 
  !!paste0("median_", lev) := ret$dat$median[match(cpg,ret$dat$cpg)], 
  !!paste0("outlier_", lev) := ret$dat$outlier[match(cpg,ret$dat$cpg)])
  cat(paste0("Outliers for time point ", lev, ": ", sum(ret$dat$outlier, na.rm = TRUE), "\n"))
}

stats <- stats %>% mutate(
  mean = rowMeans(meth, na.rm = T),
  sd = rowSds(meth, na.rm = T),
  median = rowMedians(meth, na.rm = T),
  outlier = rowSums(stats %>% select(starts_with("outlier_")), na.rm = TRUE)
)

write.table(stats, Args$sout, row.names = FALSE, quote = FALSE, sep = "\t")
dif_pos <- which(!is.na(meth) != !is.na(meth0), arr.ind = TRUE)
cat(paste0("Number of outliers: ", nrow(dif_pos), "\n"))


###################################################################
##### Analysis #####
##################################################################
##### Matrix sizes #####
cat(paste0("Methylation size: ", nrow(meth), " x ", ncol(meth), "\n"))
cat(paste0("Genotype size: ", nrow(geno), " x ", ncol(geno), "\n"))
cat(paste0("Covariates size: ", nrow(covar), " x ", ncol(covar), "\n"))

##### Outputs to be captured #####
cat(paste0("Reading analysis: ", Args$anal, "\n"))
analysis <- fread(Args$anal, header = TRUE, sep = "\t", stringsAsFactors = F, data.table = F, showProgress = F)
analysis$beta_lmm <- NA
analysis$se_lmm <- NA
analysis$p_lmm <- NA
analysis$beta_lm <- NA
analysis$se_lm <- NA
analysis$p_lm <- NA
analysis$singular <- NA
analysis$n <- NA
analysis$mm <- NA
analysis <- analysis[(analysis$snp %in% rownames(geno)) & (analysis$cpg %in% rownames(meth)), ]
cat(paste0("Number of CpG-SNP pairs: ", nrow(analysis), "\n"))

##### Main body of analysis #####
cat(paste0("Analysis:\n"))

for (i in seq_len(nrow(analysis))) {
  ### Curate analysis dataset
  snp <- geno[rownames(geno) == analysis$snp[i], ]
  cpg <- meth[rownames(meth) == analysis$cpg[i], ]
  dataset <- data.frame(
    cpg = cpg, snp = snp, age = covar$age, sex = covar$sex, batch = covar$batch,
    CD4Tnv = covar$CD4Tnv, Baso = covar$Baso, CD4Tmem = covar$CD4Tmem, Bmem = covar$Bmem, Bnv = covar$Bnv, Treg = covar$Treg, CD8Tmem = covar$CD8Tmem, CD8Tnv = covar$CD8Tnv, Eos = covar$Eos, NK = covar$NK, Neu = covar$Neu,
    pc1 = covar$pc1, pc2 = covar$pc2, pc3 = covar$pc3, pc4 = covar$pc4, pc5 = covar$pc5, pc6 = covar$pc6, pc7 = covar$pc7, pc8 = covar$pc8, pc9 = covar$pc9, pc10 = covar$pc10, ind_id = covar$ind_id
  )
  dataset$age_snp <- dataset$age * dataset$snp

  dataset <- na.omit(dataset)
  analysis$n[i] <- nrow(dataset)
  analysis$mm[i] <- round(median(dataset$cpg, na.rm = TRUE), 4) # median methylation

  ### The formulas implemented
  lm_filter <- formula(
    "cpg ~ snp + age + as.factor(sex) + as.factor(batch) + CD4Tnv + Baso + CD4Tmem + Bmem + Bnv + Treg + CD8Tmem + CD8Tnv + Eos + NK + Neu+ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
  )
  lmm_formula <- formula(
    "lm(cpg ~ as.factor(batch) + CD4Tnv + Baso + CD4Tmem + Bmem + Bnv + Treg + CD8Tmem + CD8Tnv + Eos + NK + Neu, data = dataset)$resid ~
        age + snp + age_snp + as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (age | ind_id)"
  )

  if (length(unique(dataset$sex)) == 1) {
    lm_filter <- update(lm_filter, . ~ . - as.factor(sex))
    lmm_formula <- update(lmm_formula, . ~ . - as.factor(sex))
  }

  ### mQTLs
  model <- lm(lm_filter, data = dataset)
  coef_lm <- rownames(summary(model)$coefficients)
  if ("snp" %in% coef_lm) {
    analysis$beta_lm[i] <- round(summary(model)$coefficients["snp", 1], 6)
    analysis$se_lm[i] <- round(summary(model)$coefficients["snp", 2], 6)
    analysis$p_lm[i] <- summary(model)$coefficients["snp", 4]
  } else {
    cat(paste("line", i, ": no coefficients for snp\n"))
  }

  ### longitudinal mQTLs
  model <- lmer(lmm_formula, data = dataset, control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  coef_lmm <- dimnames(summary(model)$coefficients)[[1]]
  if ("age_snp" %in% coef_lmm) {
    analysis$beta_lmm[i] <- round(summary(model)$coefficients["age_snp", 1], 6)
    analysis$se_lmm[i] <- round(summary(model)$coefficients["age_snp", 2], 6)
    analysis$p_lmm[i] <- summary(model)$coefficients["age_snp", 5]
    analysis$singular[i] <- "n"
    if (!is.null(summary(model)$optinfo$conv$lme4$messages)) {
      analysis$singular[i] <- "y"
    }
  } else {
    cat(paste("line", i, ": no coefficients for age_snp\n"))
  }

  ### Progress indicator
  if (i %% 10 == 0) {
    cat(paste0(" ", i, "--DONE\n"))
  }
}



###################################################################
##### Save results #####
###################################################################

# effect allele, other allele, and minor allele frequency
frq <- read_table("data/dataset.frq", col_names = TRUE) # "data/dataset.frq" was generated from "code/time_specific_geno.sh"
frq$snp_p <- sub("^(.*[0-9])_(.*)", "\\1", frq$SNP)
analysis$EFFECT_ALL <- frq$A1[match(analysis$snp, frq$snp_p)]
analysis$OTHER_ALL <- frq$A2[match(analysis$snp, frq$snp_p)]
analysis$maf_G <- frq$MAF[match(analysis$snp, frq$snp_p)]

cat(paste0("Saving dataset: ", Args$out))
analysis <- analysis[, c("snp", "cpg", "rsid", "beta_lmm", "se_lmm", "p_lmm", "n", "mm", "singular", "beta_lm", "se_lm", "p_lm", "EFFECT_ALL", "OTHER_ALL", "maf_G")]
write.table(analysis, Args$out, row.names = FALSE, quote = FALSE, sep = "\t")


cat(paste0("\n\nMeQTL Done"))
q("no")
