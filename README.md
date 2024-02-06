## Background
Prior research has demonstrated certain clinical traits change over time. However, the molecular mechanisms underlying these temporal changes remain inadequately explored. DNA methylation exhibits obvious time trajectories in some CpG sites. Thus, the identification of longitudinal methylation QTLs could furnish new insights into the mechanisms driving these changes. We have used the ALSPAC cohort to estimate age x snp interactions on DNA methylation levels. We restricted our discovery analysis to known additive effects previously reported through GoDMC, and only analysed influences on methylation trajectories by age amongst children. 

## Aim
These scripts aim to re-estimate the discovered age x snp interactions on DNA methylation levels in independent cohorts. The discovered SNP-CpG pairs that show an age interaction are in the `/data` directory. These scripts will take genotype and DNA methylation data + covariates, and re-estimate the age x SNP interactions for the previously discovered interacting mQTLs.

## Inclusion criteria

- Child participants (e.g. below age 20)
- DNA methylation collected at a minimum of two distinct time points
- Genotype data imputed to a recent reference panel (e.g. 1000 genomes, HRC, Topmed etc)

## Overview of datasets required

1. Individual and sample identifiers and variables

    - **ind_id**: A unique identifier assigned to each individual, e.g. 34384_A
    - **time_point**: 0 for Birth; 5 for 5y; 9 for 9y
    - **samp_id**: A unique identifier assigned to each blood sample collected from individuals at various time points, e.g. 34384_A1, 34384_A2, and 34384_A3 denote three distinct samples from the same individual, 34384_A1 corresponding to collection at Birth, 34384_A2 corresponding to collection at 5y, and 34384_A3 corresponding to collection at 9y, respectively.
    - **age**: The ages of the children when their blood samples were collected were distinct from and provided more precise temporal information than the `time_point` variable.
    - **Covariates**: sex, cell counts (12 cell counts from Salas et al. 2022 predicted with EpiDISH), methylation batch, and 10 genetic principal components of ancestry.

2. Methylation data: The analysis is capable of incorporating data from both 450k and EPIC arrays. If the IDs in the methylation data do not match the `samp_id`, it is necessary to rely on a linker file to modify them to the corresponding `samp_id`.

   Note: while all target CpG sites are present in the 450k array dataset, some are not available in the EPIC array dataset. This will not impede the analysis, although it is worth noting that the sample size for the CpG sites exclusive to the 450k array will be comparatively smaller.

3. Genotype data: Best guess genotypes imputed to Haplotype Reference Consortium (HRC) panel (hg19/GRCh37 human reference genome). File format will be bed/fam/bim. SNP ids will be generated in pipeline. If the FID and IID in the .fam file do not match the `ind_id` and `samp_id`, they will be modified to `ind_id` and `samp_id` respectively in the pipeline.



## Input
### 1. sample matrix
Each row corresponds to a blood sample (samp_id)

Column: ind_id, samp_id, time_point, age, sex, batch, cell counts, pc1-10 are mandatory. Time_num is optional.

The analysis-ready sample matrix looks like:
```r
head(samples) # Dummy data
# ind_id samp_id time_point time_num batch age sex CD4Tnv Baso CD4Tmem Bmem Bnv Treg CD8Tmem CD8Tnv Eos NK Neu Mono pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10
# 30045_A 30045_A1 0 3 BCD069 0 0 0.05 0.03 0.00 0.032 0 0.14 0.058 0.024 0.020 0.000 0.6 0.03 0.0009 0.004 -0.0077 0.0134 0.004 0.004 0.003 0.001 0.001 -0.004
# 30045_A 30045_A3 9 3 BCD078 10 0 0.08 0.08 0.065 0.076 0 0.060 0.041 0.031 0.080 0.008 0.4 0.02 0.0009 0.004 -0.0077 0.0134 0.004 0.004 0.003 0.001 0.001 -0.004
```
### 2. meth matrix
Row: 662 CpGs

Column: samp_id

The analysis-ready methylation matrix looks like:
```r
M[1:3,1:5]  # Dummy data
#           30045_A1 30045_A2 30045_A3 30046_A1 30046_A2
# cg00025357 0.2667352  0.2616095  0.2825195 0.3413778  0.3328760
# cg00695362 0.4351751  0.3275295  0.3528334 0.1353732  0.2413153
# cg01808739 0.2421042  0.2196921  0.2230272 0.7255114  0.6212163
```
### 3. geno matrix
Row: 632 target SNPs

Column: samp_id

The analysis-ready genotype matrix looks like:
```r
df[1:3,1:5] # Dummy data
#           30045_A1 30045_A2 30045_A3 30046_A1 30046_A2
# 1:2038893       1       1       1       2       2
# 1:3591143       2       2       2       1       1
# 1:10700448      2       2       2       2       2
```
### 4. analysis matrix
The target set of 673 CpG-SNP pairs: `data/analysis_matrix.txt` file.

```
snp        cpg        rsid
1:17564465 cg00025357  rs12030076
1:40481100 cg00695362  rs11207310
```

   Note: Generation R should have data for all target SNP and CpG sites. If any SNP or CpG site has no data, it suggests that the pipeline may be incompatible with the original data structure of Generation R. Please adjust the code or contact the GitHub repository owner.


## Output
The files to be shared are `results/summariseMeth.tsv`, `results/lmm.meqtl` and `results/lmm.Rout`. 
Besides, it would be advantageous to include the INFO scores for the 632 target SNPs. This could be achieved by supplying a separate file or by incorporating these scores into the existing `results/lmm.meqtl` file.


## Code and scripts
### 1. sample matrix
To prepare the analysis-ready sample matrix, use the `code/prepare_samp_matrix.R` script after modifying it according to the specific data traits of the cohort.

### 2. meth matrix
To prepare the analysis-readymethylation matrix, run:
```r
target_CpG <- read_lines("data/cpg_hits.txt") # the provided file documenting the names of 662 target CpGs
M <- M[rownames(M) %in% target_CpG, match(samples$samp_id,colnames(M))]  # trim the methylation data frame based on the selected CpGs (rows) and the samples (columns) to be included
table(samples$samp_id==colnames(M)) # check

saveRDS(M, "data/meth_matrix.rds")
```

### 3. geno matrix
Assume that the complete genotype data for Generation R is contained within the files named `gdata.bed`, `gdata.fam`, and `gdata.bim`.

If the individual identifiers (FID and IID) in the genotype data differ from the `ind_id` in the sample matrix, a linking file is necessary to associate them.
```r
genos <- read.table("data/gdata.fam", header=F, sep=" ")
# add code to load in the linker file which contains FID and IID columns and corresponding ind_id
genos$V1 <- linker$ind_id[match(genos$V1, linker$genos_id)] # FID column
genos$V2 <- genos$V1 # the IID column will be matched to samp_id in the next step
write.table(genos, file = "data/gdata_s.fam", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

To extract genotype data of participants at each time point and combine them together, run:
```
./code/time_specific_geno.sh

Rscript --no-save --no-restore --verbose code/preprocess.R \
  --meth data/meth_matrix.rds \
  --geno data/dataset.raw \
  --moutput data/meth_matrix.rds \
  --goutput data/geno_matrix.rds
```

### 4. data analysis
The Linear Mixed Model (LMM) serves as the principal analytical method. For detailed implementation, refer to the `code/lmm.R` file.

To install all the required packages, open R in the working directory and run:
```r
install.packages("renv")
renv::restore()
```

The shell command is expected to be completed in under an hour on high-performance computing systems.
```
Rscript --no-save --no-restore --verbose code/lmm.R \
    --samp=data/samp_matrix.txt \
    --geno=data/geno_matrix.rds \
    --meth=data/meth_matrix.rds \
    --anal=data/analysis_matrix.txt \
    --sout=results/summariseMeth.tsv \
    --out=results/lmm.meqtl > results/lmm.Rout 2>&1
```

