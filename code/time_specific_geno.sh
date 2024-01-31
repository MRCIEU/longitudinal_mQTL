#!/bin/bash

module add apps/plink/1.90

# Make time point specific genotypes
for i in {1..3}
do
  cp data/gdata.bed data/time/${i}_hrc.bed
  cp data/gdata.bim data/time/${i}_hrc.bim
  cp data/gdata_s.fam data/time/${i}_hrc.fam
  # change IID from ind_id to samp_id
  awk -v i="${i}" '{$2 = $2i} 1' OFS=' ' data/time/${i}_hrc.fam > temp && mv temp data/time/${i}_hrc.fam
done

# Trim the genotype data based on the samples in each time point (gsamples_t${i}.txt)
for i in {1..3}
do
    plink --bfile "data/time/${i}_hrc" \
        --set-all-var-ids @:#_$1_$2 \
        --keep "data/time/gsamples_t${i}.txt" \
        --make-bed \
        --out "data/time/${i}_hrc_1"       
done


# Merge plink files
for i in {2..3}
do 
    echo "data/time/${i}_hrc_1"
done > data/time/mergefile.txt

plink  --bfile data/time/1_hrc_1 --merge-list data/time/mergefile.txt --make-bed --out data/hrc
wc -l data/hrc.bim
wc -l data/hrc.fam

# Trim the genotype data based on the 632 target SNPs (snp_hits.txt)
plink --bfile data/hrc --recodeA --freq --extract "data/snp_hits.txt" --out data/dataset


# Remove temporary files
rm data/time/1_* data/time/2_* data/time/3_*


