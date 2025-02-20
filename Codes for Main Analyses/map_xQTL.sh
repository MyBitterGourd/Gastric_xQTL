##########################################################
# Code information
# Title: Code for mapping xQTLs
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################


## Prepare genotype file
plink --vcf ${raw_geno} \
    --geno 0.05 \
    --mind 0.05 \
    --maf 0.01 \
    --hwe 0.000001 \
    --make-bed \
    --double-id \
    --out ${geno_file}


## cis-QTL mapping: all variant-phenotype pairs
python3 -m tensorqtl ${geno_file} \
    ${pheno_file} \
    ${nominal_prefix} \
    --covariates ${covar_file} \
    --mode cis_nominal


## cis-QTL mapping: permutations
python3 -m tensorqtl ${geno_file} \
    ${pheno_file} \
    ${permutation_prefix} \
    --covariates ${covar_file} \
    --mode cis


## cis-QTL mapping: conditionally independent QTLs
python3 -m tensorqtl ${geno_file} \
    ${pheno_file} \
    ${independent_prefix} \
    --covariates ${covar_file} \
    --cis_output ${permutation_prefix}.txt.gz \
    --mode cis_independent