##########################################################
# Code information
# Title: Code for performing xTWAS by FUSION
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################


## Process gwas summary statistics
python munge_sumstats.py \
    --sumstats ${raw_gwas} \
    --snp ${MarkerName} \
    --a1 ${Allele1} \
    --a2 ${Allele2} \
    --frq ${Freq} \
    --p ${Pvalue} \
    --N ${sample_number} \
    --out ${gwas_sumstat}


## Calculate weights with own data
Rscript FUSION.compute_weights.R \
    --bfile ${OUT} \
    --tmp ${OUT}.tmp \
    --out ${FINAL_OUT} \
    --PATH_gemma ${gemma} \
    --PATH_plink ${plink} \
    --PATH_gcta ${gcta} \
    --covar ${covar_path} \
    --verbose 2 \
    --save_hsq \
    --models top1,blup,lasso,enet


## Summarize the weights
Rscript FUSION.profile_wgt.R ${wgt_list} > ${wgt_profile}


## Calculate phenotype-GC association
Rscript FUSION.assoc_test.R \
    --sumstats ${gwas_sumstat} \
    --weights ${wgt_path} \
    --weights_dir ${wgt_dir} \
    --ref_ld_chr ${fusion_ldref}. \
    --chr ${CHR} \
    --out ${asso_out}


## Perform conditional tests
Rscript FUSION.post_process.R \
    --sumstats ${gwas_sumstat} \
    --input ${top_dir} \
    --out GC.${CHR}.top.analysis \
    --ref_ld_chr ${fusion_ldref}. \
    --chr ${CHR} \
    --plot