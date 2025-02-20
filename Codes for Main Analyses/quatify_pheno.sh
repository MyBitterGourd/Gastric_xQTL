##########################################################
# Code information
# Title: Code for qualifying phenotypes
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################


## ------------------------------ MARK: Prepare expression (follow GTEx)

## List the chromosomes in the VCF
tabix --list-chroms ${vcf} > ${vcf_chr_list}


## Normalize expression between samples using TMM
eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
    ${sample_participant_lookup} ${vcf_chr_list} ${prefix} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm


## Calculate PEER factors
Rscript run_PEER.R ${prefix}.expression.bed.gz ${prefix} ${num_peer}


## Combine covariates
combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} \
    --genotype_pcs ${genotype_pcs} \
    --add_covariates ${add_covariates}



## ------------------------------ MARK: Prepare PSI (using LeafCutter)

## Quantify splicing
for i in {1..262}; do
    echo ${i}
    id=$(sed -n "${i}p" RNAseq.sample.list)
    sh leafcutter/scripts/bam2junc.sh ${id}.sort.bam ${id}.junc
    echo ${id}.junc >> juncfiles.txt
done


## Cluster introns
python leafcutter/clustering/leafcutter_cluster.py \
    -j juncfiles.list \
    --minclureads=30 \
    --mincluratio=0.001 \
    --maxintronlen=500000 \
    -o 262sample.Stomach


## Map clusters to genes
Rscript leafcutter/scripts/map_intron_to_genes_check.r \
    --count.file 262sample.Stomach_perind.counts.gz \
    --exon.file gencode.v29lift37.exon.table \
    --output_dir 262sample.Stomach.IntronToGencode.v29lift37.genes.txt


## Quality control of introns
R
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())

norm <- fread("262sample.Stomach_perind.counts.gz", h=T, data.table=F)
norm_count <- fread("262sample.Stomach_perind_numers.counts.gz", h=F, data.table=F)

# introns with few counts
# introns without any read counts in >50% of samples were filtered out
rownames(norm) <- norm$chrom
ratio <- apply(norm[,2:ncol(norm)],2,function(x){ a <- as.numeric(sapply(strsplit(x,'/'),'[',1)); b <- as.numeric(sapply(strsplit(x,'/'),'[',2)); return(a/b)})
rownames(ratio) <- rownames(norm)
pct_zero <- apply(ratio,1,function(x){length(which(x > 0))/262})
pct_zero_pass <- pct_zero[pct_zero >= 0.5]

rownames(norm_count) <- norm_count$V1
norm_count <- norm_count[,-1]
qc_1.1 <- apply(norm_count,1,function(x){length(which(x>0))})
qc_1.1pass <- qc_1.1[qc_1.1 >= 262*0.5]


# introns with low complexity
# introns with fewer than max(0.1n) unique values, where n is the sample size, were filtered out
n_unique <- apply(ratio,1,function(x){x[is.na(x)]<-0;return(length(unique(x)))})
n_unique_pass <- n_unique[n_unique >= 0.1*262] 

# introns with low complexity
# introns with insufficient variability across samples were removed (z of cluster read fractions across individuals)
zscore_df <- apply(ratio,1,function(x){x[is.na(x)]<-0; (x-mean(x))/sd(x)})
zscore_outlier <- apply(zscore_df,2,function(x){
    (sum(abs(x)<0.25)>=(262-3) ) & ( sum(abs(x)>6)<=3)
})
zscore_pass <- zscore_outlier[zscore_outlier == FALSE] 

# passed introns
int_gen <- fread("262sample.Stomach.IntronToGencode.v29lift37.genes.txt", h=T)
int_gen$clu_id <- sapply(strsplit(int_gen$clu, ":"), "[", 2)

norm <- fread("262sample.Stomach_perind.counts.gz", h=T, data.table=F)
norm_qc <- subset(norm, chrom %in% names(pct_zero_pass))
norm_qc <- subset(norm_qc, chrom %in% names(n_unique_pass))
norm_qc <- subset(norm_qc, chrom %in% names(zscore_pass))
norm_qc$clu_id <- sapply(strsplit(norm_qc$chrom, ":"), "[", 4)
norm_qc_final <- subset(norm_qc, clu_id %in% int_gen$clu_id)
norm_qc_final$chrom <- gsub("chr", "", norm_qc_final$chrom)
norm_qc_final$chr <- sapply(strsplit(norm_qc_final$chrom, ":"), "[", 1)
d <- which(norm_qc_final$chr %in% c(1:22))
fwrite(norm_qc_final[d,1:263], "intron/262sample.Stomach_perind.filtered.counts.txt", col=T, row=F, quo=F)

# normalization
# The filtered counts were normalized using the prepare_phenotype_table.py script from LeafCutter
# Then the resulting per-chromosome files merged and converted to BED format with the start/end position corresponding to the TSS of the gene
gzip 262sample.Stomach_perind.filtered.counts.txt

python leafcutter/scripts/prepare_phenotype_table.py 262sample.Stomach_perind.filtered.counts.txt.gz -p 10


# Concatenate BED files
R
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())

dat <- fread("262sample.Stomach_perind.filtered.counts.txt.gz.qqnorm_1")
for (i in 2:22) {
	temp <- fread(paste0("262sample.Stomach_perind.filtered.counts.txt.gz.qqnorm_", i))
	print(nrow(temp))
	print(all(colnames(dat)==colnames(temp)))
	dat <- rbind(dat, temp)
}

colnames(dat)[1] <- "Chr"
dat1 <- dat[order(dat$Chr, dat$start, dat$end),]
colnames(dat1)[1] <- "#Chr"
fwrite(dat1, "262sample.Stomach_perind.filtered.counts.txt.gz.qqnorm_ChrAll", col=T, row=F, quo=F, sep="\t")
system("gzip 262sample.Stomach_perind.filtered.counts.txt.gz.qqnorm_ChrAll")

dat1 <- dat[order(dat$Chr, dat$start, dat$end),]
dat1$Chr <- paste0("chr", dat1$Chr)
fwrite(dat1[,1:4], "262sample.Stomach.Intron.filtered.bed", col=T, row=F, quote=F, sep="\t") 


## Calculate PEER factors
Rscript scripts/run_PEER_check.r \
    --expr.file 262sample.Stomach_perind.filtered.counts.txt.gz.qqnorm_ChrAll.gz \
    --prefix 262sample.Stomach.Intron \
    --n 15


## Combine covariates
combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} \
    --genotype_pcs ${genotype_pcs} \
    --add_covariates ${add_covariates}



## ------------------------------ MARK: Prepare PDUI (using DaPars2)

## Prepare input data for Dapars2
bash bam2wig.sh 262sample_bam.list 8
bash extract_sequencing_depth.sh 262sample_bam.list rawdata/wig 8 wigFile_and_readDepth.txt

python DaPars_Extract_Anno.py \
    -b gencode.v29lift37.annotation.bed \
    -s hg19_gencode_transcript2geneSymbol.txt \
    -o gencode_3utr_annotation.bed
cat gencode_3utr_annotation.bed | cut -f 1 | \
    sort | uniq | grep -v 'M' | grep -v 'X' | \
    grep -v 'Y' | grep -v 'GL' > chrList.txt 
python generate_configure_for_Dapars2.py  \
    --annotation_3utr gencode_3utr_annotation.bed \
    --wigFile_depth wigFile_and_readDepth.txt \
    --coverage_threshold 15 \
    --threads 8


## Run Dapars2 analysis on each chromosome
for i in {1..22}; do
    echo $i
    python Dapars2_Multi_Sample.py Dapars2_running_configure.txt chr$i
done


## Merge Dapars2 results on individual chromosomes into one
Rscript merge_dapars2_res_by_chr.R Dapars2_out 262sample_bam.list chrList.txt


## Prepare 3'UTR location file
python extract_3UTR_location.py \
    --dapars_res Dapars2_res.all_chromosomes.txt \
    --output 3UTR_location.txt


## Calculate PEER factors
Rscript 3UTR_impute_peer.R known.cov.txt 5


## Combine covariates
combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} \
    --genotype_pcs ${genotype_pcs} \
    --add_covariates ${add_covariates}