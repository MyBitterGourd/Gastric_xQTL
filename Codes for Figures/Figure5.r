##########################################################
# Code information
# Title: Code of Figure 5
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################


R
options(stringsAsFactors=F)
library(data.table)
library(tidyverse)
library(janitor)
library(ggplot2)
library(ggrepel)
library(patchwork)
rm(list=ls())


## ------------------------------ Fig.5A top

## Load spTWAS results
twas <- fread("input/sptwas.chr.sig.Fusion.FDR10.list") %>%
    filter(gene_name == "FDPS") %>%
    filter(ID == "1:155289479:155290200:clu_10554:ENSG00000160752.14_2")


## Set the region around the spTWAS signal
upper <- twas$Start
lower <- twas$End
min.limit <- mean(c(lower,upper)) - 500000
max.limit <- mean(c(lower,upper)) + 500000


## Load results of conditional analysis
dat <- read.csv("input/sptwas.Fusion.FDR10.chr1_GWAS.condition.intron.csv") %>%
    drop_na(pv) %>%
    mutate(pos = as.numeric(tstrsplit(snp, ":")[[2]])) %>%
    filter(pos > mean(c(lower,upper)) - 500000) %>%
    filter(pos < mean(c(lower,upper)) + 500000)


## Separate the results and then combine them
rawgwas <- dat %>%
    select(pos, pv) %>%
    mutate(note = "Main GWAS")
condgwas <- dat %>%
    select(pos, pv = cond.pv) %>%
    mutate(note = "Conditioned on\npredictor")
dat1 <- rbind(rawgwas,condgwas)


## Generate the top plot of GWAS
GWASp0 <- ggplot(data = dat1, aes(x = pos/1e6, y = -log10(pv))) +
    geom_point(
        aes(color = note), pch = 16,
        size = 1.5, alpha = 1
    ) +
    scale_color_manual(
        values = rev(c("grey", "#F08E59")),
        name = ""
    ) +
    xlab("Chr1 physical position (Mb)") +
    ylab(expression("GWAS"~"-log"[10]*"("*italic("P")*")")) +
    scale_y_continuous(expand = c(0.1, 0.1))
GWASp1 <- GWASp0 + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = "white", linetype = "solid", linewidth = 0.5, colour = "black")
)



## ------------------------------ Fig.5A middle

## Load SNPs in the region and calculate MAF
dat <- fread("input/sQTL.TWAS.allpairs.txt") %>%
    filter(phenotype_id == "1:155289479:155290200:clu_10554:ENSG00000160752.14_2")
fwrite(dat[, "variant_id"], "temp/chr1.snp.list", col=F, row=F, quo=F, sep="\t")

system("plink --bfile input/geno --extract temp/chr1.snp.list --make-bed --out temp/chr1")
system("plink --bfile temp/chr1 --freq --out temp/chr1")


## Load the original sQTL results
freq <- fread("temp/chr1.frq")
all <- dat %>%
    inner_join(freq, by = c("variant_id" = "SNP")) %>%
    select(variant_id, A1, A2, MAF, slope, slope_se, pval_nominal) %>%
    rename(SNP = variant_id, freq = MAF, b = slope, se = slope_se, p = pval_nominal) %>%
    mutate(N = 262)
fwrite(all, "temp/chr1.ma", col=T, row=F, quo=F, sep="\t")
write.table("1:155285854:T:C", "temp/condition.snplist.chr1", row=F, col=F, quo=F, sep="\t")


## Perform conditional regression
system("gcta64 --bfile temp/chr1 --cojo-file temp/chr1.ma --cojo-cond temp/condition.snplist.chr1 --out temp/condition.chr1.TopsQTL")
system("rm temp/chr1.*")


## Prepare output for mPLOT!
con <- fread("temp/condition.chr1.TopsQTL.cma.cojo") %>%
    select(SNP, bC, bC_se, pC)
final <- all %>%
    left_join(con, by = "SNP")
fwrite(final, "temp/condition.chr1.TopsQTL.txt", col=T, row=F, quo=F, sep="\t")


## Load the results of spTWAS
twas <- fread("input/sptwas.chr.sig.Fusion.FDR10.list") %>%
    filter(ID == "1:155289479:155290200:clu_10554:ENSG00000160752.14_2")


## Set the region around the spTWAS signal
upper <- twas$Start
lower <- twas$End
min.limit <- mean(c(lower,upper)) - 500000
max.limit <- mean(c(lower,upper)) + 500000


## Load the results of conditional analysis
dat <- final %>%
    mutate(pos = as.numeric(tstrsplit(SNP, ":")[[2]])) %>%
    filter(pos>min.limit & pos<max.limit)


## Separate the results and then combine them
rawgwas <- dat %>%
    select(pos, p) %>%
    mutate(note = "Main sQTL")
condgwas <- dat %>%
    select(pos, p = pC) %>%
    mutate(note = "Conditioned on\ntop sQTL")
dat1 <- rbind(rawgwas,condgwas)


## top GWAS-SNP
sda <- dat %>%
    filter(SNP == "1:155285854:T:C") %>%
    mutate(rsid = "rs10908462")


## Generate the middle plot of sQTL
sQTLp0 <- ggplot(data = dat1, aes(x = pos/1e6, y = -log10(p))) +
    geom_point(
        aes(color = note), pch = 16,
        size = 1.5, alpha = 1 
    ) +
    scale_color_manual(
        values = rev(c("grey","#F08E59")),
        name = ""
    ) +
    ylab(bquote(~atop(italic(FDPS)~"sQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    geom_text(
        data = sda,
        aes(x = pos/1e6, y = -log10(p)+0.7, label = rsid)
    )
sQTLp1 <- sQTLp0 + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    legend.position = c(0.15, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = "white", linetype = "solid", linewidth = 0.5, colour = "black")
)


## ------------------------------ Fig.5A bottom

## Load the original eQTL results
dat <- arrow::read_parquet("input/262sample_eqtl.cis_qtl_pairs.1.parquet") %>%
    filter(grepl("ENSG00000160752", phenotype_id)) %>%
    select(variant_id, pval_nominal) %>%
    mutate(pos = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    select(pos, p = pval_nominal) %>%
    mutate(note = "Main eQTL")


## Load the results of spTWAS
twas <- fread("input/sptwas.chr.sig.Fusion.FDR10.list") %>%
    filter(ID == "1:155289479:155290200:clu_10554:ENSG00000160752.14_2")


## Set the region around the spTWAS signal
upper <- twas$Start
lower <- twas$End
min.limit <- mean(c(lower,upper)) - 500000
max.limit <- mean(c(lower,upper)) + 500000
dat <- dat %>%
    filter(pos>min.limit & pos<max.limit)


## Generate the bottom plot of eQTL
eQTLp0 <- ggplot(data = dat, aes(x = pos/1e6, y = -log10(p))) +
    geom_point(
        aes(color = note), pch = 16,
        size = 1.5, alpha = 1
    ) +
    scale_color_manual(values = c("black"), name = "") +
    xlab("Chr1 physical position (Mb)") +
    ylab(bquote(~atop(italic(FDPS)~"eQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(
        breaks = seq(0, 25, 5),
        labels = seq(0, 25, 5),
        limits = c(0, 26),
        expand=c(0.1, 0.1)
    )
eQTLp1 <- eQTLp0 + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.line = element_blank(),
    legend.position = c(0.15, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = "white", linetype = "solid", linewidth = 0.5, colour = "black")
)



## ------------------------------ Function for eQTL boxplot

## Set alleles
REF <- tstrsplit("1:155289545:T:G", ":")[[4]]
ALT <- tstrsplit("1:155289545:T:G", ":")[[3]]
WT <- paste0(REF, REF)
HET <- paste0(REF, ALT)
HOM <- paste0(ALT, ALT)


## Load expression data
exp <- fread("input/262sample_Stomach.expression.bed.gz") %>%
    filter(gene_id == "ENSG00000160752.14") %>%
    select(-c(1:3)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    rename(exp = 2)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr1.fst") %>%
    tidyfst::filter_fst(ID == "1:155289545:T:G") %>%
    select(-c(1,2,4:9)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    rename(snp = 2) %>%
    mutate(
        snp = case_when(
            snp == "0/0" ~ WT,
            snp %in% c("0/1", "1/0") ~ HET,
            snp == "1/1" ~ HOM,
            .default = NA
        )
    )


## Load covariate data
covar <- fread("input/262sample_Stomach.combined_covariates.txt") %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID")


## Merge data and calculate residuals
final <- exp %>%
    inner_join(covar, by = "ID")
geno <- geno[match(final$ID, geno$ID),]
final <- final %>%
    select(-1) %>%
    mutate_at(c(1:51,53), as.numeric)
res <- summary(lm(exp~., data=final))
geno$residuals <- res$residuals


## Set genotype levels
geno$snp <- factor(geno$snp, level = c(WT, HET, HOM))
WTnum = paste0(WT, "\n(n=", table(geno$snp)[[1]], ")")
HETnum = paste0(HET, "\n(n=", table(geno$snp)[[2]], ")")
HOMnum = paste0(HOM, "\n(n=", table(geno$snp)[[3]], ")")
geno$snp <- ifelse(
    geno$snp == WT,
    WTnum,
    ifelse(
        geno$snp == HET,
        HETnum,
        HOMnum
    )
)
geno$snp <- factor(geno$snp, level = c(WTnum, HETnum, HOMnum))


## Load QTL results
temp <- tidyfst::parse_fst("input/262sample_eqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::filter_fst(phenotype_id == "ENSG00000160752.14") %>%
    tidyfst::filter_fst(variant_id == "1:155289545:T:G")
Slope <- round(temp$slope, 3)
Pvalue <- format(temp$pval_nominal, digits = 3)


## Generate the boxplot of eQTL
eQTLboxp0 <- ggplot(data = geno, aes(x = snp, y = residuals)) +
    geom_boxplot(
        aes(color = snp), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = snp), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#7098BD", "#3F76A7", "#145390"), name = "") +
    xlab("rs11264361 genotype") +
    ylab("Normalized expression") +
    ggtitle(
        bquote(italic(FDPS)),
        bquote(italic(beta)*" = "*.(Slope)~", "~italic("P")*" = "*.(Pvalue))
    )
eQTLboxp1 <- eQTLboxp0 +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none"
    )



## ------------------------------ Function for sQTL boxplot

## Set intron ID
Chr <- tstrsplit("1:155289479:155290200:clu_10554:ENSG00000160752.14_2", ":")[[1]]
upper <- tstrsplit("1:155289479:155290200:clu_10554:ENSG00000160752.14_2", ":")[[2]]
lower <- tstrsplit("1:155289479:155290200:clu_10554:ENSG00000160752.14_2", ":")[[3]]
IntronID <- paste0("chr", Chr, ":", upper, "-", lower)


## Set alleles
REF <- tstrsplit("1:155289545:T:G", ":")[[4]]
ALT <- tstrsplit("1:155289545:T:G", ":")[[3]]
WT <- paste0(REF, REF)
HET <- paste0(REF, ALT)
HOM <- paste0(ALT, ALT)


## Load splicing data
splicing <- fread("input/262sample_Stomach.splicing.bed.gz") %>%
    filter(ID == "1:155289479:155290200:clu_10554:ENSG00000160752.14_2") %>%
    select(-c(1:3)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    rename(splicing = 2)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr1.fst") %>%
    tidyfst::filter_fst(ID == "1:155289545:T:G") %>%
    select(-c(1,2,4:9)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    rename(snp = 2) %>%
    mutate(
        snp = case_when(
            snp == "0/0" ~ WT,
            snp %in% c("0/1", "1/0") ~ HET,
            snp == "1/1" ~ HOM,
            .default = NA
        )
    )


## Load covariate data
covar <- fread("input/262sample_Stomach.combined_covariates.txt") %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID")


## Merge data and calculate residuals
final <- splicing %>%
    inner_join(covar, by = "ID")
geno <- geno[match(final$ID, geno$ID),]
final <- final %>%
    select(-1) %>%
    mutate_at(c(1:21,23), as.numeric)
res <- summary(lm(splicing~., data=final))
geno$residuals <- res$residuals

## Set genotype levels
geno$snp <- factor(geno$snp, level = c(WT, HET, HOM))
WTnum <- paste0(WT, "\n(n=", table(geno$snp)[[1]], ")")
HETnum <- paste0(HET, "\n(n=", table(geno$snp)[[2]], ")")
HOMnum <- paste0(HOM, "\n(n=", table(geno$snp)[[3]], ")")
geno$snp <- ifelse(
    geno$snp == WT,
    WTnum,
    ifelse(
        geno$snp == HET,
        HETnum,
        HOMnum
    )
)
geno$snp <- factor(geno$snp, level = c(WTnum, HETnum, HOMnum))


## Load QTL results
temp <- arrow::read_parquet("input/262sample_sqtl.cis_qtl_pairs.1.parquet") %>%
    filter(phenotype_id == "1:155289479:155290200:clu_10554:ENSG00000160752.14_2") %>%
    filter(variant_id == "1:155289545:T:G")
Slope <- round(temp$slope, 3)
Pvalue <- format(temp$pval_nominal, digits = 3)


## Generate the boxplot of sQTL
sQTLboxp0 <- ggplot(data = geno, aes(x = snp, y = residuals)) +
    geom_boxplot(
        aes(color = snp), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = snp), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#F7BC9B", "#F5A57A", "#F08E59"), name = "") +
    xlab("rs11264361 genotype") +
    ylab("Normalized intron usage") +
    ggtitle(
        bquote(italic(FDPS)~.(IntronID)),
        bquote(italic(beta)*" = "*.(Slope)~", "~italic("P")*" = "*.(Pvalue))
    )
sQTLboxp1 <- sQTLboxp0 +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none"
    )



##------------------------------Fig.5A\B\D

integPlot <- cowplot::plot_grid(GWASp1, sQTLp1, eQTLp1, ncol = 1, align = "v")
integPlot2 <- cowplot::plot_grid(sQTLboxp1, eQTLboxp1, ncol = 1, align = "v")
finalPlot <- cowplot::plot_grid(
    integPlot, integPlot2,
    ncol = 2, nrow = 1, align = "none",
    rel_widths = c(2.5, 1)
)
ggsave("output/Figure5ABD.pdf", finalPlot, width = 8, height = 5)



##------------------------------Fig.5C

## Load list of samples
list <- openxlsx::read.xlsx("input/262sample.info.xlsx") %>%
    select(ID, bamID)


## Load splicing data and calculate the percentage of each intron
psin <- fread("input/262sample.Stomach_perind.filtered.counts.txt.gz", data.table = F) %>%
    filter(grepl("clu_10554_NA", chrom)) %>%
    mutate(chrom = gsub(":clu_10554_NA", "", chrom)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "bamID") %>%
    mutate(`1:155289479:155289585` = tstrsplit(`1:155289479:155289585`, "/")[[1]]) %>%
    mutate(`1:155289479:155289648` = tstrsplit(`1:155289479:155289648`, "/")[[1]]) %>%
    mutate(`1:155289479:155290200` = tstrsplit(`1:155289479:155290200`, "/")[[1]]) %>%
    mutate(`1:155289719:155290200` = tstrsplit(`1:155289719:155290200`, "/")[[1]]) %>%
    mutate_at(2:5, as.numeric) %>%
    mutate(all = `1:155289479:155289585`+`1:155289479:155289648`+`1:155289479:155290200`+`1:155289719:155290200`) %>%
    mutate(`1:155289479:155289585` = `1:155289479:155289585`/all) %>%
    mutate(`1:155289479:155289648` = `1:155289479:155289648`/all) %>%
    mutate(`1:155289479:155290200` = `1:155289479:155290200`/all) %>%
    mutate(`1:155289719:155290200` = `1:155289719:155290200`/all) %>%
    select(-6)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr1.fst") %>%
    tidyfst::filter_fst(grepl("1:155289545", ID)) %>%
    select(-c(1,2,4:9)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    rename(snp = 2) %>%
    mutate(
        snp = case_when(
            snp == "0/0" ~ "GG",
            snp %in% c("0/1", "1/0") ~ "GT",
            snp == "1/1" ~ "TT",
            .default = NA
        )
    )

psd <- list %>%
    inner_join(psin, by = "bamID") %>%
    inner_join(geno, by = "ID") %>%
    select(7,3:6)

psda <- data.frame()
for (g in unique(psd$snp)) {
	sub_psd <- subset(psd, snp == g)
	sub_psda <- data.frame(
        rs11264361 = g, 
		chr1.155289479.155289585 = mean(sub_psd$`1:155289479:155289585`),
        chr1.155289479.155289648 = mean(sub_psd$`1:155289479:155289648`),
		chr1.155289479.155290200 = mean(sub_psd$`1:155289479:155290200`),
        chr1.155289719.155290200 = mean(sub_psd$`1:155289719:155290200`))
	psda <- rbind(psda, sub_psda)
}

tpsda <- t(psda)
write.table(tpsda, "temp/FDPS.subclu.meanExp.txt", col.names=F, quote=F, sep="\t")

psd <- read.delim("temp/FDPS.subclu.meanExp.txt")
colnames(psd)[1] <- "id"

cluster <- data.frame(id = colnames(psda)[2:5]) %>%
    mutate(Start = tstrsplit(id, "\\.")[[2]]) %>%
    mutate(End = tstrsplit(id, "\\.")[[3]]) %>%
    mutate_at(2:3, as.numeric)

clu1 <- merge(cluster, psd[,c("id", "GT")], by = "id")
clu2 <- merge(cluster, psd[,c("id", "GG")], by = "id")
clu3 <- merge(cluster, psd[,c("id", "TT")], by = "id")

colnames(clu1)[4] <- "MeanExp"
clu1$rs11264361 <- "GT"
clu1$yorder <- 3
colnames(clu2)[4] <- "MeanExp"
clu2$rs11264361 <- "GG"
clu2$yorder <- 5
colnames(clu3)[4] <- "MeanExp"
clu3$rs11264361 <- "TT"
clu3$yorder <- 1
clu_use <- rbind(clu1, clu2, clu3)

clu_use$up <- ifelse(clu_use$id=="chr1.155289479.155289648", "Yes", "No")
clu_use$N <- ifelse(
    clu_use$rs11264361 == "TT", 8,
    ifelse(clu_use$rs11264361 == "GT", 81,
        ifelse(clu_use$rs11264361 == "GG", 173, NA)
    )
)
clu_use$ymin <- clu_use$yorder - 0.1
clu_use$ymax <- clu_use$yorder + 0.1

theme_use <- ggplot() +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    ) + 
    xlab("") + ylab("") +
    theme(legend.position = "none")

exonData <- data.frame(
    Start = rep(c(min(clu_use$Start)-(max(clu_use$Start)-min(clu_use$Start))*0.1, 155289585, 155290200), 4),
    End = rep(c(155289479, 155289719, max(clu_use$End)+(max(clu_use$End)-min(clu_use$Start))*0.1), 4),
    ymin = clu_use$ymin,
    ymax = clu_use$ymax
)
genomeData <- data.frame(
    Start = rep(min(clu_use$Start)-(max(clu_use$Start)-min(clu_use$Start))*0.1, 3),
    End = rep(max(clu_use$End)+(max(clu_use$End)-min(clu_use$Start))*0.1, 3),
    yorder = unique(clu_use$yorder)
)


## Generate plot
# structure
p0 <- theme_use + 
    geom_rect(
        data = exonData,
	    aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax)
    ) +
	geom_segment(
        data = genomeData,
	    aes(x = Start, y = yorder, xend = End, yend = yorder),
        linewidth = 1
    ) +
    geom_curve(
        data = subset(clu_use, id == "chr1.155289479.155289585"), 
		aes(x = Start, y = yorder, xend = End, yend = yorder, color = MeanExp, size = MeanExp),
		curvature = -0.4, alpha = 1, lineend = "butt"
    ) +
    geom_curve(
        data = subset(clu_use, id == "chr1.155289479.155289648"),
        aes(x = Start, y = yorder, xend = End, yend = yorder, color = MeanExp, size = MeanExp),
        curvature = 0.6, alpha = 1, lineend = "butt"
    ) +
	geom_curve(
        data = subset(clu_use, id == "chr1.155289479.155290200"),
        aes(x = Start, y = yorder, xend = End, yend = yorder, color = MeanExp, size = MeanExp),
        curvature = -0.6, alpha = 1, lineend = "butt"
    ) +
    geom_curve(
        data = subset(clu_use, id == "chr1.155289719.155290200"),
        aes(x = Start, y = yorder, xend = End, yend = yorder, color = MeanExp, size = MeanExp),
        curvature = -0.3, alpha = 1, lineend = "butt"
    ) +
    ylim(0, 6.5) +
    scale_size(range = c(0.5, 1.5)) +
    scale_color_gradient(low = "#F9D2BD", high = "#F08E59", name = "")

text_x <- min(clu_use$Start) - (max(clu_use$End)-min(clu_use$Start))*0.1
p1 <- p0 +
    geom_text(
        data = subset(clu_use, !duplicated(rs11264361)),
        aes(x = text_x-50, y = yorder, label = paste0(rs11264361, "\n(n=", N, ")")),
        size = 8
    ) +
    geom_text(
        data = subset(clu_use, id == "chr1.155289479.155289585"),
		aes(x = 0.5*(Start+End)-40, y = yorder+0.45, label = paste0(round(MeanExp*100,2), "%")),
        size = 8
    ) +
    geom_text(
        data = subset(clu_use, id == "chr1.155289479.155289648"),
        aes(x = 0.5*(Start+End)+40, y = yorder-0.45, label = paste0(round(MeanExp*100,2), "%")),
        size = 8
    ) +
    geom_text(
        data = subset(clu_use, id == "chr1.155289479.155290200"),
		aes(x = 0.5*(Start+End), y = yorder+0.85, label = paste0(round(MeanExp*100,2), "%")),
        size = 8
    ) +
    geom_text(
        data = subset(clu_use, id == "chr1.155289719.155290200"),
		aes(x = 0.5*(Start+End), y = yorder+0.45, label = paste0(round(MeanExp*100,2), "%")),
        size = 8
    ) +
    xlim(text_x-100, max(clu_use$End)+(max(clu_use$End)-min(clu_use$Start))*0.1)

# SNP
p2 <- p1 + 
    geom_segment(
	    aes(x = 155289545, y = 0.6, xend = 155289545, yend = 5.5),
        color = "#666666", size = 1, linetype = "dotted"
    ) +		
	geom_segment(
	    aes(x = 155289545, y = 0.6, xend = 155289545+80, yend = 0.3),
        color = "#666666", size = 1, linetype = "dotted"
    ) +
    geom_text(
        aes(x = 155289545+80, y = 0.05, label = "Associated SNP:\nrs11264361"),
        size = 7
    ) +
    geom_segment(
        aes(x = 155290200-40, y = 1.6, xend = 155290200-80, yend = 1.4),
        arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ) +			  
    geom_text(
        aes(x = 155290200+20, y = 1.75, label = "Associated intron"),
        size = 7
    )

# junction
p3 <- p2 +
    geom_segment(
	    aes(x = 155289479, y = 5, xend = 155289479, yend = 5.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 155289479 + 60, y = 6.15, label = "chr1:155,289,479"),
        size = 6, angle = 45
    ) +
    geom_segment(
	    aes(x = 155289585, y = 5, xend = 155289585, yend = 5.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 155289585 + 60, y = 6.15, label = "chr1:155,289,585"),
        size = 6, angle = 45
    ) +
    geom_segment(
        aes(x = 155289648, y = 5, xend = 155289648, yend = 5.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 155289648 + 60, y = 6.15, label = "chr1:155,289,648"),
        size = 6, angle = 45
    ) +
    geom_segment(
	    aes(x = 155289719, y = 5, xend = 155289719, yend = 5.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 155289719 + 60, y = 6.15, label = "chr1:155,289,719"),
        size = 6, angle = 45
    ) +
    geom_segment(
	    aes(x = 155290200, y = 5, xend = 155290200, yend = 5.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 155290200 + 60, y = 6.15, label = "chr1:155,290,200"),
        size = 6, angle = 45
    )
ggsave("output/Figure5C.pdf", p3, width = 12, height = 15)