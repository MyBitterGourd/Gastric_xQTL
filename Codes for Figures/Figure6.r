##########################################################
# Code information
# Title: Code of Figure 6
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


## ------------------------------ Fig.6A top

## Load apaTWAS results
twas <- fread("input/apatwas.chr.sig.Fusion.FDR10.list") %>%
    filter(gene_name == "UBE2L3")


## Set the region around the apaTWAS signal
upper <- twas$Start
lower <- twas$End
min.limit <- mean(c(lower,upper)) - 68000
max.limit <- mean(c(lower,upper)) + 68000


## Load conditional apaTWAS results
dat <- read.csv("input/apatwas.Fusion.FDR10.chr22_GWAS.condition.apa.csv") %>%
    drop_na(pv) %>%
    mutate(pos = as.numeric(tstrsplit(snp, ":")[[2]])) %>%
    filter(pos > min.limit) %>%
    filter(pos < max.limit)


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
        values = rev(c("grey", "#68B69F")),
        name = ""
    ) +
    xlab("Chr22 physical position (Mb)") +
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



## ------------------------------ Fig.6A middle

## Load SNPs in the region and calculate MAF
dat <- fread("input/apaQTL.TWAS.allpairs.txt") %>%
    filter(phenotype_id == "ENST00000458578.6_1|ENSG00000185651.14_2|chr22|+")
fwrite(dat[, "variant_id"], "temp/chr22.snp.list", col=F, row=F, quo=F, sep="\t")

system("plink --bfile input/geno --extract temp/chr22.snp.list --make-bed --out temp/chr22")
system("plink --bfile temp/chr22 --freq --out temp/chr22")


## Load apaQTL results
freq <- fread("temp/chr22.frq")
all <- dat %>%
    inner_join(freq, by = c("variant_id" = "SNP")) %>%
    select(variant_id, A1, A2, MAF, slope, slope_se, pval_nominal) %>%
    rename(SNP = variant_id, freq = MAF, b = slope, se = slope_se, p = pval_nominal) %>%
    mutate(N = 262)
fwrite(all, "temp/chr22.ma", col=T, row=F, quo=F, sep="\t")
write.table("22:21977047:C:T", "temp/condition.snplist.chr22", row=F, col=F, quo=F, sep="\t")


## Perform conditional regression
system("gcta64 --bfile temp/chr22 --cojo-file temp/chr22.ma --cojo-cond temp/condition.snplist.chr22 --out temp/condition.chr22.TopaQTL")
system("rm temp/chr22.*")


## prepare output for mPLOT!
con <- fread("temp/condition.chr22.TopaQTL.cma.cojo") %>%
    select(SNP, bC, bC_se, pC)
final <- all %>%
    left_join(con, by = "SNP")
fwrite(final, "temp/condition.chr22.TopaQTL.txt", col=T, row=F, quo=F, sep="\t")


## Load apaTWAS results
twas <- fread("input/apatwas.chr.sig.Fusion.FDR10.list") %>%
    filter(ID == "ENST00000458578.6_1|ENSG00000185651.14_2|chr22|+")


## Set the region around the apaTWAS signal
upper <- twas$Start
lower <- twas$End
min.limit <- mean(c(lower,upper)) - 68000
max.limit <- mean(c(lower,upper)) + 68000
dat <- final %>%
    mutate(pos = as.numeric(tstrsplit(SNP, ":")[[2]])) %>%
    filter(pos>min.limit & pos<max.limit)


## Separate the results and then combine them
rawgwas <- dat %>%
    select(pos, p) %>%
    mutate(note = "Main apaQTL")
condgwas <- dat %>%
    select(pos, p = pC) %>%
    mutate(note = "Conditioned on\ntop apaQTL")
dat1 <- rbind(rawgwas,condgwas)


## Set top GWAS-SNP
sda <- dat %>%
    filter(SNP == "22:21980257:A:G") %>%
    mutate(rsid = "rs5754422")


## Generate the middle plot of apaQTL
apaQTLp0 <- ggplot(data = dat1, aes(x = pos/1e6, y = -log10(p))) +
    geom_point(
        aes(color = note), pch = 16,
        size = 1.5, alpha = 1 
    ) +
    scale_color_manual(
        values = rev(c("grey","#68B69F")),
        name = ""
    ) +
    ylab(bquote(~atop(italic(UBE2L3)~"apaQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    geom_text(
        data = sda,
        aes(x = pos/1e6, y = -log10(p)+0.7, label = rsid)
    )
apaQTLp1 <- apaQTLp0 + theme(
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



## ------------------------------ Fig.6A bottom

## Load eQTL results
dat <- arrow::read_parquet("input/262sample_eqtl.cis_qtl_pairs.22.parquet") %>%
    filter(grepl("ENSG00000185651", phenotype_id)) %>%
    select(variant_id, pval_nominal) %>%
    mutate(pos = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    select(pos, p = pval_nominal) %>%
    mutate(note = "Main eQTL")


## Load apaTWAS results
twas <- fread("input/apatwas.chr.sig.Fusion.FDR10.list") %>%
    filter(ID == "ENST00000458578.6_1|ENSG00000185651.14_2|chr22|+")


## Set the region around the apaTWAS signal
upper <- twas$Start
lower <- twas$End
min.limit <- mean(c(lower,upper)) - 68000
max.limit <- mean(c(lower,upper)) + 68000
dat <- dat %>%
    filter(pos>min.limit & pos<max.limit)


## Generate the bottom plot of eQTL
eQTLp0 <- ggplot(data = dat, aes(x = pos/1e6, y = -log10(p))) +
    geom_point(
        aes(color = note), pch = 16,
        size = 1.5, alpha = 1
    ) +
    scale_color_manual(values = c("black"), name = "") +
    xlab("Chr22 physical position (Mb)") +
    ylab(bquote(~atop(italic(UBE2L3)~"eQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(
        breaks = seq(0, 15, 5),
        labels = seq(0, 15, 5),
        limits = c(0, 17),
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

p <- GWASp1 / apaQTLp1 / eQTLp1
ggsave("output/Figure6A.pdf", p, width = 8, height = 8)



## ------------------------------ Fig.6B left

## Set transcript ID
TranscriptID <- tstrsplit("ENST00000458578.6_1|ENSG00000185651.14_2|chr22|+", "\\|")[[1]]
TranscriptID <- tstrsplit(TranscriptID, "\\.")[[1]]


## Set alleles
REF <- tstrsplit("22:21977047:C:T", ":")[[3]]
ALT <- tstrsplit("22:21977047:C:T", ":")[[4]]
WT <- paste0(REF, REF)
HET <- paste0(REF, ALT)
HOM <- paste0(ALT, ALT)


## Load PDUI data
apa <- fread("input/262sample_Stomach.apa.bed.gz") %>%
    filter(Gene == "ENST00000458578.6_1|ENSG00000185651.14_2|chr22|+") %>%
    select(-c(1:3)) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    rename(apa = 2)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr22.fst") %>%
    tidyfst::filter_fst(ID == "22:21977047:C:T") %>%
    select(-c(1,2,4:9)) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
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
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID")


## Merge data and calculate residuals
final <- apa %>%
    inner_join(covar, by = "ID")
geno <- geno[match(final$ID, geno$ID),]
final <- final %>%
    select(-1) %>%
    mutate_at(c(1:6,8,10:44), as.numeric)
res <- summary(lm(apa~., data=final))
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


## Load apaQTL results
temp <- tidyfst::parse_fst("input/262sample_aqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::filter_fst(phenotype_id == "ENST00000458578.6_1|ENSG00000185651.14_2|chr22|+") %>%
    tidyfst::filter_fst(variant_id == "22:21977047:C:T")
Slope <- round(temp$slope, 3)
Pvalue <- format(temp$pval_nominal, digits = 3)


## Generate boxplot
apaQTLboxp0 <- ggplot(data = geno, aes(x = snp, y = residuals)) +
    geom_boxplot(
        aes(color = snp), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = snp), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#A3D4C6", "#84C5B2", "#68B69F"), name = "") +
    xlab("rs7445 genotype") +
    ylab("Normalized PDUI") +
    ggtitle(
        bquote(italic(UBE2L3)~.(TranscriptID)),
        bquote(italic(beta)*" = "*.(Slope)~", "~italic("P")*" = "*.(Pvalue))
    )
apaQTLboxp1 <- apaQTLboxp0 +
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
ggsave("output/Figure6Bleft.pdf", apaQTLboxp1, width = 4, height = 4)



## ------------------------------ Fig.6B right

## Load annotation data for gene
library(ggtranscript)
gtf <- rtracklayer::import("input/gencode.v29lift37.annotation.gtf")

UBE2L3_gtf <- Repitools::annoGR2DF(gtf) %>%
    dplyr::rename(seqnames = chr) %>%
    filter(transcript_id == "ENST00000458578.6_1")
min(UBE2L3_gtf$start)
# [1] 21903736
max(UBE2L3_gtf$end)
# [1] 21978323
UBE2L3_exons <- UBE2L3_gtf %>%
    filter(type == "exon")
UBE2L3_exons_cod <- UBE2L3_exons %>%
    filter(transcript_type == "protein_coding")
UBE2L3_cds <- UBE2L3_gtf %>%
    filter(type == "CDS")


## Generate gene plot
genePlot <- UBE2L3_exons_cod %>%
    ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_name
    )) +
    geom_range(
        fill = "black",
        height = 0.25
    ) +
    geom_range(
        data = UBE2L3_cds,
        fill = "black"
    ) +
    geom_intron(
        data = to_intron(UBE2L3_exons_cod, "gene_name"),
        aes(strand = strand),
        arrow.min.intron.length = 1000
    ) +
    geom_vline(xintercept = 21977047, color = "#BD514A", linewidth = 0.5, linetype = "dashed") +
    xlim(21965146, 21978323) +
    xlab("Chr22 physical position (Mb)") +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none",
        plot.title = element_blank()
    )


## Generate track plot
region <- as("chr22:21965146-21978323","GRanges")

# WT
WTdata <- rtracklayer::import(
    con = "input/WT.chr22.bigwig",
    which = as("chr22:21965146-21978323","GRanges"),
    as = "NumericList"
)[[1]]
WTdata <- data.frame(
    position = start(x = region):end(x = region),
    score = WTdata,
    stringsAsFactors = FALSE
)
WTplot <- ggplot(
        data = WTdata,
        mapping = aes_string(x = "position", y = "score")
    ) +
    geom_bar(stat="identity",fill="#68B69F",color="#68B69F") +
    ylab("CC") +
    theme(
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
        legend.position = "none",
        plot.title = element_blank()
    )

# HET
HETdata <- rtracklayer::import(
    con = "input/HET.chr22.bigwig",
    which = as("chr22:21965146-21978323","GRanges"),
    as = "NumericList"
)[[1]]
HETdata <- data.frame(
    position = start(x = region):end(x = region),
    score = HETdata,
    stringsAsFactors = FALSE
)
HETplot <- ggplot(
        data = HETdata,
        mapping = aes_string(x = "position", y = "score")
    ) +
    geom_bar(stat="identity",fill="#68B69F",color="#68B69F") +
    ylab("CT") +
    theme(
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
        legend.position = "none",
        plot.title = element_blank()
    )

# HOM
HOMdata <- rtracklayer::import(
    con = "input/HOM.chr22.bigwig",
    which = as("chr22:21965146-21978323","GRanges"),
    as = "NumericList"
)[[1]]
HOMdata <- data.frame(
    position = start(x = region):end(x = region),
    score = HOMdata,
    stringsAsFactors = FALSE
)
HOMplot <- ggplot(
        data = HOMdata,
        mapping = aes_string(x = "position", y = "score")
    ) +
    geom_bar(stat="identity",fill="#68B69F",color="#68B69F") +
    ylab("TT") +
    theme(
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
        legend.position = "none",
        plot.title = element_blank()
    )


## Combine plots
p <- WTplot / HETplot / HOMplot / genePlot + plot_layout(heights = c(1.8, 1.8, 1.8, 1))
ggsave("output/Figure6Bright", p, width = 4, height = 4)



## ------------------------------ Fig.6C

## Load GWAS results
gwas_geno <- fread("input/GC.summary.txt")


## Load apaQTL results
gene_list <- fread("input/apatwas.chr.sig.Fusion.results_adj") %>%
    filter(TWAS.FDR < 0.1)
apaqtl_geno <- tidyfst::parse_fst("input/262sample_aqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::filter_dt(phenotype_id %in% gene_list$ID)


## Combine GWAS and apaQTL data
coloc_data <- apaqtl_geno %>%
    inner_join(gwas_geno, by = "SNP") %>%
    mutate(chr = tstrsplit(SNP, ":")[[1]]) %>%
    mutate(bp = tstrsplit(SNP, ":")[[2]]) %>%
    mutate(chrbp = paste0(chr, ":", bp)) %>%
    inner_join(snp, by = "chrbp") %>%
    mutate(pval_sum = pval_gwas + pval_apaqtl)


## Generate plot
library(locuscomparer)
library(ggtext)

cortheme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.ticks.length.y = unit(-0.15, "cm"),
    axis.title.y = element_markdown(size = 8, color = "black"),
    axis.title.x = element_markdown(size = 8, color = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black")
)

geneName <- "UBE2L3"
gene <- subset(gene_list, gene_name == geneName)$ID

sub <- coloc_data %>%
    filter(phenotype_id == gene) %>%
    mutate_at(5:6, as.integer)
rs <- sub[which.min(sub$pval_sum),]$rsID[1]
bpmid <- sub[which.min(sub$pval_sum),]$bp[1]
chr <- sub$chr[1]
ld <- retrieve_LD(chr, rs, "EAS") %>% mutate(rsID = SNP_B)
sub <- sub %>%
    left_join(ld, by = "rsID") %>%
    mutate(R2 = ifelse(is.na(R2), 0, R2)) %>%
    mutate(R2 = ifelse(rsID == rs, NA, R2)) %>%
    mutate(label = ifelse(rsID == rs, rsID, "")) %>%
    mutate(size = ifelse(rsID == rs, 6, 4)) %>%
    mutate(
        R2level = case_when(
            R2 <= 0.2 ~ "level 1",
            R2 > 0.2 & R2 <= 0.4 ~ "level 2",
            R2 > 0.4 & R2 <= 0.6 ~ "level 3",
            R2 > 0.6 & R2 <= 0.8 ~ "level 4",
            R2 > 0.8 ~ "level 5"
        )
    ) %>%
    mutate(R2level = factor(R2level, levels = c("level 1", "level 2", "level 3", "level 4", "level 5")))
RValue <- cor.test(-log10(sub$pval_apaqtl), -log10(sub$pval_gwas), exact=F)$estimate
PValue <- cor.test(-log10(sub$pval_apaqtl), -log10(sub$pval_gwas), exact=F)$p.value
labels <- paste("Pearson correlation", "\nr =", format(RValue, digits = 3),"\nP =", format(PValue, scientific = TRUE, digits = 3))

p <- ggplot(sub, aes(x = -log10(pval_gwas), y = -log10(pval_apaqtl))) +
    geom_point(
        data = subset(sub, label != rs),
        aes(x = -log10(pval_gwas), y = -log10(pval_apaqtl), color = R2level, size = size),
        alpha = 1, shape = 19
    ) +
    scale_color_manual(
        values = c("#2B5599", "#9EC4E6", "#47C3C6", "#FDC107", "#F06598"),
        guide = "none"
    ) +
    scale_size(range = c(1,2), guide = "none") +
    geom_point(
        data = subset(sub, label == rs),
        aes(x = -log10(pval_gwas), y = -log10(pval_apaqtl), size = size),
        alpha = 1, shape = 21, color = "black", fill = "#EF3A1E"
    ) +
    labs(
        x = "GWAS -log<sub>10</sub>(P)",
        y = paste0("<i>", geneName, "</i> apaQTL<br>-log<sub>10</sub>(P)")
    ) +
    ggrepel::geom_text_repel(
        aes(label = label), size = 3, direction = "both", segment.color = NA,
        point.padding = 1, max.overlaps = getOption("ggrepel.max.overlaps", default = 50)
    ) +
    annotate(
        "text", x = 0, y = Inf, label = labels,
        hjust = 0, vjust = 1.3, size = 3
    ) +
    cortheme
ggsave("output/Figure6C.pdf", p, width = 4, height = 4)