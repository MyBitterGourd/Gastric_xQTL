##########################################################
# Code information
# Title: Code of Figure 1
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################


## Load environment
R
options(stringsAsFactors=F)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggvenn)
library(ggrastr)
library(ggprism)
library(cowplot)
rm(list=ls())


## ------------------------------ MARK: Fig.1A

## Load data
eGene <- fread("input/262sample_eqtl.cis_qtl.perm.list") %>%
    count(phenotype_id)

sGene <- fread("input/262sample_sqtl.cis_qtl.perm.list") %>%
    mutate(phenotype_id = tstrsplit(group_id, "_")[[1]]) %>%
    count(phenotype_id)

apaGene <- fread("input/262sample_aqtl.cis_qtl.perm.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\|")[[2]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "_")[[1]]) %>%
    count(phenotype_id)

Gene <- list(eGene$phenotype_id, sGene$phenotype_id, apaGene$phenotype_id)
names(Gene) <- c(
    paste0("eGene\n(", nrow(eGene), ")"),
    paste0("sGene\n(", nrow(sGene), ")"),
    paste0("apaGene\n(", nrow(apaGene), ")")
)


## Generate Venn plot
Vennplot <- ggvenn(
    Gene, show_percentage = F,
    stroke_color = "white",
    fill_color = c("#145390", "#F08E59", "#68B69F"),
    fill_alpha = 0.8,
    set_name_color = c("#145390", "#F08E59", "#68B69F"),
    text_color = "black"
)
ggsave("output/QTL.Genes.vennplot.pdf", Vennplot, width = 4, height = 4)



##------------------------------Fig.1B

## Load data and calculate the percentage of independent signals
rm(list=ls())
eqtl <- fread("input/262sample_eqtl.cis_independent_qtl.txt.gz") %>%
    count(phenotype_id)
eqtlStats <- data.frame(
    qtl = rep("eQTL", 2),
    group = c("single SNP", ">=2 SNPs"),
    value = c(
        nrow(eqtl[eqtl$n==1,]),
        nrow(eqtl[eqtl$n>=2,])
    )
)
eqtlStats$percentage <- eqtlStats$value / nrow(eqtl) * 100

sqtl <- fread("input/262sample_sqtl.cis_independent_qtl.txt.gz") %>%
    count(group_id)
sqtlStats <- data.frame(
    qtl = rep("sQTL", 2),
    group = c("single SNP", ">=2 SNPs"),
    value = c(
        nrow(sqtl[sqtl$n==1,]),
        nrow(sqtl[sqtl$n>=2,])
    )
)
sqtlStats$percentage <- sqtlStats$value / nrow(sqtl) * 100

aqtl <- fread("input/262sample_aqtl.cis_independent_qtl.txt.gz") %>%
    count(group_id)
aqtlStats <- data.frame(
    qtl = rep("apaQTL", 2),
    group = c("single SNP", ">=2 SNPs"),
    value = c(
        nrow(aqtl[aqtl$n==1,]),
        nrow(aqtl[aqtl$n>=2,])
    )
)
aqtlStats$percentage <- aqtlStats$value / nrow(aqtl) * 100

Stats <- eqtlStats %>%
    rbind(sqtlStats) %>%
    rbind(aqtlStats) %>%
    mutate(qtl = factor(qtl, levels = c("eQTL", "sQTL", "apaQTL")))

mytheme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.ticks.length.y = unit(-0.15, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none"
)


## Generate barplot
p <- ggplot(Stats, aes(qtl, value), position = "stack") +
    scale_x_discrete(limits = c("eQTL", "sQTL", "apaQTL")) +
    geom_bar(
        aes(fill = group), stat = "identity", color = "white",
        size = 0.4, position = "fill", width = 0.6, data = Stats
    ) +
    scale_fill_manual(values = c("#A1BAD3", "#4375A6")) +
    ylab("Percentage of QTL-Genes") +
    mytheme
ggsave("output/independent.QTLnumber.Barplot.pdf", p, width = 4, height = 3)


## Chi-square test
mytable <- matrix(Stats$value[c(1:4)], nrow = 2, ncol = 2)
chisq.test(mytable)

mytable <- matrix(Stats$value[c(1:2, 5:6)], nrow = 2, ncol = 2)
chisq.test(mytable)

mytable <- matrix(Stats$value[c(3:6)], nrow = 2, ncol = 2)
chisq.test(mytable)




##------------------------------Fig.1C

## Load data and calculate the percentage of each annotation
rm(list=ls())

# eQTL
eQTL <- fread("input/eSNP.severe.txt") %>%
    mutate(Count = gsub(",","",Count)) %>%
    mutate_at(2, as.numeric) %>%
    mutate(`Consequence type` = ifelse(
        `Consequence type` %in% c("splice_donor_variant", "splice_acceptor_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_polypyrimidine_tract_variant", "splice_donor_region_variant"),
        "Splice\nregion",
        ifelse(
            `Consequence type` == "5_prime_UTR_variant",
            "5'UTR",
            ifelse(
                `Consequence type` == "3_prime_UTR_variant",
                "3'UTR",
                ifelse(
                    `Consequence type` == "intron_variant",
                    "Intron",
                    ifelse(
                        `Consequence type` == "upstream_gene_variant",
                        "Upstream\ngene",
                        ifelse(
                            `Consequence type` == "downstream_gene_variant",
                            "Downstream\ngene",
                            ifelse(
                                `Consequence type` == "intergenic_variant",
                                "Intergenic",
                                "Other"
                            )
                        )
                    )
                )
            )
        )
    )) %>%
    group_by(`Consequence type`) %>%
    reframe(N = sum(Count)) %>%
    mutate(QTL = "eQTL") %>%
    reframe(QTL, `Consequence type`, N, Percentage = N / sum(N) * 100, Total = sum(N))

# sQTL
sQTL <- fread("input/sSNP.severe.txt") %>%
    mutate(Count = gsub(",","",Count)) %>%
    mutate_at(2, as.numeric) %>%
    mutate(`Consequence type` = ifelse(
        `Consequence type` %in% c("splice_donor_variant", "splice_acceptor_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_polypyrimidine_tract_variant", "splice_donor_region_variant"),
        "Splice\nregion",
        ifelse(
            `Consequence type` == "5_prime_UTR_variant",
            "5'UTR",
            ifelse(
                `Consequence type` == "3_prime_UTR_variant",
                "3'UTR",
                ifelse(
                    `Consequence type` == "intron_variant",
                    "Intron",
                    ifelse(
                        `Consequence type` == "upstream_gene_variant",
                        "Upstream\ngene",
                        ifelse(
                            `Consequence type` == "downstream_gene_variant",
                            "Downstream\ngene",
                            ifelse(
                                `Consequence type` == "intergenic_variant",
                                "Intergenic",
                                "Other"
                            )
                        )
                    )
                )
            )
        )
    )) %>%
    group_by(`Consequence type`) %>%
    reframe(N = sum(Count)) %>%
    mutate(QTL = "sQTL") %>%
    reframe(QTL, `Consequence type`, N, Percentage = N / sum(N) * 100, Total = sum(N))

# apaQTL
apaQTL <- fread("input/aSNP.severe.txt") %>%
    mutate(Count = gsub(",","",Count)) %>%
    mutate_at(2, as.numeric) %>%
    mutate(`Consequence type` = ifelse(
        `Consequence type` %in% c("splice_donor_variant", "splice_acceptor_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_polypyrimidine_tract_variant", "splice_donor_region_variant"),
        "Splice\nregion",
        ifelse(
            `Consequence type` == "5_prime_UTR_variant",
            "5'UTR",
            ifelse(
                `Consequence type` == "3_prime_UTR_variant",
                "3'UTR",
                ifelse(
                    `Consequence type` == "intron_variant",
                    "Intron",
                    ifelse(
                        `Consequence type` == "upstream_gene_variant",
                        "Upstream\ngene",
                        ifelse(
                            `Consequence type` == "downstream_gene_variant",
                            "Downstream\ngene",
                            ifelse(
                                `Consequence type` == "intergenic_variant",
                                "Intergenic",
                                "Other"
                            )
                        )
                    )
                )
            )
        )
    )) %>%
    group_by(`Consequence type`) %>%
    reframe(N = sum(Count)) %>%
    mutate(QTL = "apaQTL") %>%
    reframe(QTL, `Consequence type`, N, Percentage = N / sum(N) * 100, Total = sum(N))

result <- eQTL %>% rbind(sQTL) %>% rbind(apaQTL)


## Fisher test
fisherResult <- data.frame()
for (testQTL in c("sQTL", "apaQTL")) {
    for (testAnno in c("Intergenic", "Downstream\ngene", "Upstream\ngene", "Intron", "3'UTR", "5'UTR", "Splice\nregion")) {
        print(paste0("Processing: ", testQTL, " in ", testAnno))

        testSub <- subset(result, QTL == testQTL & `Consequence type` == testAnno)
        controlSub <- subset(result, QTL == "eQTL" & `Consequence type` == testAnno)

        mytable <- matrix(
            c(
                testSub$N,
                controlSub$N,
                testSub$Total-testSub$N,
                controlSub$Total-controlSub$N
            ),
            nrow = 2
        )
        colnames(mytable) <- c("assigned", "not assigned")
        rownames(mytable) <- c(testQTL, "eQTL")
        fisherTest <- fisher.test(mytable)
        Pvalue <- fisherTest$p.value
        OR <- unname(fisherTest$estimate)
        Lower <- fisherTest$conf.int[1]
        Upper <- fisherTest$conf.int[2]

        temp <- c(testQTL, testAnno, OR, Lower, Upper, Pvalue)
        fisherResult <- rbind(fisherResult, temp)
    }
}
names(fisherResult) <- c("QTL", "annotation", "OR", "Lower", "Upper", "Pvalue")
fisherResult <- fisherResult %>%
    mutate_at(3:6, as.numeric) %>%
    mutate(OR = round(OR, 2), Lower = round(Lower, 2), Upper = round(Upper, 2))


## set theme
mytheme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(b = 6)),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.line = element_blank()
)


## sQTL vs eQTL
dat <- fisherResult %>%
    filter(QTL == "sQTL") %>%
    arrange(OR)
dat$annotation <- factor(dat$annotation, levels = dat$annotation)

p0 <- ggplot(dat, aes(annotation, OR)) +
    geom_bar(stat = "identity", width = 0.5, fill = "#F08E59") +
    geom_errorbar(
        ymin = dat$Lower,
        ymax = dat$Upper,
        width = 0.2, size = 0.5
    ) +
    scale_y_continuous(limits = c(0, 1.6), position = "right") +
    xlab("") + ylab("OR (sQTLs compared with eQTLs)") +
    geom_hline(yintercept = 1, color = "#BD514A", linetype = 2, size = 0.5) +
    coord_flip()
p1 <- p0 + mytheme
ggsave("output/annotation.sSNP.vs.eSNP.pdf", p1, width = 3, height = 4)


## apaQTL vs eQTL
dat <- fisherResult %>%
    filter(QTL == "apaQTL") %>%
    arrange(OR, Upper)
dat$annotation <- factor(dat$annotation, levels = dat$annotation)

p0 <- ggplot(dat, aes(annotation, OR)) +
    geom_bar(stat = "identity", width = 0.5, fill = "#68B69F") +
    geom_errorbar(
        ymin = dat$Lower,
        ymax = dat$Upper,
        width = 0.2, size = 0.5
    ) +
    scale_y_continuous(limits = c(0, 1.6), position = "right") +
    xlab("") + ylab("OR (apaQTLs compared with eQTLs)") +
    geom_hline(yintercept = 1, color = "#BD514A", linetype = 2, size = 0.5) +
    coord_flip()
p1 <- p0 + mytheme
ggsave("output/annotation.apaSNP.vs.eSNP.pdf", p1, width = 3, height = 4)



##------------------------------Fig.1D

## Load data
rm(list=ls())

# eQTL
eQTLtop <- fread("input/262sample_eqtl.cis_qtl.perm.top.list") %>%
    arrange(phenotype_id, pval_nominal, abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, eSNP = variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]])

# sQTL
sQTLtop <- fread("input/262sample_sqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(group_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, pval_nominal, abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, sSNP = variant_id)

# apaQTL
apaQTLtop <- fread("input/262sample_aqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\|")[[2]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, pval_nominal, abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, apaSNP = variant_id)


## set theme
mytheme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5)
)


## eSNP & sSNP
eSNPsSNPshared <- eQTLtop %>%
    inner_join(sQTLtop, by = "phenotype_id") %>%
    mutate(eSNPbp = tstrsplit(eSNP, ":")[[2]]) %>%
    mutate(sSNPbp = tstrsplit(sSNP, ":")[[2]]) %>%
    mutate_at(4:5, as.numeric) %>%
    mutate(distance = abs(eSNPbp - sSNPbp))

p0 <- ggplot(aes(x = log10(distance + 1)), data = eSNPsSNPshared) + 
	geom_histogram(bins = 25, fill = "#F08E59", color = "white") + 
	xlab("Distance (Kb)") + ylab("Number of Genes") +
    scale_x_continuous(
        breaks = log10(c(0,10,100,1000,10000,100000,1000000) + 1),
        labels = c(0,10,100,1000,10000,100000,1000000)/1e3, expand = c(1/50, 1/50)
    ) +
    ggtitle("Distance in base piars between\nthe lead sQTL and eQTL")
p1 <- p0 + mytheme
ggsave("output/eGene.sGene.shared.distance.pdf", p1, width = 2.5, height = 2.5)


## eSNP & apaSNP
eSNPapaSNPshared <- eQTLtop %>%
    inner_join(apaQTLtop, by = "phenotype_id") %>%
    mutate(eSNPbp = tstrsplit(eSNP, ":")[[2]]) %>%
    mutate(apaSNPbp = tstrsplit(apaSNP, ":")[[2]]) %>%
    mutate_at(4:5, as.numeric) %>%
    mutate(distance = abs(eSNPbp - apaSNPbp))

p0 <- ggplot(aes(x = log10(distance + 1)), data = eSNPapaSNPshared) + 
	geom_histogram(bins = 25, fill = "#68B69F", color = "white") + 
	theme_classic() + xlab("Distance (Kb)") + ylab("Number of Genes") +
    scale_x_continuous(
        breaks = log10(c(0,10,100,1000,10000,100000,1000000) + 1),
        labels = c(0,10,100,1000,10000,100000,1000000)/1e3, expand = c(1/50, 1/50)
    ) +
    ggtitle("Distance in base piars between\nthe lead apaQTL and eQTL")
p1 <- p0 + mytheme
ggsave("output/eGene.apaGene.shared.distance.pdf", p1, width = 2.5, height = 2.5)



##------------------------------Fig.1E

## eSNP & sSNP
eSNPsSNP <- data.frame(c(eSNPsSNPshared$eSNP, eSNPsSNPshared$sSNP))
fwrite(eSNPsSNP, "temp/eGene.sGene.shared.SNP.list", col=F, row=F, quo=F, sep="\t")
fwrite(data.frame(unique(eSNPsSNPshared$eSNP)), "temp/eGene.sGene.shared.eSNP.list", col=F, row=F, quo=F, sep="\t")

# Calculate LD relationship
system("plink --vcf input/genotype.vcf.gz --extract temp/eGene.sGene.shared.SNP.list --make-bed --out temp/eGene.sGene.shared.SNP") # 1184 variants
system("plink --bfile temp/eGene.sGene.shared.SNP --r2 --ld-window-kb 99999999 --ld-window 99999999 --ld-window-r2 0 --ld-snp-list temp/eGene.sGene.shared.eSNP.list --out temp/eGene.sGene.shared.eSNP")
system("rm temp/*.bed")
system("rm temp/*.bim")
system("rm temp/*.fam")
system("rm temp/*.log")
system("rm temp/*.nosex")

eSNPsSNPLD <- fread("temp/eGene.sGene.shared.eSNP.ld") %>%
    select(SNP_A, SNP_B, R2) %>%
    mutate(snp_pair = paste(SNP_A, SNP_B, sep = "_")) %>%
    select(snp_pair, R2)
eSNPsSNPshared <- eSNPsSNPshared %>%
    mutate(snp_pair = paste(eSNP, sSNP, sep = "_")) %>%
    left_join(eSNPsSNPLD, by = "snp_pair")
fwrite(eSNPsSNPshared, "temp/eGene.sGene.shared.list", col=T, row=F, quo=F, sep="\t")

p0 <- ggplot(aes(x = R2), data = eSNPsSNPshared) +
    geom_histogram(bins = 25, fill = "#F08E59", color = "white") +
    xlab(expression("Linkage disequilibrium (R"^2*")")) +
    ylab("Number of Genes") +
    ggtitle("Linkage disequilibrium between\nthe lead sQTL and eQTL")
p1 <- p0 + mytheme
ggsave("output/eGene.sGene.shared.LD.pdf", p1, width = 2.5, height = 2.5)


## eSNP & apaSNP
eSNPapaSNP <- data.frame(c(eSNPapaSNPshared$eSNP, eSNPapaSNPshared$apaSNP))
fwrite(eSNPapaSNP, "temp/eGene.apaGene.shared.SNP.list", col=F, row=F, quo=F, sep="\t")
fwrite(data.frame(unique(eSNPapaSNPshared$eSNP)), "temp/eGene.apaGene.shared.eSNP.list", col=F, row=F, quo=F, sep="\t")

# Calculate LD relationship
system("plink --vcf input/genotype.vcf.gz --extract temp/eGene.apaGene.shared.SNP.list --make-bed --out temp/eGene.apaGene.shared.SNP") # 485 variants
system("plink --bfile temp/eGene.apaGene.shared.SNP --r2 --ld-window-kb 99999999 --ld-window 99999999 --ld-window-r2 0 --ld-snp-list temp/eGene.apaGene.shared.eSNP.list --out temp/eGene.apaGene.shared.eSNP")
system("rm temp/*.bed")
system("rm temp/*.bim")
system("rm temp/*.fam")
system("rm temp/*.log")
system("rm temp/*.nosex")

eSNPapaSNPLD <- fread("temp/eGene.apaGene.shared.eSNP.ld") %>%
    select(SNP_A, SNP_B, R2) %>%
    mutate(snp_pair = paste(SNP_A, SNP_B, sep = "_")) %>%
    select(snp_pair, R2)
eSNPapaSNPshared <- eSNPapaSNPshared %>%
    mutate(snp_pair = paste(eSNP, apaSNP, sep = "_")) %>%
    left_join(eSNPapaSNPLD, by = "snp_pair")
fwrite(eSNPapaSNPshared, "temp/eGene.apaGene.shared.list", col=T, row=F, quo=F, sep="\t")

p0 <- ggplot(aes(x = R2), data = eSNPapaSNPshared) +
    geom_histogram(bins = 25, fill = "#68B69F", color = "white") +
    xlab(expression("Linkage disequilibrium (R"^2*")")) +
    ylab("Number of Genes") +
    ggtitle("Linkage disequilibrium between\nthe lead apaQTL and eQTL")
p1 <- p0 + mytheme
ggsave("output/eGene.apaGene.shared.LD.pdf", p1, width = 2.5, height = 2.5)



##------------------------------Fig.1F

## set theme
mytheme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
    axis.ticks.length.x = unit(0, "cm"),
    axis.ticks.length.y = unit(0, "cm"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.line = element_blank(),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5)
)


## sQTL in eQTL
topsQTL <- fread("input/262sample_sqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, ":")[[5]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, -abs(slope), abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, variant_id, slope_sqtl = slope)

# Match the results of the specified gene-variant pair in all eQTLs
eQTL <- tidyfst::parse_fst("input/262sample_eqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topsQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_eqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_eqtl) %>%
    arrange(pair, -abs(slope_eqtl)) %>%
    distinct(pair, .keep_all = TRUE)

topsQTLineQTL <- topsQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_sqtl) %>%
    inner_join(eQTL, by = "pair") %>%
    rename(eqtl = slope_eqtl, sqtl = slope_sqtl)
nrow(topsQTLineQTL) # 1239 pairs
fwrite(topsQTLineQTL, "temp/top.sQTL.in.eQTL.slope.list", col=T, row=F, quo=F, sep="\t")

model.lm <- lm(formula = eqtl ~ sqtl, data = topsQTLineQTL)
r <- round(summary(model.lm)$coefficients[2,1], 3)
eq <- substitute(italic(r) ~ "=" ~ x, list(x = r))

p0 <- ggplot(aes(x = sqtl, y = eqtl), data = topsQTLineQTL) +
    ggrastr::geom_point_rast(
        alpha = 0.6, colour = "#F08E59",
        shape = 16, size = 2,
        raster.dpi = getOption("ggrastr.default.dpi", 300)
    ) + 
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "black") +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "black") +
    xlab("sQTL effect size") + ylab("eQTL effect size") +
    geom_smooth(
        method = lm, formula = y ~ x,
        colour = "red", se = FALSE, lwd = 1
    ) +
    annotate(
        "text", x = -2.5, y = 1.5,
        label =  as.character(as.expression(eq)),
        parse = TRUE, size = 3
    )
p1 <- p0 + mytheme
ggsave("output/top.sQTL.in.eQTL.slope.pdf", p1, width = 3, height = 3)


## apaQTL in eQTL
topapaQTL <- fread("input/262sample_aqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\|")[[2]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, -abs(slope), abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, variant_id, slope_apaqtl = slope)

# Match the results of the specified gene-variant pair in all eQTLs
eQTL <- tidyfst::parse_fst("input/262sample_eqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topapaQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_eqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_eqtl) %>%
    arrange(pair, -abs(slope_eqtl)) %>%
    distinct(pair, .keep_all = TRUE)

topapaQTLineQTL <- topapaQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_apaqtl) %>%
    inner_join(eQTL, by = "pair") %>%
    rename(eqtl = slope_eqtl, apaqtl = slope_apaqtl)
nrow(topapaQTLineQTL) # 492 pairs
fwrite(topapaQTLineQTL, "temp/top.apaQTL.in.eQTL.slope.list", col=T, row=F, quo=F, sep="\t")

model.lm <- lm(formula = eqtl ~ apaqtl, data = topapaQTLineQTL)
r <- round(summary(model.lm)$coefficients[2,1], 3)
eq <- substitute(italic(r) ~ "=" ~ x, list(x = r))

p0 <- ggplot(aes(x = apaqtl, y = eqtl), data = topapaQTLineQTL) +
    ggrastr::geom_point_rast(
        alpha = 0.3,  colour = "#68B69F",
        shape = 16, size = 2,
        raster.dpi = getOption("ggrastr.default.dpi", 300)
    ) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "black") +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "black") +
    xlab("apaQTL effect size") + ylab("eQTL effect size") +
    geom_smooth(
        method = lm, formula = y ~ x,
        colour = "red", se = FALSE, lwd = 1.25
    ) +
    annotate(
        "text", x = -1.3, y = 1.3,
        label =  as.character(as.expression(eq)),
        parse = TRUE, size = 5
    )
p1 <- p0 + mytheme
ggsave("output/top.apaQTL.in.eQTL.slope.pdf", p1, width = 3, height = 3)


## Pvalue heatmap
library(psych)
library(ComplexHeatmap)
rm(list=ls())

# top sQTL in eQTL
topsQTL <- fread("input/262sample_sqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, ":")[[5]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, -abs(slope), abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, variant_id, slope_sqtl = slope)
eQTL <- tidyfst::parse_fst("input/262sample_eqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topsQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_eqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_eqtl) %>%
    arrange(pair, -abs(slope_eqtl)) %>%
    distinct(pair, .keep_all = TRUE)
topsQTLineQTL <- topsQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_sqtl) %>%
    inner_join(eQTL, by = "pair") %>%
    rename(eqtl = slope_eqtl, sqtl = slope_sqtl)

# top eQTL in sQTL
topeQTL <- fread("input/262sample_eqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, -abs(slope), abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, variant_id, slope_eqtl = slope)
sQTL <- tidyfst::parse_fst("input/262sample_sqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topeQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, ":")[[5]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_sqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_sqtl) %>%
    arrange(pair, -abs(slope_sqtl)) %>%
    distinct(pair, .keep_all = TRUE)
topeQTLinsQTL <- topeQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_eqtl) %>%
    inner_join(sQTL, by = "pair") %>%
    rename(eqtl = slope_eqtl, sqtl = slope_sqtl)

# top apaQTL in eQTL
topapaQTL <- fread("input/262sample_aqtl.cis_qtl.perm.top.list") %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\|")[[2]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    arrange(phenotype_id, -abs(slope), abs(tss_distance)) %>%
    distinct(phenotype_id, .keep_all=T) %>%
    select(phenotype_id, variant_id, slope_apaqtl = slope)
eQTL <- tidyfst::parse_fst("input/262sample_eqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topapaQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_eqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_eqtl) %>%
    arrange(pair, -abs(slope_eqtl)) %>%
    distinct(pair, .keep_all = TRUE)
topapaQTLineQTL <- topapaQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_apaqtl) %>%
    inner_join(eQTL, by = "pair") %>%
    rename(eqtl = slope_eqtl, apaqtl = slope_apaqtl)

# top eQTL in apaQTL
apaQTL <- tidyfst::parse_fst("input/262sample_aqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topeQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\|")[[2]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_apaqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_apaqtl) %>%
    arrange(pair, -abs(slope_apaqtl)) %>%
    distinct(pair, .keep_all = TRUE)
topeQTLinapaQTL <- topeQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_eqtl) %>%
    inner_join(apaQTL, by = "pair") %>%
    rename(eqtl = slope_eqtl, apaqtl = slope_apaqtl)

# top sQTL in apaQTL
apaQTL <- tidyfst::parse_fst("input/262sample_aqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topsQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\|")[[2]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_apaqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_apaqtl) %>%
    arrange(pair, -abs(slope_apaqtl)) %>%
    distinct(pair, .keep_all = TRUE)
topsQTLinapaQTL <- topsQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_sqtl) %>%
    inner_join(apaQTL, by = "pair") %>%
    rename(sqtl = slope_sqtl, apaqtl = slope_apaqtl)

# top apaQTL in sQTL
sQTL <- tidyfst::parse_fst("input/262sample_sqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::select_fst(phenotype_id, variant_id, tss_distance, slope) %>%
    filter(variant_id %in% topapaQTL$variant_id) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, ":")[[5]]) %>%
    mutate(phenotype_id = tstrsplit(phenotype_id, "\\.")[[1]]) %>%
    select(phenotype_id, variant_id, slope_sqtl = slope) %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_sqtl) %>%
    arrange(pair, -abs(slope_sqtl)) %>%
    distinct(pair, .keep_all = TRUE)
topapaQTLinsQTL <- topapaQTL %>%
    mutate(pair = paste0(phenotype_id, "_", variant_id)) %>%
    select(pair, slope_apaqtl) %>%
    inner_join(sQTL, by = "pair") %>%
    rename(apaqtl = slope_apaqtl, sqtl = slope_sqtl)


## Save temp data
save(topeQTLinsQTL, topeQTLinapaQTL, topsQTLineQTL, topsQTLinapaQTL, topapaQTLineQTL, topapaQTLinsQTL, file = "temp/QTL_Pvalue.RData")


## Correlation test
dat <- matrix(
    c(
        NA,
        corr.test(topeQTLinsQTL$eqtl, topeQTLinsQTL$sqtl)$r,
        corr.test(topeQTLinapaQTL$eqtl, topeQTLinapaQTL$apaqtl)$r,
        corr.test(topsQTLineQTL$sqtl, topsQTLineQTL$eqtl)$r,
        NA,
        corr.test(topsQTLinapaQTL$sqtl, topsQTLinapaQTL$apaqtl)$r,
        corr.test(topapaQTLineQTL$apaqtl, topapaQTLineQTL$eqtl)$r,
        corr.test(topapaQTLinsQTL$apaqtl, topapaQTLinsQTL$sqtl)$r,
        NA
    ),
    nrow = 3, ncol = 3
)
colnames(dat) <- c("Top eQTL", "Top sQTL", "Top apaQTL")
rownames(dat) <- c("in eQTL", "in sQTL", "in apaQTL")


## Generate heatmap
pheatmap::pheatmap(
    mat = dat,
    scale = "none",
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    fontsize = 10,
    na_col = "grey",
    border_color = "white",
    angle_col = 0,
    color = colorRampPalette(c("#68B69F", "white", "#F08E59"))(10),
    breaks = seq(-1, 1, 0.2),
    filename = "output/Pvalue.Heatmap.pdf",
    width = 6, height = 5
)
