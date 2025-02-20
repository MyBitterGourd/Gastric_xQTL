##########################################################
# Code information
# Title: Code of Figure 4
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


## ------------------------------ Fig.4A top

## Load eQTL, sQTL and apaQTL results of PSCA
eQTLdat <- arrow::read_parquet("input/262sample_eqtl.cis_qtl_pairs.8.parquet") %>%
    filter(grepl("ENSG00000167653", phenotype_id)) %>%
    mutate(position = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    select(position, eQTLpval = pval_nominal)

sQTLdat <- arrow::read_parquet("input/262sample_sqtl.cis_qtl_pairs.8.parquet") %>%
    filter(phenotype_id == "8:143762852:143763339:clu_72986:ENSG00000167653.4_4") %>% 
    mutate(position = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    select(position, sQTLpval = pval_nominal)

apaQTLdat <- arrow::read_parquet("input/262sample_aqtl.cis_qtl_pairs.8.parquet") %>%
    filter(phenotype_id == "ENST00000513264.1_3|ENSG00000167653.4_4|chr8|+") %>%
    mutate(position = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    select(position, apaQTLpval = pval_nominal)

QTLdat <- eQTLdat %>%
    inner_join(sQTLdat, by = "position") %>%
    inner_join(apaQTLdat, by = "position") %>%
    mutate(
        class = case_when(
            position == 143761931 ~ "eQTL",
            position == 143763043 ~ "sQTL",
            position == 143762932 ~ "apaQTL",
            .default = "none"
        )
    )


## Gene region plot
eQTLplot <- ggplot(data = subset(QTLdat, class!="eQTL"), aes(x = position/1e6, y = -log10(eQTLpval))) +
    geom_point(
        color = "grey", pch = 16,
        size = 1, alpha = 0.8
    ) +
    geom_point(
        data = subset(QTLdat, class=="eQTL"),
        color = "#145390", pch = 16,
        size = 1.5, alpha = 1
    ) +
    ylab(bquote(~atop(italic(PSCA)~"eQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    geom_vline(xintercept = 143751600/1e6, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 143777000/1e6, color = "black", linewidth = 0.5, linetype = "dashed") +
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
        legend.position = c(0.15, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", linetype = "solid", linewidth = 0.5, colour = "black")
    )

sQTLplot <- ggplot(data = subset(QTLdat, class!="sQTL"), aes(x = position/1e6, y = -log10(sQTLpval))) +
    geom_point(
        color = "grey", pch = 16,
        size = 1, alpha = 0.8
    ) +
    geom_point(
        data = subset(QTLdat, class=="sQTL"),
        color = "#F08E59", pch = 16,
        size = 1.5, alpha = 1
    ) +
    ylab(bquote(~atop(italic(PSCA)~"sQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    geom_vline(xintercept = 143751600/1e6, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 143777000/1e6, color = "black", linewidth = 0.5, linetype = "dashed") +
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
        legend.position = c(0.15, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", linetype = "solid", linewidth = 0.5, colour = "black")
    )

apaQTLplot <- ggplot(data = subset(QTLdat, class!="apaQTL"), aes(x = position/1e6, y = -log10(apaQTLpval))) +
    geom_point(
        color = "grey", pch = 16,
        size = 1, alpha = 0.8
    ) +
    geom_point(
        data = subset(QTLdat, class=="apaQTL"),
        color = "#68B69F", pch = 16,
        size = 1.5, alpha = 1
    ) +
    xlab("Chr8 physical position (Mb)") +
    ylab(bquote(~atop(italic(PSCA)~"apaQTL", "-log"[10]*"("*italic("P")*")"))) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    geom_vline(xintercept = 143751600/1e6, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 143777000/1e6, color = "black", linewidth = 0.5, linetype = "dashed") +
    theme(
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


## Merge three plots vertically
p <- eQTLplot / sQTLplot / apaQTLplot
ggsave("output/Figure4Atop.pdf", p, width = 6, height = 5)



## ------------------------------ Fig.4A bottom

## Extract gene structure
library(karyoploteR)
library(BiocFileCache)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

PSCA.region <- toGRanges("chr8:143756000-143768000")
kp <- plotKaryotype(zoom = PSCA.region)
genes.data <- makeGenesDataFromTxDb(
    TxDb.Hsapiens.UCSC.hg19.knownGene,
    karyoplot = kp,
    plot.transcripts = TRUE, 
    plot.transcripts.structure = TRUE
)
genes.data <- addGeneNames(genes.data)
genes.data$transcripts$`8000` <- genes.data$transcripts$`8000`[genes.data$transcripts$`8000`$tx_id == "32868"]
genes.data <- mergeTranscripts(genes.data)


## Download chromatin state data
bfc <- BiocFileCache(cache = "temp", ask = FALSE)
K562.hmm.file <- bfcrpath(bfc, "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz")
K562.hmm <- regioneR::toGRanges(K562.hmm.file)
chromHMM <- as.vector(unique(K562.hmm$V4))
chromHMM_sort <- stringr::str_sort(chromHMM, numeric = TRUE)
K562.hmm$V4 <- factor(K562.hmm$V4, levels=chromHMM_sort)
colors <- c("#FF0000","#CD5C5C","#800080","#FFA500","#FFA500","#FFFF00","#FFFF00","#00BFFF","#006400","#006400","#90EE90","#808080","#D3D3D3","#D3D3D3","#D3D3D3")
K562.hmm$V6 <- as.character(factor(K562.hmm$V4, labels = colors))
K562.hmm <- toGRanges(K562.hmm)

kp <- plotKaryotype(zoom = PSCA.region, cex = 1)
kpPlotGenes(kp, data = genes.data, r0 = 0, r1 = 0.15, gene.name.cex = 1)
kpPlotRegions(kp, K562.hmm, col = K562.hmm$V6, r0 = 0.22, r1 = 0.3, avoid.overlapping = F)
kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0 = 0.22, r1 = 0.3, cex = 1)


## Load annotation data from ABC model
promoter <- makeGRangesFromDataFrame(data.frame(chr="chr8", start=143761563, end=143762123, score=1))
enhancer <- makeGRangesFromDataFrame(
    data.frame(
        chr = rep("chr8", 4),
        start = c(143756350, 143757428, 143758285, 143764542),
        end = c(143756850, 143757928, 143759158, 143766727),
        score = rep(1, 4)
    )
)
TSS <- makeGRangesFromDataFrame(
    data.frame(
        chr = rep("chr8", 4),
        start = rep(143761873, 4),
        end = rep(143761873, 4),
        score = rep(1, 4)
    )
)


## High LD SNPs
ldSNPs <- fread("input/topQTL.flank1Mb.1kgEAS.rs2920283.ld") %>%
    mutate(chr = paste0("chr", CHR_B)) %>%
    mutate(start = BP_B - 1) %>%
    dplyr::select(chr, start, end = BP_B, rsID = SNP_B) %>%
    filter(start>=143756000 & start<=143768000)
ldSNPs.gr <- toGRanges(ldSNPs)


## Histone Marks
annotation.marks <- data.frame(
        path = c(
            "input/1391712-N.last.shift.sort.bw",
            "input/1391712-N_H3K4me1.last.bw",
            "input/1391712-N_H3K4me3.last.bw",
            "input/1391712-N_H3K27ac_abcam.last.bw"
        ),
        signal = c("ATAC-seq (NMU)", "H3K4me1 (NMU)", "H3K4me3 (NMU)", "H3K27ac (NMU)")
    )


## Generate track plot
# Set plot parameters
pp <- getDefaultPlotParams(plot.type = 1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

set.seed(22)
pdf("output/Figure4Abottom.pdf", width = 6, height = 6)

    # Gene structure
    kp <- plotKaryotype(zoom = PSCA.region, cex = 0.5, plot.params = pp)
    kpAddBaseNumbers(
        kp, tick.dist = 10000, minor.tick.dist = 2000,
        add.units = TRUE, cex = 0.5, digits = 6
    )
    kpPlotGenes(
        kp, data = genes.data,
        r0 = 0, r1 = 0.03,
        gene.name.cex = 0.5,
        gene.name.position = "right"
    )

    # HMM model
    kpPlotRegions(kp, K562.hmm, col = K562.hmm$V6, r0 = 0.04, r1 = 0.08, avoid.overlapping = F)
    kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0 = 0.04, r1 = 0.08, cex = 0.5)

    # ABC model
    kpPlotRegions(kp, promoter, r0=0.10, r1=0.11, col="#FF8D92")
    kpPlotRegions(kp, enhancer, r0=0.10, r1=0.11, col="#8D9AFF")
    kpPlotRegions(kp, TSS, r0=0.10, r1=0.11, col="#8D9AFF")
    kpPlotLinks(kp, data=enhancer, data2=TSS, col="#fac7ffaa", r0=0.11, r1 = 0.15)
    kpAddLabels(kp, labels = "ABC Model(NMU)", r0 = 0.11, r1 = 0.15, cex = 0.5)

    # Histone marks
    total.tracks <- nrow(annotation.marks) + 1
    out.at <- autotrack(1:nrow(annotation.marks), total.tracks, margin = 0.1, r0 = 0.16)
    for (i in seq_len(nrow(annotation.marks))) {
        bigwig.file <- annotation.marks$path[i]
        at <- autotrack(i, nrow(annotation.marks), r0 = out.at$r0, r1 = out.at$r1, margin = 0.1)
        kp <- kpPlotBigWig(
            kp, data = bigwig.file, ymax = "visible.region",
            r0 = at$r0, r1 = at$r1,
            col = "cadetblue2", border = NA
        )
        computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
        kpAxis(kp, ymin = 0, ymax = computed.ymax, numticks = 2, r0 = at$r0, r1 = at$r1, cex = 0.5)
        kpAddLabels(
            kp, labels = annotation.marks$signal[i],
            r0 = at$r0, r1 = at$r1,
            cex = 0.5, label.margin = 0.035
        )
    }
    at <- autotrack(total.tracks, total.tracks, margin = 0.1, r0 = 0.16)
    kp <- kpPoints(
        kp, data = ldSNPs.gr,
        r0 = at$r0, r1 = at$r1,
        y = sample(seq(0.1, 1, 0.1), 47, replace=T)
    )
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp, ymin = 0, ymax = computed.ymax, numticks = 2, r0 = at$r0, r1 = at$r1, cex = 0.5)
    kpAddLabels(kp, labels = "High LD SNPs\n(R2>0.6)", r0 = at$r0, r1 = at$r1, cex = 0.5)
dev.off()



## ------------------------------ Fig.4B

## Set alleles
REF <- "C"
ALT <- "T"
WT <- paste0(REF, REF)
HET <- paste0(REF, ALT)
HOM <- paste0(ALT, ALT)


## Load expression data
exp <- fread("input/262sample_Stomach.expression.bed.gz") %>%
    filter(gene_id == "ENSG00000167653.4") %>%
    dplyr::select(-c(1:3)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(exp = 2)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr8.fst") %>%
    tidyfst::filter_fst(ID == "8:143761931:G:A") %>%
    dplyr::select(-c(1,2,4:9)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(snp = 2) %>%
    mutate(
        snp = case_when(
            snp == "0/0" ~ WT,
            snp %in% c("0/1", "1/0") ~ HET,
            snp == "1/1" ~ HOM,
            .default = NA
        )
    )


## Load covariate data
covar <- fread("input/262sample_Stomach.eQTL.combined_covariates.txt") %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID")


## Merge data and calculate residuals
final <- exp %>%
    inner_join(covar, by = "ID")
geno <- geno[match(final$ID, geno$ID),]
final <- final %>%
    dplyr::select(-1) %>%
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
    tidyfst::filter_fst(phenotype_id == "ENSG00000167653.4") %>%
    tidyfst::filter_fst(variant_id == "8:143761931:G:A")
Slope <- round(temp$slope, 3)
Pvalue <- format(temp$pval_nominal, digits = 3)


## Generate boxplot
eQTLboxp0 <- ggplot(data = geno, aes(x = snp, y = residuals)) +
    geom_boxplot(
        aes(color = snp), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = snp), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#7098BD", "#3F76A7", "#145390"), name = "") +
    xlab("rs2294008 genotype") +
    ylab("Normalized expression") +
    ggtitle(
        bquote(italic(PSCA)),
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
ggsave("output/Figure4B.pdf", eQTLboxp1, width = 3, height = 3)



## ------------------------------ Fig.4C

## Set intron ID
Chr <- tstrsplit("8:143762852:143763339:clu_72986:ENSG00000167653.4_4", ":")[[1]]
upper <- tstrsplit("8:143762852:143763339:clu_72986:ENSG00000167653.4_4", ":")[[2]]
lower <- tstrsplit("8:143762852:143763339:clu_72986:ENSG00000167653.4_4", ":")[[3]]
IntronID <- paste0("chr", Chr, ":", upper, "-", lower)


## Set alleles
REF <- tstrsplit("8:143763043:A:G", ":")[[3]]
ALT <- tstrsplit("8:143763043:A:G", ":")[[4]]
WT <- paste0(REF, REF)
HET <- paste0(REF, ALT)
HOM <- paste0(ALT, ALT)


## Load splicing data
splicing <- fread("input/262sample_Stomach.splicing.bed.gz") %>%
    filter(ID == "8:143762852:143763339:clu_72986:ENSG00000167653.4_4") %>%
    dplyr::select(-c(1:3)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(splicing = 2)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr8.fst") %>%
    tidyfst::filter_fst(ID == "8:143763043:A:G") %>%
    dplyr::select(-c(1,2,4:9)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(snp = 2) %>%
    mutate(
        snp = case_when(
            snp == "0/0" ~ WT,
            snp %in% c("0/1", "1/0") ~ HET,
            snp == "1/1" ~ HOM,
            .default = NA
        )
    )


## Load covariate data 
covar <- fread("input/262sample_Stomach.sQTL.combined_covariates.txt") %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID")


## Merge data and calculate residuals
final <- splicing %>%
    inner_join(covar, by = "ID")
geno <- geno[match(final$ID, geno$ID),]
final <- final %>%
    dplyr::select(-1) %>%
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
temp <- arrow::read_parquet("input/262sample_sqtl.cis_qtl_pairs.8.parquet") %>%
    filter(phenotype_id == "8:143762852:143763339:clu_72986:ENSG00000167653.4_4") %>%
    filter(variant_id == "8:143763043:A:G")
Slope <- round(temp$slope, 3)
Pvalue <- format(temp$pval_nominal, digits = 3)


## Generate boxplot
sQTLboxp0 <- ggplot(data = geno, aes(x = snp, y = residuals)) +
    geom_boxplot(
        aes(color = snp), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = snp), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#F7BC9B", "#F5A57A", "#F08E59"), name = "") +
    xlab("rs2920298 genotype") +
    ylab("Normalized intron usage") +
    ggtitle(
        bquote(italic(PSCA)~.(IntronID)),
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
ggsave("output/Figure4C.pdf", sQTLboxp1, width = 3, height = 3)



## ------------------------------ Fig.4D

## Load list of samples
list <- openxlsx::read.xlsx("input/262sample.info.xlsx") %>%
    dplyr::select(ID, bamID)


## Calculate percentage of intron usage
## 内含子矩阵计算百分比
psin <- fread("input/262sample.Stomach_perind.filtered.counts.txt.gz", data.table = F) %>%
    filter(grepl(":clu_72986_NA", chrom)) %>%
    mutate(chrom = gsub(":clu_72986_NA", "", chrom)) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "bamID") %>%
    mutate(`8:143762779:143763339` = tstrsplit(`8:143762779:143763339`, "/")[[1]]) %>%
    mutate(`8:143762852:143763339` = tstrsplit(`8:143762852:143763339`, "/")[[1]]) %>%
    mutate_at(2:3, as.numeric) %>%
    mutate(all = `8:143762779:143763339`+`8:143762852:143763339`) %>%
    mutate(`8:143762779:143763339` = `8:143762779:143763339`/all) %>%
    mutate(`8:143762852:143763339` = `8:143762852:143763339`/all) %>%
    dplyr::select(-4)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr8.fst") %>%
    tidyfst::filter_fst(grepl("8:143763043", ID)) %>%
    dplyr::select(-c(1,2,4:9)) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(snp = 2) %>%
    mutate(
        snp = case_when(
            snp == "0/0" ~ "AA",
            snp %in% c("0/1", "1/0", "1/1") ~ "AG + GG",
            .default = NA
        )
    )

psd <- list %>%
    inner_join(psin, by = "bamID") %>%
    inner_join(geno, by = "ID") %>%
    dplyr::select(5,3:4)

psda <- data.frame()
for (g in unique(psd$snp)) {
	sub_psd <- subset(psd, snp == g)
	sub_psda <- data.frame(
        rs2920298 = g, 
		chr8.143762779.143763339 = mean(sub_psd$`8:143762779:143763339`),
        chr8.143762852.143763339 = mean(sub_psd$`8:143762852:143763339`))
	psda <- rbind(psda, sub_psda)
}

tpsda <- t(psda)
write.table(tpsda, "temp/PSCA.subclu.meanExp.txt", col.names=F, quote=F, sep="\t")

psd <- read.delim("temp/PSCA.subclu.meanExp.txt")
colnames(psd) <- c("id", "AG+GG", "AA")

cluster <- data.frame(id = colnames(psda)[2:3]) %>%
    mutate(Start = tstrsplit(id, "\\.")[[2]]) %>%
    mutate(End = tstrsplit(id, "\\.")[[3]]) %>%
    mutate_at(2:3, as.numeric)

clu1 <- merge(cluster, psd[,c("id", "AG+GG")], by = "id")
clu2 <- merge(cluster, psd[,c("id", "AA")], by = "id")

colnames(clu1)[4] <- "MeanExp"
clu1$rs2920298 <- "AG+GG"
clu1$yorder <- 1
colnames(clu2)[4] <- "MeanExp"
clu2$rs2920298 <- "AA"
clu2$yorder <- 3
clu_use <- rbind(clu1, clu2)

clu_use$up <- ifelse(clu_use$id=="chr8.143762779.143763339", "No", "Yes")
clu_use$N <- ifelse(
    clu_use$rs2920298 == "AG+GG", 119,
    ifelse(
        clu_use$rs2920298 == "AA", 143, NA
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
    Start = rep(c(min(clu_use$Start)-(max(clu_use$Start)-min(clu_use$Start))*0.1, 143762779, 143763339), 2),
    End = rep(c(143762779, 143762852, max(clu_use$End)+(max(clu_use$End)-min(clu_use$Start))*0.1), 2),
    ymin = rep(unique(clu_use$ymin), each=3),
    ymax = rep(unique(clu_use$ymax), each=3)
)
genomeData <- data.frame(
    Start = rep(min(clu_use$Start)-(max(clu_use$Start)-min(clu_use$Start))*0.1, 2),
    End = rep(max(clu_use$End)+(max(clu_use$End)-min(clu_use$Start))*0.1, 2),
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
        data = subset(clu_use, id == "chr8.143762779.143763339"), 
		aes(x = Start, y = yorder, xend = End, yend = yorder, color = MeanExp, size = MeanExp),
		curvature = 0.4, alpha = 1, lineend = "butt"
    ) +
	geom_curve(
        data = subset(clu_use, id == "chr8.143762852.143763339"),
        aes(x = Start, y = yorder, xend = End, yend = yorder, color = MeanExp, size = MeanExp),
        curvature = 0.3, alpha = 1, lineend = "butt"
    ) +
    ylim(0, 6.5) +
    scale_size(range = c(0.5, 1.5)) +
    scale_color_gradient(low = "#F9D2BD", high = "#F08E59", name = "")

text_x <- min(clu_use$Start) - (max(clu_use$End)-min(clu_use$Start))*0.1
p1 <- p0 +
    geom_text(
        data = subset(clu_use, !duplicated(rs2920298)),
        aes(x = text_x-50, y = yorder, label = paste0(rs2920298, "\n(n=", N, ")")),
        size = 8
    ) +
    geom_text(
        data = subset(clu_use, id == "chr8.143762779.143763339"),
		aes(x = 0.5*(Start+End)-40, y = yorder-0.9, label = paste0(round(MeanExp*100,2), "%")),
        size = 8
    ) +		
    geom_text(
        data = subset(clu_use, id == "chr8.143762852.143763339"),
		aes(x = 0.5*(Start+End), y = yorder-0.3, label = paste0(round(MeanExp*100,2), "%")),
        size = 8
    ) +
    xlim(text_x-100, max(clu_use$End)+(max(clu_use$End)-min(clu_use$Start))*0.1)

# SNP
p2 <- p1 + 
    geom_segment(
	    aes(x = 143763043, y = 0.6, xend = 143763043, yend = 3.5),
        color = "#666666", size = 1, linetype = "dotted"
    ) +		
	geom_segment(
	    aes(x = 143763043, y = 3.5, xend = 143763043+80, yend = 3.8),
        color = "#666666", size = 1, linetype = "dotted"
    ) +
    geom_text(
        aes(x = 143763043+150, y = 3.9, label = "Associated SNP:\nrs2920298"),
        size = 6
    )

# junction
p3 <- p2 +
    geom_segment(
	    aes(x = 143762779, y = 3, xend = 143762779, yend = 3.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 143762779 + 50, y = 4.15, label = "chr8:143,762,779"),
        size = 6, angle = 45
    ) +
    geom_segment(
	    aes(x = 143762852, y = 3, xend = 143762852, yend = 3.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 143762852 + 50, y = 4.15, label = "chr8:143,762,852"),
        size = 6, angle = 45
    ) +
    geom_segment(
	    aes(x = 143763339, y = 3, xend = 143763339, yend = 3.5),
        color = "black", size = 0.5
    ) +
    geom_text(
        aes(x = 143763339 + 50, y = 4.15, label = "chr8:143,763,339"),
        size = 6, angle = 45
    )
ggsave("output/Figure4D.pdf", p3, width = 10, height = 6)



## ------------------------------ Fig.4D

## Set transcript ID
TranscriptID <- tstrsplit("ENST00000513264.1_3|ENSG00000167653.4_4|chr8|+", "\\|")[[1]]
TranscriptID <- tstrsplit(TranscriptID, "\\.")[[1]]


## Set alleles
REF <- tstrsplit("8:143762932:G:A", ":")[[3]]
ALT <- tstrsplit("8:143762932:G:A", ":")[[4]]
WT <- paste0(REF, REF)
HET <- paste0(REF, ALT)
HOM <- paste0(ALT, ALT)


## Load APA data
apa <- fread("input/262sample_Stomach.apa.bed.gz") %>%
    filter(Gene == "ENST00000513264.1_3|ENSG00000167653.4_4|chr8|+") %>%
    dplyr::select(-c(1:3)) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(apa = 2)


## Load genotype data
geno <- tidyfst::parse_fst("input/chr8.fst") %>%
    tidyfst::filter_fst(ID == "8:143762932:G:A") %>%
    dplyr::select(-c(1,2,4:9)) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    dplyr::rename(snp = 2) %>%
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
    dplyr::select(-1) %>%
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


## Load QTL results
temp <- tidyfst::parse_fst("input/262sample_aqtl.cis_qtl_pairs.all.chr.fst") %>%
    tidyfst::filter_fst(phenotype_id == "ENST00000513264.1_3|ENSG00000167653.4_4|chr8|+") %>%
    tidyfst::filter_fst(variant_id == "8:143762932:G:A")
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
    xlab("rs2976392 genotype") +
    ylab("Normalized PDUI") +
    ggtitle(
        bquote(italic(PSCA)~.(TranscriptID)),
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
ggsave("output/Figure4E.pdf", apaQTLboxp1, width = 3, height = 3)