##########################################################
# Code information
# Title: Code of Figure 3
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################

R
options(stringsAsFactors=F)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrastr)
library(patchwork)
rm(list=ls())


## ------------------------------ Fig.3A

library(circlize)

## Load significant results from three TWAS
TWASgene <- fread("input/twas.chr.sig.Fusion.results_adj", h=T) %>%
    filter(`TWAS.FDR` < 0.1) %>%
    select(Chrom, Start, End, gene_name)

spTWASgene <- fread("input/sptwas.chr.sig.Fusion.results_adj", h=T) %>%
    filter(`TWAS.FDR` < 0.1) %>%
    select(Chrom, Start, End, gene_name)

apaTWASgene <- fread("input/apatwas.chr.sig.Fusion.results_adj", h=T) %>%
    filter(`TWAS.FDR` < 0.1) %>%
    select(Chrom, Start, End, gene_name)

allTWASgene <- TWASgene %>%
    bind_rows(spTWASgene) %>%
    bind_rows(apaTWASgene) %>%
    distinct() %>%
    arrange(as.numeric(gsub("chr", "", Chrom)), Start) %>%
    select(chr = Chrom, start = Start, end = End, gene_name)

allTWASgeneValue <- allTWASgene %>%
    mutate(TWAS = ifelse(gene_name %in% TWASgene$gene_name, 3, 0)) %>%
    mutate(spTWAS = ifelse(gene_name %in% spTWASgene$gene_name, 2, 0)) %>%
    mutate(apaTWAS = ifelse(gene_name %in% apaTWASgene$gene_name, 1, 0)) %>%
    select(chr, start, end, TWAS, spTWAS, apaTWAS)

colSet = colorRamp2(c(0, 1, 2, 3), c("white", "#68B69F", "#F08E59", "#145390"))


## Generate circle plot
pdf("output/Figure3A.pdf", width = 6, height = 6)
    circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), plotType = NULL)
    circos.genomicLabels(allTWASgene, labels.column = 4, cex = 0.5, side = "outside", line_col = "black", connection_height = NULL)
    circos.genomicHeatmap(allTWASgeneValue, col = colSet, side = "outside", border = "black")
    circos.genomicIdeogram()
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        circos.genomicAxis(h = "bottom", direction = "inside")
    })
dev.off()



## ------------------------------ Fig.3B

## Load significant results from three TWAS
dat <- openxlsx::read.xlsx("input/xTWAS.sigGene.xlsx", sheet = 2) %>%
    tibble::column_to_rownames("Class") %>%
    replace(x = ., list = is.na(.), values = 0)


## Generate heatmap
plotB <- ComplexHeatmap::Heatmap(
    as.matrix(dat),                         
    col = c("white", "#F08E59"),
    show_heatmap_legend = FALSE,
    rect_gp = grid::gpar(col = "grey30", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 8, angle=45, fontface="italic"),
    column_names_side = "top",               
    row_names_gp = grid::gpar(fontsize = 8),
    row_names_side = "left",
    column_names_rot = 90
)
pdf("output/Figure3B.pdf", width = 8, height = 2)
plotB
dev.off()



## ------------------------------ Fig. 3C

## Load data
rm(list=ls())

genelist <- openxlsx::read.xlsx("input/GC.sigGene.xlsx") %>%
    rename(GeneID = `Gene.ID`, GeneName = `Gene.Name`) %>%
    select(GeneName)
model <- fread("input/Model.csv", sep=",") %>%
    filter(DepmapModelType == "STAD")
dependency <- fread("input/CRISPRGeneDependency.csv", sep=",") %>%
    filter(V1 %in% model$ModelID) %>%
    tibble::column_to_rownames("V1")
colNames <- tstrsplit(colnames(dependency), " " )[[1]]
colnames(dependency) <- colNames


## Extarct HGC-27 cell line dependency data
HGCdependency <- dependency %>%
    slice(which(rownames(dependency) == model[CellLineName=="HGC-27",]$ModelID)) %>%
    t() %>%
    as.data.frame() %>%
    rename(value = 1) %>%
    tibble::rownames_to_column(var = "GeneName") %>%
    full_join(genelist, by = "GeneName") %>%
    replace_na(list(value = 0))
HGCdependency <- HGCdependency %>%
    arrange(value) %>%
    mutate(rank = 1:nrow(HGCdependency)) %>%
    mutate(class = ifelse(value>=0.5 & GeneName%in%genelist$GeneName, "sig", "nonsig")) %>%
    mutate(class = factor(class, levels = c("nonsig", "sig")))


## Generate plot
p1 <- ggplot(HGCdependency, aes(x = rank, y = value, color = class)) +
    geom_line(col = "#BFBFBF") +
    geom_point(
        data = subset(HGCdependency, class=="sig"),
        alpha = 1, size = 1.5,
        color = "#BD514A"
    ) +
    ggrepel::geom_text_repel(
        data = subset(HGCdependency, class=="sig"),
        aes(x = rank, y = value, label = GeneName),
        size = 2, force = 0.1,
        color = "#BD514A"
    ) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    labs(x = "Gene sorted by score", y = "CRISPR dependency score\nin GC cell line") +
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
        legend.position = "none"
    )


## Calculate the results of significant genes in each cell line
dependency <- fread("input/CRISPRGeneDependency.csv", sep=",") %>%
    tibble::column_to_rownames("V1")
colNames <- tstrsplit(colnames(dependency), " \\(")[[1]]
colnames(dependency) <- colNames

model <- fread("input/Model.csv", sep=",")

results <- data.frame()
for (i in seq_len(nrow(dependency))) {
    print(paste0("Processing line ", i, "..."))
    cellLine <- model[ModelID==rownames(dependency)[i],]$CellLineName
    temp <- dependency %>%
        slice(i)
    temp <- temp %>%
        select(which(colnames(temp) %in% genelist$GeneName)) %>%
        t() %>%
        as.data.frame() %>%
        mutate_at(1, as.numeric) %>%
        drop_na() %>%
        pull(1)
    res <- c(cellLine, model[ModelID==rownames(dependency)[i],]$DepmapModelType, length(temp[temp>=0.5]))
    results <- rbind(results, res)
}
colnames(results) <- c("cellLine", "modelName", "sigCount")
results <- results %>%
    mutate(sigCount = as.integer(sigCount)) %>%
    mutate(isSTAD = ifelse(modelName=="STAD", 1, 0)) %>% # whether it belongs to GC cell line
    arrange(sigCount, isSTAD) %>%
    mutate(rank = 1:nrow(results))
fwrite(results, "temp/TWAS.genes.allcell.list", col=T, row=F, quo=F, sep="\t")


## Generate plot
p2 <- ggplot(results, aes(x = rank, y = sigCount, color = factor(isSTAD))) +
    geom_line(col = "#BFBFBF") +
    geom_point(
        data = subset(results, isSTAD==1),
        alpha = 1, size = 1.5,
        color = "#BD514A"
    ) +
    ylim(2, 13) +
    ggrepel::geom_text_repel(
        data = subset(results, cellLine=="HGC-27"),
        aes(x = rank, y = sigCount, label = cellLine),
        size = 2, force = 0.1,
        color = "#BD514A"
    ) +
    labs(x = "Cell line sorted by significant genes", y = "Number of significant\ngenes in cell line") +
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
        legend.position = "none"
    )
p <- p1 / p2
ggsave("output/Figure3C.pdf", p, width = 3.5, height = 5)



## ------------------------------ Fig. 3D

## Load data
TWASgene <- openxlsx::read.xlsx("input/xTWAS.sigGene.xlsx")
genelevels <- TWASgene$`Gene.Name`
bubbleplotData <- fread("input/TWASgene.GDSC2.AllCancers.cor.txt") %>%
    filter(TargetPathway == "WNT signaling") %>%
    mutate(GeneName = factor(GeneName, levels = genelevels))


## Generate plot
drugBubble <- ggplot(bubbleplotData, aes(x = GeneName ,y = DrugName)) + 
    geom_point(
        aes(size = -log10(Pbon), color = Correlation),
        shape = 20
    ) +
    scale_color_gradient2(
        name = "Spearman Correlation",
        low = "#145390", high = "#BD514A", 
        mid = "white", midpoint = 0,
        limit = c(-0.4,0.4),
        breaks = c(-0.4,-0.2,0.0,0.2,0.4),
    ) +
    scale_size_continuous(
        name = "-log10 FDR",
        limit = c(-0.001, 20),
        breaks = c(0, 10, 20),
        range = c(2, 6)
    ) +
    scale_x_discrete(position = "top") +
    theme(
        panel.grid.major = element_line(linewidth = 0.5, color = "#BFBFBF"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(0, "cm"),
        axis.ticks.length.y = unit(0, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_blank()
    )



## ------------------------------ Fig. 3E

## Load data
TWASgene <- openxlsx::read.xlsx("input/xTWAS.sigGene.xlsx") %>%
    rename(GeneID = `Gene.ID`, GeneName = `Gene.Name`)


## Load TIMER2 results
timer2 <- fread("input/infiltration_estimation_for_tcga.csv.gz", sep=",")
timer2 <- timer2 %>%
    select(c(1, which(grepl("_TIMER",colnames(timer2))))) %>%
    rename(sampleID = cell_type)


## Load TCGA STAD RNA-seq data
exp <- fread("input/TCGA.STAD.combined_RNAseq_TPM.txt") %>%
    filter(V1 %in% TWASgene$GeneID) %>%
    select(-c("gene_symbol", "type")) %>%
    rename(gene_symbol = V1)
colnames(exp) <- substr(colnames(exp), 1, 15)


## Set the order of genes and immune cells
genelevels <- TWASgene$GeneName
typelevels <- gsub("_TIMER", "", colnames(timer2)[2:ncol(timer2)])
results <- fread("input/TWASgene.TIMER2.STAD.cor.txt") %>%
    mutate(ImmuneCell = factor(ImmuneCell, levels = typelevels)) %>%
    mutate(GeneName = factor(GeneName, levels = genelevels))


## Generate plot
cellBubble <- ggplot(results, aes(x = GeneName ,y = ImmuneCell)) + 
    geom_point(
        aes(size = -log10(Pbon), color = Correlation),
        shape = 20
    ) +
    scale_color_gradient2(
        name = "Spearman Correlation",
        low = "#145390", high = "#BD514A", 
        mid = "white", midpoint = 0,
        limit = c(-0.4,0.8),
        breaks = c(-0.4,0.0,0.4,0.8),
    ) +
    scale_size_continuous(
        name = "-log10 FDR",
        limit = c(-0.001, 80),
        breaks = c(0, 20, 40, 60, 80),
        range = c(2, 6)
    ) +
    theme(
        panel.grid.major = element_line(linewidth = 0.5, color = "#BFBFBF"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        axis.ticks.length.y = unit(0, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_blank()
    )

p <- drugBubble / cellBubble + plot_layout(nrow = 2, heights = c(1.2, 1))
ggsave("output/Figure3DE.pdf", p, width = 8, height = 6)