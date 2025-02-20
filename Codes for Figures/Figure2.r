##########################################################
# Code information
# Title: Code of Figure 2
# Author: Hu Beiping
# Email: hubeiping@njmu.edu.cn
##########################################################


R
options(stringsAsFactors=F)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrastr)
library(cowplot)
rm(list=ls())


##------------------------------Fig.2A

## Load the genotype file to customize the axis scale
genome <- fread("input/geno.bim") %>%
    select(SNP = V2, CHR = V1, BP = V4)

chrLen <- genome %>%
    group_by(CHR) %>%
    summarise(chrLen = max(BP))

chrPos <- chrLen %>%
    mutate(total = cumsum(as.numeric(chrLen)) - chrLen) %>%
    select(-chrLen)

snpPos <- chrPos %>%
    left_join(genome, ., by="CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + total)

xAxis <-  snpPos %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)


## Load the significant part of the three QTLs
# eQTL
eqtl <- fread("input/262sample_eqtl.cis_qtl.perm.list") %>%
    mutate(CHR = as.numeric(tstrsplit(variant_id, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = pval_nominal) %>%
    arrange(CHR, BP, P) %>%
    distinct(.keep_all = TRUE)

# sQTL
sqtl <- fread("input/262sample_sqtl.cis_qtl.perm.list") %>%
    mutate(CHR = as.numeric(tstrsplit(variant_id, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = pval_nominal) %>%
    arrange(CHR, BP, P) %>%
    distinct(.keep_all = TRUE)

# apaQTL
aqtl <- fread("input/262sample_aqtl.cis_qtl.perm.list") %>%
    mutate(CHR = as.numeric(tstrsplit(variant_id, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = pval_nominal) %>%
    arrange(CHR, BP, P) %>%
    distinct(.keep_all = TRUE)


## Load the previous GC GWAS summary data
gwas <- fread("input/GC.summary.txt") %>%
    filter((HetPVal > 0.0001) & (HetISq < 75) & (N > max(N)*2/3)) %>%
    filter(`P-value` < 5e-4) %>%
    select(MarkerName, `P-value`) %>%
    mutate(CHR = as.numeric(tstrsplit(MarkerName, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(MarkerName, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = `P-value`) %>%
    arrange(CHR, BP) %>%
    distinct(SNP, .keep_all = TRUE)


## Classify the SNPs according to whether they have significant QTLs
qtl <- eqtl %>%
    rbind(sqtl) %>%
    rbind(aqtl) %>%
    distinct(SNP) %>%
    mutate(eSNP = ifelse(SNP %in% eqtl$SNP, 1, 0)) %>%
    mutate(sSNP = ifelse(SNP %in% sqtl$SNP, 1, 0)) %>%
    mutate(aSNP = ifelse(SNP %in% aqtl$SNP, 1, 0)) %>%
    mutate(class = ifelse(
        (eSNP == 1) & (sSNP == 1) & (aSNP == 1),
        "eQTL & sQTL & apaQTL",
        ifelse(
            (eSNP == 1) & (sSNP == 1) & (aSNP == 0),
            "eQTL & sQTL",
            ifelse(
                (eSNP == 1) & (sSNP == 0) & (aSNP == 1),
                "eQTL & apaQTL",
                ifelse(
                    (eSNP == 0) & (sSNP == 1) & (aSNP == 1),
                    "sQTL & apaQTL",
                    ifelse(
                        (eSNP == 1) & (sSNP == 0) & (aSNP == 0),
                        "eQTL",
                        ifelse(
                            (eSNP == 0) & (sSNP == 1) & (aSNP == 0),
                            "sQTL",
                            "apaQTL"
                        )
                    )
                )
            )
        )
    )) %>%
    select(SNP, class)


## Merge the above data and save it as a temporary file
gwasPos <- chrPos %>%
    left_join(gwas, ., by = "CHR") %>%
    left_join(qtl, by = "SNP") %>%
    replace_na(list(class = "No QTLorGWAS")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + total) %>%
    mutate(class = factor(class, levels = c("No QTLorGWAS", "eQTL", "sQTL", "apaQTL", "eQTL & sQTL", "eQTL & apaQTL", "sQTL & apaQTL", "eQTL & sQTL & apaQTL"))) %>%
    filter(!(CHR == 2 & P < 5e-8))
fwrite(gwasPos, "temp/Figure2A_data.txt", col=T, row=F, quo=F, sep="\t")

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
    legend.position = c(0.85, 0.8),
    legend.title = element_blank()
)

p0 <- ggplot(subset(gwasPos, class = "No QTLorGWAS"), aes(x = BPcum, y = -log10(P))) +
    geom_hline(yintercept = -log10(5e-8), linewidth = 0.5, linetype = "dashed") +
    geom_point_rast(color = "grey", alpha = 1, size = 1) +
    geom_point_rast(
        data = subset(gwasPos, class != "No QTLorGWAS"),
        aes(x = BPcum, y = -log10(P), color = class),
        alpha = 1, size = 1.5
    ) +
    scale_color_manual(values = c("#145390", "#F08E59", "#68B69F", "#826BA2", "#C0C05A", "#4E8872", "#BD514A")) +
    scale_x_continuous(label = xAxis$CHR, breaks = xAxis$center) +
    scale_y_continuous(limits = c(0, 45), expand = c(0, 0)) +
    xlab("Chromosome") + ylab("P of GC GWAS (-log10)")
p1 <- p0 + mytheme
ggsave("output/Figure2A.pdf", p1, width = 10, height = 4)



##------------------------------Fig.2B

## Load the above data and save it as a temporary file
library(UpSetR)
upsetData <- gwas %>%
    mutate(eQTL = ifelse(SNP %in% eqtl$SNP, 1, 0)) %>%
    mutate(sQTL = ifelse(SNP %in% sqtl$SNP, 1, 0)) %>%
    mutate(apaQTL = ifelse(SNP %in% aqtl$SNP, 1, 0)) %>%
    filter(!(eQTL==0 & sQTL==0 & apaQTL==0))

upsetPlot <- upset(
    upsetData,
    sets = c("eQTL", "sQTL", "apaQTL"),
    order.by = c("freq", "degree"),
    main.bar.color = "#F08E59",
    matrix.color = "#145390",
    shade.color = "#E7E4E7",
    sets.bar.color = "#68B69F",
    text.scale = 2
)

pdf("output/Figure2B.pdf", width = 8, height = 6)
upsetPlot
dev.off()



##------------------------------Fig.2C

qtl <- eqtl %>%
    rbind(sqtl) %>%
    rbind(aqtl) %>%
    distinct(SNP) %>%
    mutate(eSNP = ifelse(SNP %in% eqtl$SNP, 1, 0)) %>%
    mutate(sSNP = ifelse(SNP %in% sqtl$SNP, 1, 0)) %>%
    mutate(aSNP = ifelse(SNP %in% aqtl$SNP, 1, 0))

gwas <- fread("input/GC.summary.txt") %>%
    filter((HetPVal > 0.0001) & (HetISq < 75) & (N > max(N)*2/3)) %>%
    filter(`P-value` < 5e-2) %>%
    select(MarkerName, `P-value`) %>%
    mutate(CHR = as.numeric(tstrsplit(MarkerName, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(MarkerName, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = `P-value`) %>%
    arrange(CHR, BP) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    left_join(qtl, by = "SNP") %>%
    replace_na(list(eSNP = 0, sSNP = 0, aSNP = 0))

dat <- data.frame(
    QTL = c(rep("eQTL", 2), rep("sQTL", 2), rep("apaQTL", 2)),
    N = c(
        nrow(gwas[gwas$eSNP==1 & gwas$P<5e-2,]), length(table(eqtl$SNP)) - nrow(gwas[gwas$eSNP==1 & gwas$P<5e-2,]),
        nrow(gwas[gwas$sSNP==1 & gwas$P<5e-2,]), length(table(sqtl$SNP)) - nrow(gwas[gwas$sSNP==1 & gwas$P<5e-2,]),
        nrow(gwas[gwas$aSNP==1 & gwas$P<5e-2,]), length(table(aqtl$SNP)) - nrow(gwas[gwas$aSNP==1 & gwas$P<5e-2,])
    )
)

test <- matrix(dat$N, ncol=3)
rownames(test) <- c("sig","nosig")
colnames(test) <- c("eQTL", "sQTL", "apaQTL")
test <- test %>% as.data.frame() %>% tibble::rownames_to_column(var = "P")
dat <- reshape2::melt(test, id.vars="P")

barTheme <- theme_classic() +
    theme(
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)), 
    axis.title.x = element_text(size = 10, color = "black"), 
    axis.title.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm"),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5),
    axis.line = element_line(colour = "black", size = 0.5),
    legend.position = "none"
)

p <- ggplot(dat, aes(variable, value), position = "stack") +
    scale_x_discrete(limits = c("eQTL", "sQTL", "apaQTL")) +
    geom_bar(
        aes(fill = P), stat = "identity", color = "white",
        size = 0.4, position = "stack", width = 0.6, data = dat
    ) +
    scale_fill_manual(values = c("#A1BAD3", "#4375A6")) +
    ylab("Number of significant QTLs") +
    barTheme
ggsave("output/Figure2C.pdf", p, width = 4, height = 6)



##------------------------------Fig.2D

## Calculate the enrichment degree using garfield
vim calculate.sig.garfield.ld.sh

#####################################
cd /Public/hbp/software/garfield-v2

INPUTNAME=GC.sigQTL.ld
DATADIR=/Public/hbp/3QTL/garfield

PRUNETAGSDIR=$DATADIR/Tag/r01
CLUMPTAGSDIR=$DATADIR/Tag/r08
MAFTSSDDIR=$DATADIR/maftssd
PVALDIR=$DATADIR/GWAS
ANNOTDIR=$DATADIR/sigData
OUTDIR=$DATADIR/Result/$INPUTNAME
mkdir -p $OUTDIR

ANNOTLINKFILE=$ANNOTDIR/link.file.txt
PTHRESH=5e-2,5e-3,5e-4,5e-5,5e-6,5e-7,5e-8
BINNING=m5,n5,t5
CONDITION=0
SUBSET=0

F1=$OUTDIR/garfield.prep.$INPUTNAME.out
F0=$OUTDIR/garfield.Meff.$INPUTNAME.out

echo 'Prune and Clump'
echo -n > $F1
for CHR in `seq 1 22`
do
	echo 'CHR'$CHR
	./garfield-prep-chr -ptags $PRUNETAGSDIR/chr$CHR -ctags $CLUMPTAGSDIR/chr$CHR -maftss $MAFTSSDDIR/chr$CHR -pval $PVALDIR/chr$CHR -ann $ANNOTDIR/chr$CHR -excl -1 -chr $CHR -o $F1 || { echo 'Failure!'; } 
done

echo 'Calculate effective number of annotations'
Rscript garfield-Meff-Padj.R -i $F1 -o $F0
NEA=$(head -1 $F0 |awk '{print $2}')
Padj=$(tail -1 $F0 |awk '{print $2}')

echo 'Calculate Enrichment and Significance'
F2=$OUTDIR/garfield.test.$INPUTNAME.out
Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -s $SUBSET -c $CONDITION

echo 'GARFIELD Analysis Complete!'
#####################################

chmod u+x calculate.sig.garfield.ld.sh
./calculate.sig.garfield.ld.sh


## Generate plot
R
options(stringsAsFactors=F)
library(data.table)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(forestplot)
rm(list=ls())

dat <- fread("temp/garfield.test.GC.sigQTL.out", h=T) %>%
    mutate(L_95CI = round(exp(CI95_lower), 2)) %>%
    mutate(U_95CI = round(exp(CI95_upper), 2)) %>%
    mutate(OR = round(OR, 2)) %>%
    filter(PThresh %in% c(0.05, 0.0005, 0.00000005)) %>%
    arrange(-PThresh) %>%
    select(Annotation, PThresh, OR, L_95CI, U_95CI)

subgps <- c(3:6, 8:11, 13:16)
dat$Variable[subgps] <- paste0("    ", dat$Variable[subgps]) 

attach(dat)
labeltext <- as.matrix(dat[,1:2])

fn_box <- local({
  i = 0
  function(..., clr.line, clr.marker) {
    i <<- i + 1
    if(i %% 4 == 1) {
        ## First group
        fpDrawCircleCI(..., clr.line = "#000000", clr.marker = "#145390")
    }
    else if(i %% 4 == 2) {
        ## Second group
        fpDrawCircleCI(..., clr.line = "#000000", clr.marker = "#F08E59")
    }
    else if(i %% 4 == 3) {
        ## Third group
        fpDrawCircleCI(..., clr.line = "#000000", clr.marker = "#68B69F")
    }
    else(
        ## Fourth group
        fpDrawCircleCI(..., clr.line = "#000000", clr.marker = "#BD514A")
    )
  }
})

styles <- fpShapesGp(
  lines = list(
    gpar(col = "white"),
    gpar(col = "white"),
    gpar(col = "#145390"),
    gpar(col = "#F08E59"),
    gpar(col = "#68B69F"),
    gpar(col = "#BD514A"),
    gpar(col = "white"),
    gpar(col = "#145390"),
    gpar(col = "#F08E59"),
    gpar(col = "#68B69F"),
    gpar(col = "#BD514A"),
    gpar(col = "white"),
    gpar(col = "#145390"),
    gpar(col = "#F08E59"),
    gpar(col = "#68B69F"),
    gpar(col = "#BD514A")
  )
)

pdf("output/Figure2D.pdf", onefile = FALSE, width = 10, height = 6)
forestplot(
    labeltext,        
    mean = OR, lower = LowerCI, upper = UpperCI,
    zero = 1, lwd.zero = 2, lwd.ci = 3,
    align = "l", xlog = TRUE,
    is.summary = c(T,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F),
    boxsize = 0.2, fn.ci_norm = fn_box, shapes_gp = styles
)
dev.off()



##------------------------------Supplementary Fig.8C

## Load data
rm(list=ls())

# eQTL
eQTL <- fread("input/262sample_eqtl.cis_qtl.perm.list") %>%
    mutate(CHR = as.numeric(tstrsplit(variant_id, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = pval_nominal) %>%
    arrange(CHR, BP, P) %>%
    distinct(.keep_all = TRUE)

# sQTL
sQTL <- fread("input/262sample_sqtl.cis_qtl.perm.list") %>%
    mutate(CHR = as.numeric(tstrsplit(variant_id, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = pval_nominal) %>%
    arrange(CHR, BP, P) %>%
    distinct(.keep_all = TRUE)

# apaQTL
apaQTL <- fread("input/262sample_aqtl.cis_qtl.perm.list") %>%
    mutate(CHR = as.numeric(tstrsplit(variant_id, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(variant_id, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = pval_nominal) %>%
    arrange(CHR, BP, P) %>%
    distinct(.keep_all = TRUE)

# GWAS
gwas <- fread("input/GC.summary.txt") %>%
    filter((HetPVal > 0.0001) & (HetISq < 75) & (N > max(N)*2/3)) %>%
    select(MarkerName, `P-value`) %>%
    mutate(CHR = as.numeric(tstrsplit(MarkerName, ":")[[1]])) %>%
    mutate(BP = as.numeric(tstrsplit(MarkerName, ":")[[2]])) %>%
    mutate(SNP = paste0(CHR, ":", BP)) %>%
    select(SNP, CHR, BP, P = `P-value`) %>%
    arrange(CHR, BP) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    mutate(eSNP = ifelse(SNP %in% eQTL$SNP, 1, 0)) %>%
    mutate(sSNP = ifelse(SNP %in% sQTL$SNP, 1, 0)) %>%
    mutate(aSNP = ifelse(SNP %in% apaQTL$SNP, 1, 0))

eSNP <- gwas %>%
    filter(eSNP == 1) %>%
    mutate(model = "eQTL") %>%
    select(model, value = P)

sSNP <- gwas %>%
    filter(sSNP == 1) %>%
    mutate(model = "sQTL") %>%
    select(model, value = P)

aSNP <- gwas %>%
    filter(aSNP == 1) %>%
    mutate(model = "apaQTL") %>%
    select(model, value = P)

dat <- gwas %>%
    mutate(model = "Genome-wide") %>%
    select(model, value = P) %>%
    rbind(eSNP) %>%
    rbind(sSNP) %>%
    rbind(aSNP) %>%
    drop_na(value) %>%
    mutate(model = factor(model, levels = rev(c("eQTL", "sQTL", "apaQTL", "Genome-wide")))) %>%
    arrange(model, value)

n <- table(dat$model)
dat$expected.p <- unlist(lapply(n, ppoints))
df <- data.frame(
    observed = -log10(dat$value),
    expected = -log10(dat$expected.p),
    model = dat$model
)

log10Pe <- expression(Expected~~-log[10](P))
log10Po <- expression(Observed~~-log[10](P))

save(df, log10Pe, log10Po, file = "temp/Figure2E.RData")


## Generate Q-Q plot
color_set <- rev(c("#145390", "#F08E59", "#68B69F", "grey"))
p <- ggplot(df) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.8, color = "red") +
    geom_point_rast(aes(expected, observed, color = model), shape = 16, alpha = 0.8, size = 1) +
    scale_discrete_manual(
        aesthetics = "color", values = color_set,
        guide = guide_legend(reverse = T, byrow = T,nrow = 4)
    ) +
    labs(x = log10Pe, y = log10Po, title = "") +
    theme(
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
        legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_blank()
    )
ggsave("output/Figure2E.pdf", p, width = 4, height = 6)