# Not used.
# Based on https://f1000research.com/articles/5-1384/v1

library(Biobase)
library(oligoClasses)
library(knitr)
library(BiocStyle)
library(oligo)
library(geneplotter)
library(ggplot2)
library(dplyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(ArrayExpress)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(topGO)
library(genefilter)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pheatmap)
library(mvtnorm)
library(DAAG)
library(multcomp)
library(limma)
library(ReactomePA)
library(clusterProfiler)
library(openxlsx)
library(devtools)
library(biomaRt)
library(EnrichmentBrowser)

setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/analysis")
## Load data: 2
SDRF <- read.delim("/home/gtan/analysis/p1688-Ula/ImperialCollege/data/BAL_expt_design_fixed.txt")
SDRF$Array.Data.File <- SDRF$X
SDRF$X <- NULL
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

raw_data <- read.celfiles(filenames = file.path(celfiles, SDRF$Array.Data.File),
                          verbose = FALSE, phenoData = SDRF)

exp_raw <- log2(exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Time = pData(raw_data)$TIME,
                     Group = pData(raw_data)$GROUP,
                     SAMPLE_ID = pData(raw_data)$SAMPLE_ID,
                     GT = pData(raw_data)$GT,
                     labels = paste(pData(raw_data)$TIME, 
                                    pData(raw_data)$SAMPLE_ID, sep="_"))

p <- qplot(PC1, PC2, data = dataGG, color = Group,
           main = "PCA plot of the raw data", size = I(5),
           asp = 1.0, geom = "text",
           label = labels) + 
  scale_colour_brewer(palette = "Set2") + theme_bw()
ggsave(filename="PCA_rawData.pdf", p, device="pdf", width=12, height=10)

pdf("boxplot_rawData.pdf", width=10, height=7)
boxplot(raw_data, target = "core",
        main = "Boxplots of log2-intensities for the raw data")
dev.off()

palmieri_eset <- oligo::rma(raw_data, target="core")

exp_palmieri <- exprs(palmieri_eset)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Time = pData(palmieri_eset)$TIME,
                     Group = pData(palmieri_eset)$GROUP,
                     SAMPLE_ID = pData(palmieri_eset)$SAMPLE_ID,
                     GT = pData(palmieri_eset)$GT,
                     labels = paste(pData(palmieri_eset)$TIME, 
                                    pData(palmieri_eset)$SAMPLE_ID, sep="_"))

p <- qplot(PC1, PC2, data = dataGG, color = Group,
           main = "PCA plot of the calibrated data", size = I(5),
           asp = 1.0, geom = "text",
           label = labels) + 
  scale_colour_brewer(palette = "Set2") + theme_bw()
ggsave(filename="PCA_rma.pdf", p, device="pdf", width=12, height=10)

pdf("boxplot_rma.pdf", width=10, height=7)
boxplot(palmieri_eset, target = "core",
        main = "Boxplots of log2-intensities for the raw data")
dev.off()
