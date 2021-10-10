# Not used.
# Bioconductor arrays workflow
# https://www.bioconductor.org/packages/devel/workflows/vignettes/arrays/inst/doc/arrays.html

## Check samples
# library(readxl)
# anno <- read_xlsx("/home/gtan/analysis/p1688-Ula/ImperialCollege/data/David's BAL samples for Switzerland.xlsx")
# setdiff(anno$`Patient ID`, pData(phenoData)$SAMPLE_ID)

setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BAL")

# Bioconductor pipeline
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
library(ggplot2)
library(oligo)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(readr)

## Load data:
phenoData <- 
  read.AnnotatedDataFrame("/home/gtan/analysis/p1688-Ula/ImperialCollege/data/BAL_expt_design_fixed.txt")
celfiles <- "/home/gtan/analysis/p1688-Ula/ImperialCollege/data/BAL"
eset <- justRMA(phenoData=phenoData, celfile.path=celfiles)

## Output
write_tsv(as_tibble(exprs(eset), rownames="PROBEID"), 
          path="exprs_normalised.tsv")

## QC
expValue <- exprs(eset)
PCA <- prcomp(t(expValue), scale = FALSE)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Time = pData(eset)$TIME,
                     Group = pData(eset)$GROUP,
                     SAMPLE_ID = pData(eset)$SAMPLE_ID,
                     GT = pData(eset)$GT,
                     labels = paste(pData(eset)$TIME, 
                                    pData(eset)$SAMPLE_ID, sep="_"))

p <- qplot(PC1, PC2, data = dataGG, color = Group,
       main = "PCA plot of the calibrated data", size = I(5),
       asp = 1.0, geom = "text",
       label = labels) + 
  scale_colour_brewer(palette = "Set2") + theme_bw()
ggsave(filename="PCA.pdf", p, device="pdf", width=12, height=10)

pdf("boxplot.pdf", width=10, height=7)
boxplot(eset, target = "core",
        main = "Boxplots of log2-intensities for the normalised data")
dev.off()

dists <- as.matrix(dist(t(expValue), method = "manhattan"))
diag(dists) <- NA
rownames(dists) <- colnames(dists) <- paste(pData(eset)$GT, pData(eset)$SAMPLE_ID, sep="_")
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)

pheatmap(dists, col = rev(hmcol),
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         main="Distance between samples",
         filename="heatmap.pdf",
         height=10, width=10)

palmieri_filtered <- eset

# Add annotation
anno_palmieri  <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                        keys=(featureNames(palmieri_filtered)),
                                        columns = c("SYMBOL", "GENENAME"),
                                        keytype="PROBEID")
sum(!is.na(anno_palmieri$SYMBOL))

probe_stats <- anno_palmieri %>%
  group_by(PROBEID) %>%
  summarize(no_of_matches = n_distinct(SYMBOL)) %>%
  filter(no_of_matches > 1)

probe_stats

### Exclude probes with multiple mapping or without gene name
ids_to_exlude <- ((featureNames(palmieri_filtered) %in% probe_stats$PROBEID) |
                    featureNames(palmieri_filtered)  %in% subset(anno_palmieri ,
                                                                 is.na(SYMBOL))$PROBEID)
table(ids_to_exlude)

palmieri_final <- subset(palmieri_filtered, !ids_to_exlude)
validObject(palmieri_final)

fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))
fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID
validObject(palmieri_final)

## Export normalised signal
write_tsv(as_tibble(exprs(palmieri_final), rownames="PROBEID"), 
          path="exprs_normalised_filtered.tsv")

## Limma
sampleID <- factor(pData(phenoData)$SAMPLE_ID)
group <- factor(pData(phenoData)$GT)
design <- model.matrix(~0+group+sampleID)

fit <- lmFit(palmieri_final, design)  # fit each probeset to model
contrasts <- makeContrasts(groupHEA_DAY_4 - groupHEA_DAY_14, 
                           groupMIA_DAY_4 - groupMIA_DAY_14, 
                           groupMOA_DAY_4 - groupMOA_DAY_14, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
efit <- eBayes(fit2)        # empirical Bayes adjustment
for(i in 1:ncol(contrasts)){
  ans <- topTable(efit, num=Inf, coef=i)
  write_tsv(ans, path=gsub(" ", "", 
                           paste0("results_", colnames(contrasts)[i], ".tsv")))
}
