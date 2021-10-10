# Bioconductor arrays workflow
# Comparing ICS allergic to HDM vs ICS not allergic to HDM

# Bioconductor pipeline
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
library(ggplot2)
library(oligo)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(readr)
library(hugene10sttranscriptcluster.db)
library(readxl)
library(cowplot)
library(matrixStats)
library(tibble)

# Load inflammasome genes
infGenes <- read_xlsx("/home/gtan/analysis/p1688-Ula/inflammasomeGenes/data/Inflammasomes_and_IL1signaling_gene_list_v2_280217.xlsx")
infGenes <- infGenes$`Gene Symbol`

# grouping of samples
grouping <- list("ICSHDM"=c("DJ113","DJ75","DJ61","DJ136","DJ86"),
                 "ICSHDMNeg"=c("DJ60","DJ125","DJ73","DJ134","DJ126"))
grouping <- stack(grouping) %>% rename(SAMPLE_ID=values,
                                       GroupUla=ind)

samples <- c("BAL", "BEC")
dataRoot <- "/home/gtan/analysis/p1688-Ula/ImperialCollege"
for(sample in samples){
  message(sample)
  wdDir <- file.path(dataRoot, paste0(sample, "_", Sys.Date(), "_ICS"))
  dir.create(wdDir)
  setwd(wdDir)
  phenoData <-
    read.AnnotatedDataFrame(file.path(dataRoot, "data", 
                                      paste0(sample, "_expt_design_fixed.txt")))
  celfiles <- file.path(dataRoot, "data", sample)
  

  sampleAnno <- pData(phenoData)
  sampleAnno <- left_join(sampleAnno, grouping) %>%
    mutate(TIME=sub("DAY_", "DAY", TIME),
           GROUP=GroupUla,
           GT=paste(GROUP, TIME, sep="_"))
  
  pData(phenoData) <- sampleAnno
  phenoData <- phenoData[!is.na(pData(phenoData)$GROUP), ]
  
  eset <- justRMA(phenoData=phenoData,
                  celfile.path=celfiles,
                  filenames=paste0(pData(phenoData)$EXPT, ".CEL"))
  
  ## Output
  exprs <- as_tibble(exprs(eset), rownames="PROBEID")
  exprs <- rename_at(exprs, vars(paste0(pData(phenoData)$EXPT, ".CEL")), ~paste(pData(phenoData)$GT, 
                                                                                pData(phenoData)$SAMPLE_ID, sep="-"))
  write_tsv(exprs, path="exprs_normalised.tsv")
  
  ## QC
  ### hierarchical clusters
  exprsClustering <- exprs
  exprsClustering <- data.frame(exprsClustering, stringsAsFactors = FALSE,
                                check.names = FALSE)
  rownames(exprsClustering) <- exprsClustering$PROBEID
  exprsClustering$PROBEID <- NULL
  exprsClustering <- as.matrix(exprsClustering)
  exprsClustering <- exprsClustering[order(rowVars(exprsClustering), decreasing = TRUE), ]
  exprsClustering <- head(exprsClustering, 2e3)
  annotation_col = data.frame(
    Time = pData(eset)$TIME,
    Group = pData(eset)$GROUP,
    GT = pData(eset)$GT
  )
  rownames(annotation_col) = colnames(exprsClustering)
  exprsClustering <- t(scale(t(exprsClustering), scale=FALSE))
  pheatmap(exprsClustering, color = colorRampPalette(c("blue", "white", "red"))(100),
           main="Clustering of 2000 most variable genes",
           breaks=seq(from=-4, to=4, length.out=101),
           show_rownames=FALSE, annotation_col=annotation_col,
           cluster_rows=FALSE,
           filename="ClusteringOfMostVariableGenes.pdf",
           width=10, height=10)
  ### PCA with all probes
  expValue <- exprs(eset)
  allGroups <- list("ICSHDM"="ICSHDM", "ICSHDMNeg"="ICSHDMNeg", 
                    "All"=c("ICSHDM", "ICSHDMNeg"))
  onegroup <- names(allGroups)[4]
  for(onegroup in names(allGroups)){
    samplesToUse <- pData(eset)$GROUP %in% allGroups[[onegroup]]
    PCA <- prcomp(t(expValue[ ,samplesToUse]),
                  scale = FALSE)
    dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                         Time = pData(eset)$TIME[samplesToUse],
                         Group = pData(eset)$GROUP[samplesToUse],
                         SAMPLE_ID = pData(eset)$SAMPLE_ID[samplesToUse],
                         GT = pData(eset)$GT[samplesToUse],
                         labels = paste(pData(eset)$TIME[samplesToUse], 
                                        pData(eset)$SAMPLE_ID[samplesToUse],
                                        sep="_"))
    p <- qplot(PC1, PC2, data = dataGG, color = Time,
               main = paste("PCA plot of the calibrated data:", onegroup),
               size = I(5),
               asp = 1.0, geom = "text",
               label = labels) + 
      scale_colour_brewer(palette = "Set2") + theme_bw()
    save_plot(filename=paste0("PCA_", onegroup, ".pdf"), p,
              base_height=10, base_aspect_ratio=1.2)
  }
  p <- qplot(PC1, PC2, data = dataGG, color = Group,
             main = paste("PCA plot of the calibrated data:", onegroup),
             size = I(5),
             asp = 1.0, geom = "text",
             label = labels) + 
    scale_colour_brewer(palette = "Set2") + theme_bw()
  save_plot(filename=paste0("PCA.pdf"), p,
            base_height=10, base_aspect_ratio=1.2)
  ### PCA with most variable probes
  
  pdf("boxplot.pdf", width=10, height=7)
  boxplot(eset, target = "core",
          main = "Boxplots of log2-intensities for the normalised data")
  dev.off()
  
  dists <- as.matrix(dist(t(expValue), method = "manhattan"))
  diag(dists) <- NA
  rownames(dists) <- colnames(dists) <- paste(pData(eset)$GT,
                                              pData(eset)$SAMPLE_ID, sep="_")
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
    dplyr::summarize(no_of_matches = n_distinct(SYMBOL)) %>%
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
  exprs <- as_tibble(exprs(palmieri_final), rownames="PROBEID")
  exprs <- rename_at(exprs, vars(paste0(pData(phenoData)$EXPT, ".CEL")), ~paste(pData(phenoData)$GT, 
                                                                                pData(phenoData)$SAMPLE_ID, sep="-"))
  exprsAnno <- merge(fData(palmieri_final), exprs)
  write_tsv(exprsAnno, path="exprs_normalised_filtered.tsv")
  write_tsv(exprsAnno[toupper(exprsAnno$SYMBOL) %in% toupper(infGenes), ],
            path="exprs_normalised_filtered_inflammasomeGenes.tsv")
  
  ## Limma
  sampleID <- factor(pData(phenoData)$SAMPLE_ID)
  group <- factor(pData(phenoData)$GT)
  design <- model.matrix(~0+group+sampleID)
  
  fit <- lmFit(palmieri_final, design)  # fit each probeset to model
  contrasts <- makeContrasts(groupICSHDM_DAY4-groupICSHDM_DAY14,
                             groupICSHDMNeg_DAY4-groupICSHDMNeg_DAY14,
                             
                             groupICSHDM_DAY14-groupICSHDMNeg_DAY14,
                             groupICSHDM_DAY4-groupICSHDMNeg_DAY4,
                             
                             (groupICSHDM_DAY4-groupICSHDM_DAY14)-(groupICSHDMNeg_DAY4-groupICSHDMNeg_DAY14),
                             levels=design)
  fit2 <- contrasts.fit(fit, contrasts)
  efit <- eBayes(fit2)        # empirical Bayes adjustment
  for(i in 1:ncol(contrasts)){
    ans <- topTable(efit, num=Inf, coef=i)
    ansExpr <- merge(ans, exprs)
    write_tsv(ansExpr, path=gsub(" ", "",
                             paste0("results_", colnames(contrasts)[i], ".tsv")))
    write_tsv(ansExpr[toupper(ansExpr$SYMBOL) %in% toupper(infGenes), ],
              path=gsub(" ", "", paste0("results_", colnames(contrasts)[i], "_inflammasomeGenes.tsv")))
  }

}
