# clusterprofile for upregulated and downregulated genes.
library(clusterProfiler)
library(readxl)
library(readr)
library(dplyr)
library(org.Hs.eg.db)
library(cowplot)
library(ReactomePA)

pValue <- 0.05
log2Ratio <- 0

## BEC
# setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BEC_20181105/GOEnrichment")

## BAL
# setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BAL_20181105/GOEnrichment")

comparisonFns <- c("results_groupHEA_DAY4-groupHEA_DAY14.tsv",
                   "results_groupMIA_DAY14-groupHEA_DAY14.tsv",
                   "results_groupMIA_DAY4-groupHEA_DAY4.tsv",
                   "results_groupMIA_DAY4-groupMIA_DAY14.tsv",
                   "results_groupMOA_DAY14-groupHEA_DAY14.tsv",
                   "results_groupMOA_DAY4-groupHEA_DAY4.tsv",
                   "results_groupMOA_DAY4-groupMOA_DAY14.tsv",
                   "results_(groupMIA_DAY4-groupMIA_DAY14)-(groupHEA_DAY4-groupHEA_DAY14).tsv",
                   "results_(groupMIA_DAY4-groupMIA_DAY14)-(groupMOA_DAY4-groupMOA_DAY14).tsv",
                   "results_(groupMOA_DAY4-groupMOA_DAY14)-(groupHEA_DAY4-groupHEA_DAY14).tsv")
for(i in 1:length(comparisonFns)){
message(i)
# foreach(i=1:length(comparisonFns)) %dopar%{
  resultsComp <- read_tsv(file.path("..", comparisonFns[i]))
  resultsComp$ENTREZID <- mapIds(org.Hs.eg.db, keys=resultsComp$SYMBOL,
                                 column="ENTREZID", keytype = "SYMBOL")
  ## Up-regulated 
  resultsCompSig <- filter(resultsComp, `P.Value` <= pValue, `logFC` > log2Ratio)
  ego <- enrichGO(gene          = resultsCompSig$ENTREZID,
                  universe      = resultsComp$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable      = TRUE)
  pUp <- dotplot(ego, title="GOBP up-reglated genes")
  ego <- enrichKEGG(gene          = resultsCompSig$ENTREZID,
                    universe      = resultsComp$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
  pKEGGUp <- dotplot(ego, title="KEGG up-reglated genes")
  ego <- enrichPathway(gene          = resultsCompSig$ENTREZID,
                       universe      = resultsComp$ENTREZID,
                       organism = "human",
                       pvalueCutoff = 0.05)
  pReactomeUp <- dotplot(ego, title="Reactome up-reglated genes")
  ## Down-regulated
  resultsCompSig <- filter(resultsComp, `P.Value` <= pValue, `logFC` < -log2Ratio)
  ego <- enrichGO(gene          = resultsCompSig$ENTREZID,
                  universe      = resultsComp$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable      = TRUE)
  pDown <- dotplot(ego, title="GOBP down-reglated genes")
  ego <- enrichKEGG(gene          = resultsCompSig$ENTREZID,
                    universe      = resultsComp$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
  pKEGGDown <- dotplot(ego, title="KEGG down-reglated genes")
  ego <- enrichPathway(gene          = resultsCompSig$ENTREZID,
                       universe      = resultsComp$ENTREZID,
                       organism = "human",
                       pvalueCutoff = 0.05)
  pReactomeDown <- dotplot(ego, title="Reactome down-reglated genes")
  ## All differential expressed genes
  resultsCompSig <- filter(resultsComp, `P.Value` <= pValue, abs(`logFC`) >= log2Ratio)
  ego <- enrichGO(gene          = resultsCompSig$ENTREZID,
                  universe      = resultsComp$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable      = TRUE)
  pAll <- dotplot(ego, title="GOBP all significant genes")
  ego <- enrichKEGG(gene          = resultsCompSig$ENTREZID,
                    universe      = resultsComp$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
  pKEGGAll <- dotplot(ego, title="KEGG all significant genes")
  ego <- enrichPathway(gene          = resultsCompSig$ENTREZID,
                       universe      = resultsComp$ENTREZID,
                       organism = "human",
                       pvalueCutoff = 0.05)
  pReactomeAll <- dotplot(ego, title="Reactome all significant genes")
  
  p <- plot_grid(pUp, pDown, pAll, labels="AUTO", nrow=3, align="v")
  save_plot(filename=sub("\\.tsv$", "_GOBP.pdf", comparisonFns[i]), plot=p,
            ncol=1, nrow=3, base_height=4, base_aspect_ratio=2.5)
  
  p <- plot_grid(pKEGGUp, pReactomeDown, pKEGGAll, labels="AUTO", nrow=3, align="v")
  save_plot(filename=sub("\\.tsv$", "_KEGG.pdf", comparisonFns[i]), plot=p,
            ncol=1, nrow=3, base_height=4, base_aspect_ratio=2.5)
  
  p <- plot_grid(pReactomeUp, pKEGGDown, pReactomeAll, labels="AUTO", nrow=3, align="v")
  save_plot(filename=sub("\\.tsv$", "_Reactome.pdf", comparisonFns[i]), plot=p,
            ncol=1, nrow=3, base_height=4, base_aspect_ratio=2.5)
}
