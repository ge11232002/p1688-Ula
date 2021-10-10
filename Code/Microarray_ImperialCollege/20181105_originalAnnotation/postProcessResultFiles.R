# Post process the result files for each gene group
library(readxl)
library(TFBSTools)
library(JASPAR2018)
library(readr)
library(dplyr)

# Load annotation
## Geneappend3.xlsx
geneGroups <- read_excel("/home/gtan/analysis/p1688-Ula/data/Geneappend3 from ImmPort_downloaded on 2018-02-27.xlsx")
geneGroups <- split(geneGroups$Symbol, geneGroups$Category)

## TF
opts <- list()
opts[["species"]] <- 9606
opts[["all_versions"]] <- FALSE
pfmList <- getMatrixSet(JASPAR2018, opts)
geneGroups$TF <- unname(name(pfmList))

## KEGG
keggFns <- list.files(path="/home/gtan/analysis/p1688-Ula/data",
                      pattern="^KEGG.*xlsx$", full.names = TRUE)
for(keggFn in keggFns){
  kegg <- read_excel(keggFn)
  geneGroups[[colnames(kegg)]] <- unname(grep("^>", unlist(kegg), value=TRUE, invert = TRUE))
}

## REACTOME_TIGHT_JUNCTION_INTERACTIONS.xlsx
reactomeFns <- list.files(path="/home/gtan/analysis/p1688-Ula/data",
                          pattern="^REACTOME.*xlsx$", full.names=TRUE)
for(reactomeFn in reactomeFns){
  reactome <- read_excel(reactomeFn)
  geneGroups[[colnames(reactome)]] <- unname(grep("^>", unlist(reactome), value=TRUE, invert = TRUE))
}

## Lipid
lipidsGenes <- read_excel("/home/gtan/analysis/p1688-Ula/data/HUGO GENE_LIPID LIKE RECEPTORS_ALL.xlsx")
geneGroups[["lipidGene"]] <- pull(lipidsGenes)
  
## BEC
# setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BEC_20181105/resultsGenesets")

## BAL
setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BAL_20181105/resultsGenesets")

exprs <- read_tsv("../exprs_normalised_filtered.tsv")
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
  resultsComp <- read_tsv(file.path("..", comparisonFns[i]))
  ans <- left_join(resultsComp, exprs)
  for(group in names(geneGroups)){
    write_tsv(ans[toupper(ans$SYMBOL) %in% toupper(geneGroups[[group]]), ],
              path=gsub(" ", "", paste0("results_", comparisonFns[i], "_", group, ".txt")))
  }
}
