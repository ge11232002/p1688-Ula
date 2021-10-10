# overlaps of up and down regulated genes between asthma, healthy from d-14 to d4
library(GOplot)
library(readr)
library(dplyr)
library(cowplot)

pValue <- 0.05
log2Ratio <- 0

## BEC
# setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BEC_20181105/VenDiagram")

## BAL
setwd("/home/gtan/analysis/p1688-Ula/ImperialCollege/BAL_20181105/VenDiagram")


comparisonFns <- c(HEA="results_groupHEA_DAY4-groupHEA_DAY14.tsv",
                   MIA="results_groupMIA_DAY4-groupMIA_DAY14.tsv",
                   MOA="results_groupMOA_DAY4-groupMOA_DAY14.tsv")
ans <- list()
i <- 1
for(i in 1:length(comparisonFns)){
  results <- read_tsv(file.path("..", comparisonFns[i]))
  results <- filter(results, `P.Value` <= pValue) %>%
    dplyr::select(SYMBOL, logFC) %>%
    dplyr::rename(ID=SYMBOL) %>%
    dplyr::mutate(ID=make.unique(ID))
  results <- as.data.frame(results)
  ans[[names(comparisonFns)[i]]] <- results
}

p <- GOVenn(ans[[1]], ans[[2]], ans[[3]], label=names(ans),
            title="Venn diagram of differentially expressed genes \nbefore and after infection",
            plot=FALSE)
save_plot(filename="VennDiagram_HEAMIAMOA_D4D14.pdf", plot=p$plot,
          base_height=7)
for(eachTable in names(p$table)){
  resTable <- as_tibble(p$table[[eachTable]], rownames="GeneName")
  write_tsv(resTable, path=paste0("Overlaps_", eachTable, ".tsv"))
}
