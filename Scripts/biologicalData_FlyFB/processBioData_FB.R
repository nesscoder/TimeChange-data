library(devtools)
library(tidyverse)
library(matrixStats)
library(readr)
library(MASS)
library(circacompare)
library(clusterProfiler)
library(ggnewscale)
library(enrichplot)
library(GOSemSim)
library(DOSE)
library(org.Dm.eg.db)
library(RcisTarget)
library(DT)
library(ggpubr)
library(viridis)
library(ggupset)
library(RColorBrewer)
library(europepmc)
library(fasano.franceschini.test)
source("../TimeChange/TimeChange.R")
source("../TimeChange/internal_functions.R")

# BiocManager::install("org.Dm.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("RcisTarget")
# BiocManager::install("DOSE")

# Functions Definitions ------------
source("~/Desktop/TimeChange-data/Scripts/biologicalData_FlyFB/JTK_CYCLEv3.1.R")

loadData <- function(dataFile) {
  data <- read.csv(dataFile, sep = "\t")
  rownames(data) <- data[, 1]
  data <- data[, -1]
  orderByRep <- colnames(data)[order(as.numeric(gsub(pattern = ".*ZT",replacement = "", colnames(data))))]
  data <- data[, orderByRep]
  return(data)
}

runJTK <- function(data, timepoints, reps, start, end, interval) {
  annot <- rownames(data)
  jtkdist(as.numeric(timepoints), as.numeric(reps))

  periods <- as.numeric(start):as.numeric(end)
  jtk.init(periods, as.numeric(interval))

  cat("JTK analysis started on", date(), "\n")
  flush.console()

  st <- system.time({
    res <- apply(data, 1, function(z) {
      jtkx(z)
      c(JTK.ADJP, JTK.PERIOD, JTK.LAG, JTK.AMP)
    })
    res <- as.data.frame(t(res))
    bhq <- p.adjust(unlist(res[, 1]), "BH")
    res <- cbind(bhq, res)
    colnames(res) <- c("BH.Q", "ADJ.P", "PER", "LAG", "AMP")
    results <- cbind(annot, res)
    # results <- results[order(res$ADJ.P,-res$AMP),]
  })
  print(st)
  return(results)
}

# Process HC18 vs HC25 with JTK Cycle ----------------------------------------------------
HC18_data <- loadData("~/Desktop/TimeChange-data/Data/biologicalData_FlyFB/FB_HC18_Iso_abundance.tsv")
HC25_data <- loadData("~/Desktop/TimeChange-data/Data/biologicalData_FlyFB/FB_HC25_Iso_abundance.tsv")
HC25_JRK_data <- loadData("~/Desktop/TimeChange-data/Data/biologicalData_FlyFB/FB_HC25_Jrk_abundance.tsv")

# set seed for reproducibility
set.seed(2020)

## pre-filter using row median
HC18_data <- HC18_data[rowMedians(as.matrix(HC18_data)) >= 2, ] ## 7884
HC25_data <- HC25_data[rowMedians(as.matrix(HC25_data)) >= 2, ] ## 7885
HC25_JRK_data <- HC25_JRK_data[rowMedians(as.matrix(HC25_JRK_data)) >= 2, ] ## 7458
cmn_names <- Reduce(intersect, list(rownames(HC18_data), rownames(HC25_data), rownames(HC25_JRK_data))) ## 7333
HC18_data <- HC18_data[cmn_names, ]
HC25_data <- HC25_data[cmn_names, ]
HC25_JRK_data <- HC25_JRK_data[cmn_names, ]

# Run TimeChange on pairwise interaction of data sets
start <- Sys.time()
results_HC18vsHC25 <- TimeChange(
  sample1 = HC18_data,
  sample2 = HC25_data,
  repLabelS1 = rep(2, ncol(HC18_data) / 2),
  repLabelS2 = rep(2, ncol(HC25_data) / 2)
)
end <- Sys.time()
print(end - start)

start <- Sys.time()
results_HC18vsHC25_JRK <- TimeChange(
  sample1 = HC18_data,
  sample2 = HC25_JRK_data,
  repLabelS1 = rep(2, ncol(HC18_data) / 2),
  repLabelS2 = rep(2, ncol(HC25_JRK_data) / 2),
  maxLag = 3
)
end <- Sys.time()
print(end - start)

start <- Sys.time()
results_HC25vsHC25_JRK <- TimeChange(
  sample1 = HC25_data,
  sample2 = HC25_JRK_data,
  repLabelS1 = rep(2, ncol(HC18_data) / 2),
  repLabelS2 = rep(2, ncol(HC25_JRK_data) / 2),
  maxLag = 3
)
end <- Sys.time()
print(end - start)

# FDR correct for multiple hypothesis testing
results_HC18vsHC25$FDR <- p.adjust(results_HC18vsHC25$p.value, method = "fdr")
results_HC18vsHC25_JRK$FDR <- p.adjust(results_HC18vsHC25_JRK$p.value, method = "fdr")
results_HC25vsHC25_JRK$FDR <- p.adjust(results_HC25vsHC25_JRK$p.value, method = "fdr")

# Get significant deferentially dynamic genes ----------------------------------------------
hist(as.numeric(results_HC18vsHC25$FDR), breaks = seq(0.0, 1, 0.025)) # Bimodal distribution. set break at 0.1
threshold <- 0.1
HC18vs25_sigGenes <- unlist(results_HC18vsHC25[results_HC18vsHC25$FDR < threshold, 1]) # 997
HC18vs25_nonSigGenes <- unlist(results_HC18vsHC25[results_HC18vsHC25$FDR > threshold, 1]) # 6336
HC18vsJRK_sigGenes <- unlist(results_HC18vsHC25_JRK[results_HC18vsHC25_JRK$FDR < threshold, 1]) # 2449
HC25vsJRK_sigGenes <- unlist(results_HC25vsHC25_JRK[results_HC25vsHC25_JRK$FDR < threshold, 1]) # 2181

# Run JTK_CYCLE on HC18 vs HC25 for deferentially Cycling Genes
HC18vs25_JTK_results <- runJTK(
  data = HC18_data[HC18vs25_sigGenes, ],
  timepoints = 12,
  reps = 2,
  start = 12,
  end = 12,
  interval = 2
)
HC25vs18_JTK_results <- runJTK(
  data = HC25_data[HC18vs25_sigGenes, ],
  timepoints = 12,
  reps = 2,
  start = 12,
  end = 12,
  interval = 2
)

# Run JTK_CYCLE on HC18 vs HC25 for non-differentially Cycling Genes
HC18vs25_JTK_nonDiff_results <- runJTK(
  data = HC18_data[HC18vs25_nonSigGenes, ],
  timepoints = 12,
  reps = 2,
  start = 12,
  end = 12,
  interval = 2
)
HC25vs18_JTK_nonDiff_results <- runJTK(
  data = HC25_data[HC18vs25_nonSigGenes, ],
  timepoints = 12,
  reps = 2,
  start = 12,
  end = 12,
  interval = 2
)

# Run JTK_CYCLE on HC18 vs HC25 JRK for Cycling Genes
HC18vsJRK_JTK_results <- runJTK(
  data = HC18_data[HC18vsJRK_sigGenes, ],
  timepoints = 12,
  reps = 2,
  start = 12,
  end = 12,
  interval = 2
)
JRKvsHC18_JTK_results <- runJTK(
  data = HC25_JRK_data[HC18vsJRK_sigGenes, ],
  timepoints = 6,
  reps = 2,
  start = 6,
  end = 6,
  interval = 4
)

# Run JTK_CYCLE on HC25 vs HC25 JRK for Cycling Genes
HC25vsJRK_JTK_results <- runJTK(
  data = HC25_data[HC25vsJRK_sigGenes, ],
  timepoints = 12,
  reps = 2,
  start = 12,
  end = 12,
  interval = 2
)
JRKvsHC25_JTK_results <- runJTK(
  data = HC25_JRK_data[HC25vsJRK_sigGenes, ],
  timepoints = 6,
  reps = 2,
  start = 6,
  end = 6,
  interval = 4
)

cycThreshold <- 0.1
diffCyclingGenes_HC18vsHC25 <- HC18vs25_JTK_results[HC18vs25_JTK_results$BH.Q < cycThreshold, ] # 59
diffCyclingGenes_HC25vsHC18 <- HC25vs18_JTK_results[HC25vs18_JTK_results$BH.Q < cycThreshold, ] # 11
nondiffCyclingGenes_HC18vsHC25 <- HC18vs25_JTK_nonDiff_results[HC18vs25_JTK_nonDiff_results$BH.Q < cycThreshold, ] # 76
nondiffCyclingGenes_HC25vsHC18 <- HC25vs18_JTK_nonDiff_results[HC25vs18_JTK_nonDiff_results$BH.Q < cycThreshold, ] # 23
diffCyclingGenes_HC18vsHCJRK <- HC18vsJRK_JTK_results[HC18vsJRK_JTK_results$BH.Q < cycThreshold, ] # 114
diffCyclingGenes_HC25vsHCJRK <- HC25vsJRK_JTK_results[HC25vsJRK_JTK_results$BH.Q < cycThreshold, ] # 22
# no genes cycling in JRK mutants
sum(JRKvsHC18_JTK_results$BH.Q < cycThreshold) # 0
sum(JRKvsHC25_JTK_results$BH.Q < cycThreshold) # 0

# check overlap between sets
# get significantly differentially cycling genes only present at HC18, ie genes that have no rhythm at 25 and gain rhythms at 18
gainofCycling_25to18 <- setdiff(diffCyclingGenes_HC18vsHC25$annot, diffCyclingGenes_HC25vsHC18$annot) # 54
# get significantly differentially cycling genes only present at HC25, ie genes that gain rhythm at 25 and have no rhythm at 18
lossofCycling_25to18 <- setdiff(diffCyclingGenes_HC25vsHC18$annot, diffCyclingGenes_HC18vsHC25$annot) # 6
# get significantly differentially cycling genes in both conditions
common_cyclers <- intersect(diffCyclingGenes_HC25vsHC18$annot, diffCyclingGenes_HC18vsHC25$annot) # 5
# get significantly cycling genes in both conditions
common_nondiffcyclers <- intersect(nondiffCyclingGenes_HC18vsHC25$annot, nondiffCyclingGenes_HC25vsHC18$annot) # 12

save.image("~/Desktop/TimeChange-data/Results/biologicalData_FlyFB/Results_HC18vsHC25vsJRK_JTK_TimeChange.RData")
