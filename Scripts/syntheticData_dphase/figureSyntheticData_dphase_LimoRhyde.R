library(pROC)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)

nLevel <- c(0.1, 0.2, 0.4)
reps <- c(2, 3)
gridOfparams <- expand.grid(nLevel, reps)

apply(gridOfparams, 1, function(params){

  #get and load file
  fileName <- paste0("results_syntheticData_dphase_2_48_rep_", params[2], "_nl_", params[1])
  load(paste0("~/Desktop/TimeChange-data/Results/syntheticData_dphase/limoRhyde/", fileName, ".RData"))

  #calculate percentage of correct classification
  perCorrectClassification <- lapply(1:24, function(i) {
    selectOutput <- output[[i]]
    selectOutput$results$FDR <- p.adjust(selectOutput$results$P.Value,method = 'fdr')
    splitByWaveform <- split(unlist(selectOutput$results$FDR), ceiling(seq_along(unlist(selectOutput$results$FDR)) / 500))
    youdenThreshold <- unname(unlist(coords(selectOutput$ROC, "best", ret = "threshold", transpose = F, best.method = "youden")))
    results <- lapply(splitByWaveform, function(j) {
      sum(j < youdenThreshold[1]) / 500
    })
  })

  x <- rep(paste0("Phase ", 0:23, "h"), 24)
  y <- c(
    rep("Phase 0h", 24), rep("Phase 1h", 24), rep("Phase 2h", 24), rep("Phase 3h", 24),
    rep("Phase 4h", 24), rep("Phase 5h", 24), rep("Phase 6h", 24), rep("Phase 7h", 24),
    rep("Phase 8h", 24), rep("Phase 9h", 24), rep("Phase 10h", 24), rep("Phase 11h", 24),
    rep("Phase 12h", 24), rep("Phase 13h", 24), rep("Phase 14h", 24), rep("Phase 15h", 24),
    rep("Phase 16h", 24), rep("Phase 17h", 24), rep("Phase 18h", 24), rep("Phase 19h", 24),
    rep("Phase 20h", 24), rep("Phase 21h", 24), rep("Phase 22h", 24), rep("Phase 23h", 24)
  )
  z <- unname(unlist(perCorrectClassification))

  plotData <- data.frame(x, y, z)

  ouputFile <- paste0("~/Desktop/TimeChange-data/Results/syntheticData_dphase/limoRhyde/figs/", fileName, ".png")
  png(filename = ouputFile, width = 6, height = 6, units = "in", res = 300)
  p1 <- plotData %>%
    mutate(y = fct_relevel(
      y,
      rev(paste0("Phase ", 0:23, "h"))
    )) %>%
    mutate(x = fct_relevel(
      x,
      paste0("Phase ", 0:23, "h")
    )) %>%
    ggplot(aes(x, y)) +
    geom_tile(aes(fill = z)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust=1, face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )

  print(p1)
  dev.off()
})

