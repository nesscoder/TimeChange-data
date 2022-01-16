library(pROC)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)

nLevel <- c(0.1, 0.2, 0.4)
reps <- c(2, 3)
gridOfparams <- expand.grid(nLevel, reps)

apply(gridOfparams, 1, function(params) {

  #get and load file
  fileName <- paste0("results_syntheticData_2_48_rep_", params[2], "_nl_", params[1])
  load(paste0("~/Desktop/TimeChange-data/Results/syntheticData_general/timeChange/", fileName, ".RData"))

  #calculate percentage of correct classification
  perCorrectClassification <- lapply(1:6, function(i) {
    selectOutput <- output[[i]]
    selectOutput$results$p.value <- p.adjust(selectOutput$results$p.value, method = "fdr")
    splitByWaveform <- split(unlist(selectOutput$results$p.value), ceiling(seq_along(unlist(selectOutput$results$p.value)) / 1000))
    youdenThreshold <- unname(unlist(coords(selectOutput$ROC, "best", ret = "threshold", transpose = F, best.method = "youden")))
    results <- lapply(splitByWaveform, function(j) {
      sum(j < youdenThreshold[1]) / 1000
    })
  })

  x <- rep(c("Low-Amp Sin", "High-Amp Sin", "Sawtooth", "Linear Trend", "Damped", "Contractile"), 6)
  y <- c(rep("Low-Amp Sin", 6), rep("High-Amp Sin", 6), rep("Sawtooth", 6), rep("Linear Trend", 6), rep("Damped", 6), rep("Contractile", 6))
  z <- unname(unlist(perCorrectClassification))
  plotData <- data.frame(x, y, z)

  ouputFile <- paste0("~/Desktop/TimeChange-data/Results/syntheticData_general/timeChange/figs/", fileName, ".png")
  png(filename = ouputFile, width = 4, height = 4, units = "in", res = 300)
  p1 <- plotData %>%
    mutate(y = fct_relevel(
      y,
      "Contractile", "Damped", "Linear Trend", "Sawtooth", "High-Amp Sin", "Low-Amp Sin"
    )) %>%
    mutate(x = fct_relevel(
      x,
      "Low-Amp Sin", "High-Amp Sin", "Sawtooth", "Linear Trend", "Damped", "Contractile"
    )) %>%
    ggplot(aes(x, y)) +
    geom_tile(aes(fill = z)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust=1, face="bold", size = 12),
      axis.text.y = element_text(face="bold", size = 12),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  print(p1)
  dev.off()

})
