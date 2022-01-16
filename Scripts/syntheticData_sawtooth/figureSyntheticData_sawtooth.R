library(pROC)
library(tidyverse)


nLevel <- c(0.1, 0.2, 0.4)
reps <- c(2, 3)
gridOfparams <- expand.grid(nLevel, reps)

apply(gridOfparams, 1, function(params) {

  # get and load file
  fileName <- paste0("results_syntheticData_sawtooth_2_48_rep_", params[2], "_nl_", params[1])
  load(paste0("~/Desktop/TimeChange-data/Results/syntheticData_sawtooth/timeChange/", fileName, ".RData"))


  perCorrectClassification <- lapply(1:11, function(i) {
    selectOutput <- output[[i]]
    selectOutput$results$p.value <- p.adjust(selectOutput$results$p.value, method = "fdr")
    splitByWaveform <- split(unlist(selectOutput$results$p.value), ceiling(seq_along(unlist(selectOutput$results$p.value)) / 500))
    youdenThreshold <- unname(unlist(coords(selectOutput$ROC, "best", ret = "threshold", transpose = F, best.method = "youden")))
    results <- lapply(splitByWaveform, function(j) {
      sum(j < youdenThreshold) / 500
    })
  })


  x <- rep(paste0("Peak ", seq(2,22,2), "h"), 11)
  y <- c(rep("Peak 2h", 11), rep("Peak 4h", 11), rep("Peak 6h", 11), rep("Peak 8h", 11), rep("Peak 10h", 11), rep("Peak 12h", 11), rep("Peak 14h", 11), rep("Peak 16h", 11), rep("Peak 18h", 11), rep("Peak 20h", 11), rep("Peak 22h", 11))
  z <- unname(unlist(perCorrectClassification))

  plotData <- data.frame(x, y, z)

  ouputFile <- paste0("~/Desktop/TimeChange-data/Results/syntheticData_sawtooth/timeChange/figs/", fileName, ".png")
  png(filename = ouputFile, width = 4, height = 4, units = "in", res = 300)
  p1 <- plotData %>%
    mutate(y = fct_relevel(
      y,
      rev(paste0("Peak ", seq(2,22,2), "h"))
    )) %>%
    mutate(x = fct_relevel(
      x,
      paste0("Peak ", seq(2,22,2), "h")
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
