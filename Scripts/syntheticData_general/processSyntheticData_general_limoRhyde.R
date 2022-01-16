#get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
print(args)

library('annotate')
library('tidyverse')
library('data.table')
library('foreach')
library('ggplot2')
library('knitr')
library('limma')
library('limorhyde')
library('readr')
library('stringr')
library('pROC')


source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))

# load Data
file <- args[1]
load(args[2])

mutExpression <- syntheticData[[2]]
reps <- parse_number(str_extract(string = file, pattern = "rep_._"))
waveforms <- c("sin", "sin2X", "saw", "LT", "Damped", "Contract")


output <- lapply(X = 1:6, FUN = function(contWaveform) {
  contIdx <- (1000 * (contWaveform - 1) + 1):(1000 * contWaveform)
  pairwiseComp <- replicate(n = 6, syntheticData[[1]][contIdx, ], simplify = FALSE)
  contExpression <- do.call(rbind.data.frame, pairwiseComp)
  rownames(contExpression) <- paste0(waveforms[contWaveform], "_", 1:nrow(contExpression))

  ## --- Run LimoRhyde
  start <- Sys.time()
  emat <- cbind(contExpression, mutExpression)
  timeZT <- parse_number(colnames(emat))
  sm <- data.frame(cond = c(rep("wt", 25*reps), rep("mut", 25*reps)), time = timeZT, limorhyde(timeZT, 'time_'))
  design <- model.matrix(~ cond * (time_cos + time_sin), data = sm)
  fit <- lmFit(emat, design)
  fit <- eBayes(fit, trend = TRUE)
  results <- data.table(topTable(fit, coef = 5:6, number = Inf,sort.by = "none"), keep.rownames = TRUE)
  end <- Sys.time()
  totalTime <- end - start

  ## --- ROC & AUC Check
  expected <- rep(1, 6000)
  expected[contIdx] <- 0
  pred <- as.numeric(as.vector(unlist(results$P.Value)))
  ROC <- roc(expected, pred)

  fileName <- paste0(waveforms[contWaveform], "results_", file)
  return(list(fileName = fileName,
              inputData = list(control = contExpression, mutant = mutExpression),
              results = results,
              totalTime = totalTime,
              ROC = ROC))
})

names(output) <- waveforms
rm(list = setdiff(ls(), c("args", "output", "syntheticData")))

save.image(file = args[3])
