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
waveforms <- paste0('peak', 2:12)

# use each waveform in turn as control
output <- lapply(X = 1:length(waveforms), FUN = function(contWaveform) {
  contIdx <- (500 * (contWaveform - 1) + 1):(500 * contWaveform)
  pairwiseComp <- replicate(n = length(waveforms), syntheticData[[1]][contIdx, ], simplify = FALSE)
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
  expected <- rep(1, 5500)
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
