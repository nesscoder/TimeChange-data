#get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
print(args)

#load functions and required packages for TimeChange
library(tidyverse)
library(pROC)
library(fasano.franceschini.test)
library(MASS)
source("/home/emn6548/TimeChange-data/Scripts/TimeChange/TimeChange.R")
source("/home/emn6548/TimeChange-data/Scripts/TimeChange/internal_functions.R")

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


  ## --- Run TimeChange
  start <- Sys.time()
  results <- TimeChange(
    sample1 = contExpression,
    sample2 = mutExpression,
    repLabelS1 = rep(reps, ncol(contExpression) / reps)
  )
  end <- Sys.time()
  totalTime <- end - start

  ## --- ROC & AUC Check
  expected <- rep(1, 5500)
  expected[contIdx] <- 0
  pred <- as.numeric(as.vector(unlist(results$p.value)))
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
