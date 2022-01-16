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

contExpression <- syntheticData[[1]]
mutExpression <- syntheticData[[2]]

reps <- parse_number(str_extract(string = file, pattern = "rep_._"))
phaseShift <- paste0('phase', 0:23)

period <- 24
startZT <- 0
len <- 48
interval <- 2

ZTtp <- seq(startZT,len, interval)
baseTPs <- length(ZTtp)
repLabel <- rep(reps, ncol(contExpression) / reps)
idxNum <- table(sort(ZTtp%%period))
tpReps <- split(colnames(contExpression), rep(1:length(repLabel), repLabel))
df <- data.frame(idx = 1:baseTPs, rep = (ZTtp%%period/interval+1)) %>%
  arrange(rep)
tpIdx <- as.vector(sapply(df$idx, FUN = function(i){
  tpReps[[i]]
}))
repLabel <- as.vector(unlist(lapply(split(repLabel, rep(names(idxNum),unname(idxNum))),sum)))
contExpression <- contExpression[,tpIdx]
mutExpression <- mutExpression[,tpIdx]

# use each waveform in turn as control
output <- lapply(X = 1:length(phaseShift), FUN = function(phase) {
  contIdx <- (500 * (phase - 1) + 1):(500 * phase)
  pairwiseComp <- replicate(n = length(phaseShift), contExpression[contIdx, ], simplify = FALSE)
  contExpression <- do.call(rbind.data.frame, pairwiseComp)
  rownames(contExpression) <- paste0(phaseShift[phase], "_", 1:nrow(contExpression))


  ## --- Run TimeChange
  start <- Sys.time()
  results <- TimeChange(
    sample1 = contExpression,
    sample2 = mutExpression,
    repLabelS1 = repLabel
  )
  end <- Sys.time()
  totalTime <- end - start

  ## --- ROC & AUC Check
  expected <- rep(1, 12000)
  expected[contIdx] <- 0
  pred <- as.numeric(as.vector(unlist(results$p.value)))
  ROC <- roc(expected, pred)

  fileName <- paste0(phaseShift[phase], "results_", file)
  return(list(fileName = fileName,
              inputData = list(control = contExpression, mutant = mutExpression),
              results = results,
              totalTime = totalTime,
              ROC = ROC))
})

names(output) <- phaseShift
rm(list = setdiff(ls(), c("args", "output", "syntheticData")))

save.image(file = args[3])
